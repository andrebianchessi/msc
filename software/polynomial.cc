#include "polynomial.h"

#include <boost/math/special_functions/binomial.hpp>
#include <queue>

#include "maybe.h"
#include "utils.h"

Node::Node(int nChildren) {
    this->children = std::vector<Node*>(nChildren);
    this->childrenI = 0;
    this->exp = 0;
    this->parentsExpSum = 0;
}

Node::Node(int exp, int parentsExpSum, int nChildren) {
    this->exp = exp;
    this->parentsExpSum = parentsExpSum;
    this->children = std::vector<Node*>(nChildren);
    this->childrenI = 0;
}

Node* Node::AddChild(int exp, int parentsExpSum, int nChildren) {
    Node* child = new Node(exp, parentsExpSum, nChildren);
    this->children[this->childrenI] = child;
    this->childrenI += 1;
    return child;
}

int Poly::nTerms() {
    // Source:
    // https://mathoverflow.net/questions/225953/number-of-polynomial-terms-for-certain-degree-and-certain-number-of-variables/225963#225963?newreg=2a0208ceb740461d8eaa21e304b0e341
    return int(boost::math::binomial_coefficient<double>(this->n + this->order,
                                                         this->order));
}

int nChilds(int nodeParentsExpSum, int nodeExp, int polyOrder,
            int nodeTreeDepth, int nVars) {
    if (nodeTreeDepth == nVars) {
        return 0;
    }
    int currentExpSum = nodeParentsExpSum + nodeExp;
    int childMaxExp = polyOrder - currentExpSum;
    int childMinExp = 0;
    return childMaxExp - childMinExp + 1;
}

Maybe<Poly> Poly::NewPoly(int n, int order, double coefficients) {
    Maybe<Poly> r;
    if (n < 0 || order < 0) {
        r.isError = true;
        r.errMsg = "n and order must be >=0";
        return r;
    }

    Poly p;
    p.k = 1;
    p.n = n;
    p.order = order;

    // create root node
    // the root node will have order+1 children
    Node* root = new Node(order + 1);
    p.coefficients = root;

    // We now perform BFS, "layer by layer"
    // each layer is a depth in the coefficients tree.
    std::queue<Node*> q;
    q.push(p.coefficients);
    for (int treeDepth = 0; treeDepth <= n; treeDepth++) {
        for (int qI = q.size(); qI > 0; qI--) {
            Node* currentNode = q.front();
            q.pop();

            if (treeDepth == n) {
                currentNode->a = coefficients;
            } else {
                // Number of children current node will have
                // their exps will go from 0 to `currentNodeChildren-1`
                int currentNodeTreeDepth = treeDepth;
                int currentNodeChildren =
                    nChilds(currentNode->parentsExpSum, currentNode->exp, order,
                            currentNodeTreeDepth, n);

                for (int childExp = currentNodeChildren - 1; childExp >= 0;
                     childExp--) {
                    int childNodeParentsExpSum =
                        currentNode->parentsExpSum + currentNode->exp;
                    int childNodeTreeDepth = treeDepth + 1;
                    int childNodeChildren =
                        nChilds(childNodeParentsExpSum, childExp, order,
                                childNodeTreeDepth, n);
                    auto childRef = currentNode->AddChild(
                        childExp, childNodeParentsExpSum, childNodeChildren);
                    q.push(childRef);
                }
            }
        }
    }

    r.val = p;
    return r;
};

// Auxiliary DFS recursive function used on operator()
double dfs(Node* node, int nodeTreeDepth, std::vector<double>* X) {
    if (node->IsLeaf()) {
        return node->a * std::pow(X->at(nodeTreeDepth - 1), node->exp);
    }
    double sumOfChildren = 0;
    for (int c = 0; c < int(node->children.size()); c++) {
        sumOfChildren += dfs(node->children[c], nodeTreeDepth + 1, X);
    }
    return std::pow(X->at(nodeTreeDepth - 1), node->exp) * sumOfChildren;
}
Maybe<double> Poly::operator()(std::vector<double>* X) {
    Maybe<double> r;
    if (int(X->size()) != this->n) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }

    Node* root = this->coefficients;
    double sumOfChildren = 0;
    for (int c = 0; c < int(root->children.size()); c++) {
        sumOfChildren += dfs(root->children[c], 1, X);
    }
    r.val = this->k * sumOfChildren;
    return r;
};

// Auxiliary function that returns x^pow if pow>=0 and 0 otherwise.
// This makes it easier to write the functions such as dfsDxi
double powOrZero(double x, int pow) {
    if (pow == 0) {
        return 1;
    }
    if (pow > 0) {
        return std::pow(x, pow);
    }
    return 0;
}

// Auxiliary DFS recursive function used on Dxi()
double dfsDxi(int i, Node* node, int nodeTreeDepth, std::vector<double>* X) {
    // indicates if the derivative is with respect to the variable of this node
    bool derivateThisNode = i == (nodeTreeDepth - 1);
    if (node->IsLeaf()) {
        if (derivateThisNode) {
            return node->exp * node->a *
                   powOrZero(X->at(nodeTreeDepth - 1), node->exp - 1);
        }
        return node->a * std::pow(X->at(nodeTreeDepth - 1), node->exp);
    }
    double sumOfChildren = 0;
    for (int c = 0; c < int(node->children.size()); c++) {
        sumOfChildren += dfsDxi(i, node->children[c], nodeTreeDepth + 1, X);
    }
    if (derivateThisNode) {
        return node->exp * powOrZero(X->at(nodeTreeDepth - 1), node->exp - 1) *
               sumOfChildren;
    }
    return std::pow(X->at(nodeTreeDepth - 1), node->exp) * sumOfChildren;
}
Maybe<double> Poly::Dxi(int i, std::vector<double>* X) {
    Maybe<double> r;
    if (int(X->size()) != this->n) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }
    if (i < 0 || i >= this->n) {
        r.isError = true;
        r.errMsg = "invalid i";
        return r;
    }

    Node* root = this->coefficients;
    double sumOfChildren = 0;
    for (int c = 0; c < int(root->children.size()); c++) {
        sumOfChildren += dfsDxi(i, root->children[c], 1, X);
    }
    r.val = this->k * sumOfChildren;
    return r;
}

// Auxiliary DFS recursive function used on D2xi()
double dfsD2xi(int i, Node* node, int nodeTreeDepth, std::vector<double>* X) {
    // indicates if the derivative is with respect to the variable of this node
    bool derivateThisNode = i == (nodeTreeDepth - 1);
    if (node->IsLeaf()) {
        if (derivateThisNode) {
            return node->exp * (node->exp - 1) * node->a *
                   powOrZero(X->at(nodeTreeDepth - 1), node->exp - 2);
        }
        return node->a * std::pow(X->at(nodeTreeDepth - 1), node->exp);
    }
    double sumOfChildren = 0;
    for (int c = 0; c < int(node->children.size()); c++) {
        sumOfChildren += dfsD2xi(i, node->children[c], nodeTreeDepth + 1, X);
    }
    if (derivateThisNode) {
        return node->exp * (node->exp - 1) *
               powOrZero(X->at(nodeTreeDepth - 1), node->exp - 2) *
               sumOfChildren;
    }
    return std::pow(X->at(nodeTreeDepth - 1), node->exp) * sumOfChildren;
}
Maybe<double> Poly::D2xi(int i, std::vector<double>* X) {
    Maybe<double> r;
    if (int(X->size()) != this->n) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }
    if (i < 0 || i >= this->n) {
        r.isError = true;
        r.errMsg = "invalid i";
        return r;
    }

    Node* root = this->coefficients;
    double sumOfChildren = 0;
    for (int c = 0; c < int(root->children.size()); c++) {
        sumOfChildren += dfsD2xi(i, root->children[c], 1, X);
    }
    r.val = this->k * sumOfChildren;
    return r;
}

// Auxiliary DFS recursive function used on Da
void dfsDa(std::vector<double>* X, double parentsProduct, Node* thisNode,
           int thisNodeTreeDepth, int* aIndex, std::vector<double>* target) {
    const double currentProduct =
        parentsProduct * std::pow(X->at(thisNodeTreeDepth - 1), thisNode->exp);
    if (thisNode->IsLeaf()) {
        target->at(*aIndex) = currentProduct;
        *aIndex = *aIndex + 1;
        return;
    }
    for (int c = 0; c < int(thisNode->children.size()); c++) {
        dfsDa(X, currentProduct, thisNode->children[c], thisNodeTreeDepth + 1,
              aIndex, target);
    }
}
Maybe<Void> Poly::Da(std::vector<double>* X, std::vector<double>* target) {
    Maybe<Void> r;
    if (int(target->size()) != this->nTerms()) {
        r.isError = true;
        r.errMsg = "target must have same length as the number of terms";
        return r;
    }
    if (int(X->size()) != this->n) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }

    int aIndex = 0;
    for (int c = 0; c < int(this->coefficients->children.size()); c++) {
        dfsDa(X, 1, this->coefficients->children[c], 1, &aIndex, target);
    }
    return r;
}

// Auxiliary DFS recursive function used on DaDxi
void dfsDaDxi(int i, std::vector<double>* X, double parentsProduct,
              Node* thisNode, int thisNodeTreeDepth, int* aIndex,
              std::vector<double>* target) {
    bool derivateThisNode = i == (thisNodeTreeDepth - 1);
    double currentProduct;
    if (derivateThisNode) {
        currentProduct =
            parentsProduct * thisNode->exp *
            powOrZero(X->at(thisNodeTreeDepth - 1), thisNode->exp - 1);
    } else {
        currentProduct = parentsProduct *
                         std::pow(X->at(thisNodeTreeDepth - 1), thisNode->exp);
    }
    if (thisNode->IsLeaf()) {
        target->at(*aIndex) = currentProduct;
        *aIndex = *aIndex + 1;
        return;
    }
    for (int c = 0; c < int(thisNode->children.size()); c++) {
        dfsDaDxi(i, X, currentProduct, thisNode->children[c],
                 thisNodeTreeDepth + 1, aIndex, target);
    }
}
Maybe<Void> Poly::DaDxi(int i, std::vector<double>* X,
                        std::vector<double>* target) {
    Maybe<Void> r;
    if (int(target->size()) != this->nTerms()) {
        r.isError = true;
        r.errMsg = "target must have same length as the number of terms";
        return r;
    }
    if (int(X->size()) != this->n) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }
    if (i < 0 || i >= this->n) {
        r.isError = true;
        r.errMsg = "invalid i";
        return r;
    }

    int aIndex = 0;
    for (int c = 0; c < int(this->coefficients->children.size()); c++) {
        dfsDaDxi(i, X, 1, this->coefficients->children[c], 1, &aIndex, target);
    }
    return r;
}

// Auxiliary DFS recursive function used on DaDxi
void dfsDaD2xi(int i, std::vector<double>* X, double parentsProduct,
               Node* thisNode, int thisNodeTreeDepth, int* aIndex,
               std::vector<double>* target) {
    bool derivateThisNode = i == (thisNodeTreeDepth - 1);
    double currentProduct;
    if (derivateThisNode) {
        currentProduct =
            parentsProduct * thisNode->exp * (thisNode->exp - 1) *
            powOrZero(X->at(thisNodeTreeDepth - 1), thisNode->exp - 2);
    } else {
        currentProduct = parentsProduct *
                         std::pow(X->at(thisNodeTreeDepth - 1), thisNode->exp);
    }
    if (thisNode->IsLeaf()) {
        target->at(*aIndex) = currentProduct;
        *aIndex = *aIndex + 1;
        return;
    }
    for (int c = 0; c < int(thisNode->children.size()); c++) {
        dfsDaD2xi(i, X, currentProduct, thisNode->children[c],
                  thisNodeTreeDepth + 1, aIndex, target);
    }
}
Maybe<Void> Poly::DaD2xi(int i, std::vector<double>* X,
                         std::vector<double>* target) {
    Maybe<Void> r;
    if (int(target->size()) != this->nTerms()) {
        r.isError = true;
        r.errMsg = "target must have same length as the number of terms";
        return r;
    }
    if (int(X->size()) != this->n) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }
    if (i < 0 || i >= this->n) {
        r.isError = true;
        r.errMsg = "invalid i";
        return r;
    }

    int aIndex = 0;
    for (int c = 0; c < int(this->coefficients->children.size()); c++) {
        dfsDaD2xi(i, X, 1, this->coefficients->children[c], 1, &aIndex, target);
    }
    return r;
}

// Auxiliary DFS recursive function used on GetCoefficients
void dfsGetCoefficients(double polynomialK, Node* thisNode, int* i,
                        std::vector<double>* target) {
    if (thisNode->IsLeaf()) {
        target->at(*i) = polynomialK * thisNode->a;
        *i = *i + 1;
        return;
    }
    for (int c = 0; c < int(thisNode->children.size()); c++) {
        dfsGetCoefficients(polynomialK, thisNode->children[c], i, target);
    }
}
Maybe<Void> Poly::GetCoefficients(std::vector<double>* target) {
    Maybe<Void> r;
    if (int(target->size()) != this->nTerms()) {
        r.isError = true;
        r.errMsg = "target must have same length as the number of terms";
        return r;
    }
    int i = 0;
    for (int c = 0; c < int(this->coefficients->children.size()); c++) {
        dfsGetCoefficients(this->k, this->coefficients->children[c], &i,
                           target);
    }
    return r;
}

void Poly::Multiply(double k) { this->k = k; }