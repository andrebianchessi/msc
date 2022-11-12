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

    Poly lr;
    lr.n = n;
    lr.order = order;

    // create root node
    // the root node will have order+1 children
    Node* root = new Node(order + 1);
    lr.coefficients = root;

    // We now perform BFS, "layer by layer"
    // each layer is a depth in the coefficients tree.
    std::queue<Node*> q;
    q.push(lr.coefficients);
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

    r.val = lr;
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
    r.val = sumOfChildren;
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
    r.val = sumOfChildren;
    return r;
}

// Auxiliary DFS recursive function used on DDxi()
double dfsDDxi(int i, Node* node, int nodeTreeDepth, std::vector<double>* X) {
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
        sumOfChildren += dfsDDxi(i, node->children[c], nodeTreeDepth + 1, X);
    }
    if (derivateThisNode) {
        return node->exp * (node->exp - 1) *
               powOrZero(X->at(nodeTreeDepth - 1), node->exp - 2) *
               sumOfChildren;
    }
    return std::pow(X->at(nodeTreeDepth - 1), node->exp) * sumOfChildren;
}
Maybe<double> Poly::DDxi(int i, std::vector<double>* X) {
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
        sumOfChildren += dfsDDxi(i, root->children[c], 1, X);
    }
    r.val = sumOfChildren;
    return r;
}