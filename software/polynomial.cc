#include "polynomial.h"

#include <boost/math/special_functions/binomial.hpp>
#include <cassert>
#include <memory>
#include <queue>

#include "maybe.h"
#include "utils.h"

Node::Node(int nChildren) {
    this->children = std::vector<std::shared_ptr<Node>>(nChildren);
    this->childrenI = 0;
    this->exp = 0;
    this->parentsExpSum = 0;
}

Node::Node(int exp, int parentsExpSum, int nChildren) {
    this->exp = exp;
    this->parentsExpSum = parentsExpSum;
    this->children = std::vector<std::shared_ptr<Node>>(nChildren);
    this->childrenI = 0;
}

std::shared_ptr<Node> Node::AddChild(int exp, int parentsExpSum,
                                     int nChildren) {
    std::shared_ptr<Node> child =
        std::make_shared<Node>(exp, parentsExpSum, nChildren);
    this->children[this->childrenI] = child;
    this->childrenI += 1;
    return child;
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

Maybe<Void> Poly::Build(int n, int order, double coefficients) {
    Maybe<Void> r;
    if (n < 0 || order < 0) {
        r.isError = true;
        r.errMsg = "n and order must be >=0";
        return r;
    }

    this->n = n;
    this->order = order;
    // Source:
    // https://mathoverflow.net/questions/225953/number-of-polynomial-terms-for-certain-degree-and-certain-number-of-variables/225963#225963?newreg=2a0208ceb740461d8eaa21e304b0e341
    this->nTerms = int(boost::math::binomial_coefficient<double>(
        this->n + this->order, this->order));

    // create root node
    // the root node will have order+1 children
    this->coefficients = std::make_shared<Node>(order + 1);

    this->leafNodes = std::vector<std::shared_ptr<Node>>(this->nTerms);
    int leafNodesI = 0;

    // We now perform BFS, "layer by layer"
    // each layer is a depth in the coefficients tree.
    std::queue<std::shared_ptr<Node>> q;
    q.push(this->coefficients);
    for (int treeDepth = 0; treeDepth <= n; treeDepth++) {
        for (int qI = q.size(); qI > 0; qI--) {
            std::shared_ptr<Node> currentNode = q.front();
            q.pop();

            if (treeDepth == n) {
                currentNode->a = coefficients;
                this->leafNodes[leafNodesI] = currentNode;
                leafNodesI += 1;
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
    this->isZero = false;
    return r;
};

// Auxiliary DFS recursive function used on operator()
double dfs(std::shared_ptr<Node> node, int nodeTreeDepth,
           std::vector<double>* X) {
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

    std::shared_ptr<Node> root = this->coefficients;
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

// Auxiliary DFS recursive function used on Dxi
double dfsDxi(int i, std::shared_ptr<Node> node, int nodeTreeDepth,
              std::vector<double>* X) {
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

// Auxiliary DFS recursive function used on GetDxi
// This goes through all the nodes and, if the derivative is with respect to
// the variable that node represents, it is differentiated.
// Note that this doesn't remove "empty" nodes from the tree. Just the exponents
// of the nodes and the coefficients (values at leaf nodes) are adjusted. Thus,
// if you call this function multiple times you'll end up with a tree in which
// all exponents (exp fields) are zero, and the coefficients are are zero; which
// works but might be considered inefficient. IMO this is totally fine
// especially if we only take one or two derivatives;
void dfsSetAsDxi(int i, std::shared_ptr<Node> node, int nodeTreeDepth,
                 int multiplyChildLeafsBy) {
    // indicates if the derivative is with respect to the variable of this node
    bool derivateThisNode = i == (nodeTreeDepth - 1);
    if (derivateThisNode) {
        assert(multiplyChildLeafsBy == 1);
        if (node->exp == 0) {
            multiplyChildLeafsBy = 0;
        } else {
            multiplyChildLeafsBy = node->exp;
            node->exp -= 1;
        }
    }
    if (node->IsLeaf()) {
        node->a *= multiplyChildLeafsBy;
        return;
    }

    for (int c = 0; c < int(node->children.size()); c++) {
        dfsSetAsDxi(i, node->children[c], nodeTreeDepth + 1,
                    multiplyChildLeafsBy);
    }
}
void dfsFixDxiTree(std::shared_ptr<Node> node, std::shared_ptr<Node> dxiNode) {
    if (node->IsLeaf()) {
        node->a += dxiNode->a;
        return;
    }

    std::shared_ptr<Node> dxiNodeChild;
    std::shared_ptr<Node> nodeChild;
    for (int dc = 0; dc < int(dxiNode->children.size()); dc++) {
        dxiNodeChild = dxiNode->children[dc];
        for (int c = 0; c < int(node->children.size()); c++) {
            nodeChild = node->children[c];
            if (nodeChild->exp == dxiNodeChild->exp) {
                dfsFixDxiTree(nodeChild, dxiNodeChild);
            }
        }
    }
}
Maybe<Void> Poly::Dxi(int i) {
    Maybe<Void> r;
    if (i < 0 || i >= this->n) {
        r.isError = true;
        r.errMsg = "Invalid i";
        return r;
    }

    Poly dxi = (*this);
    for (int c = 0; c < int(dxi.coefficients->children.size()); c++) {
        dfsSetAsDxi(i, dxi.coefficients->children[c], 1, 1);
    }

    for (int i = 0; i < (this->leafNodes.size()); i++) {
        this->leafNodes[i]->a = 0;
    }
    dfsFixDxiTree(this->coefficients, dxi.coefficients);

    return r;
}

// Auxiliary DFS recursive function used on D2xi()
double dfsD2xi(int i, std::shared_ptr<Node> node, int nodeTreeDepth,
               std::vector<double>* X) {
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

// Auxiliary DFS recursive function used on Da
void dfsDa(std::vector<double>* X, double parentsProduct,
           std::shared_ptr<Node> thisNode, int thisNodeTreeDepth, int* aIndex,
           std::vector<double>* target) {
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
    if (int(target->size()) != this->nTerms) {
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
              std::shared_ptr<Node> thisNode, int thisNodeTreeDepth,
              int* aIndex, std::vector<double>* target) {
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

// Auxiliary DFS recursive function used on DaDxi
void dfsDaD2xi(int i, std::vector<double>* X, double parentsProduct,
               std::shared_ptr<Node> thisNode, int thisNodeTreeDepth,
               int* aIndex, std::vector<double>* target) {
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

Maybe<Void> Poly::GetCoefficients(std::vector<double>* target) {
    Maybe<Void> r;
    if (int(target->size()) != this->nTerms) {
        r.isError = true;
        r.errMsg = "target must have same length as the number of terms";
        return r;
    }
    for (int i = 0; i < this->nTerms; i++) {
        (*target)[i] = this->leafNodes[i]->a;
    }
    return r;
}

Maybe<Void> Poly::SetCoefficients(std::vector<double>* target) {
    Maybe<Void> r;
    if (int(target->size()) != this->nTerms) {
        r.isError = true;
        r.errMsg = "target must have same length as the number of terms";
        return r;
    }
    for (int i = 0; i < this->nTerms; i++) {
        this->leafNodes[i]->a = target->at(i);
    }
    return r;
}

Poly::Poly() { this->isZero = true; }

Poly::Poly(int k) { this->isZero = true; }

// Copy constructor
Poly::Poly(const Poly& right) {
    if (right.isZero) {
        this->isZero = true;
        return;
    }
    this->isZero = false;
    assert(!this->Build(right.n, right.order).isError);
    for (int i = 0; i < int(right.leafNodes.size()); i++) {
        this->leafNodes[i]->a = right.leafNodes[i]->a;
    }
}
// Assignment
Poly& Poly::operator=(const Poly& right) {
    if (right.isZero) {
        this->isZero = true;
        return (*this);
    }
    this->isZero = false;
    assert(!this->Build(right.n, right.order).isError);
    for (int i = 0; i < int(right.leafNodes.size()); i++) {
        this->leafNodes[i]->a = right.leafNodes[i]->a;
    }
    return (*this);
}

Poly operator+(Poly const& left, Poly const& right) {
    if (left.isZero) {
        return right;
    }
    if (right.isZero) {
        return left;
    }
    assert(left.n == right.n);
    assert(left.order == right.order);
    Poly newP;
    assert(!newP.Build(left.n, left.order).isError);
    for (int i = 0; i < int(newP.leafNodes.size()); i++) {
        newP.leafNodes[i]->a = right.leafNodes[i]->a + left.leafNodes[i]->a;
    }
    return newP;
};
Poly operator-(Poly const& left, Poly const& right) {
    if (left.isZero) {
        return right;
    }
    if (right.isZero) {
        return left;
    }
    assert(left.n == right.n);
    assert(left.order == right.order);
    Poly newP;
    assert(!newP.Build(left.n, left.order).isError);
    for (int i = 0; i < int(newP.leafNodes.size()); i++) {
        newP.leafNodes[i]->a = left.leafNodes[i]->a - right.leafNodes[i]->a;
    }
    return newP;
};
Poly& Poly::operator+=(const Poly& right) {
    if (right.isZero) {
        return *this;
    }
    if (this->isZero) {
        this->Build(right.n, right.order);
    }
    assert(this->order == right.order);
    for (int i = 0; i < int(this->leafNodes.size()); i++) {
        this->leafNodes[i]->a += right.leafNodes[i]->a;
    }
    return *this;
}
Poly operator*(double x, const Poly& p) {
    for (int i = 0; i < int(p.leafNodes.size()); i++) {
        p.leafNodes[i]->a *= x;
    }
    return p;
}
Poly operator*(const Poly& p, double x) { return x * p; }

Poly operator+(double x, const Poly& p) {
    Poly newP = p;
    newP.leafNodes[newP.leafNodes.size() - 1]->a += x;
    return newP;
}
Poly operator+(const Poly& p, double x) {
    Poly newP = p;
    newP.leafNodes[newP.leafNodes.size() - 1]->a += x;
    return newP;
}
Poly operator-(const Poly& p, double x) {
    Poly newP = p;
    newP.leafNodes[newP.leafNodes.size() - 1]->a -= x;
    return newP;
}
Poly operator-(double x, const Poly& p) {
    Poly newP = p;
    newP = newP * (-1);
    newP.leafNodes[newP.leafNodes.size() - 1]->a += x;
    return newP;
}

bool Poly::operator==(Poly const& right) {
    if (right.isZero && this->isZero) {
        return true;
    }
    if (right.isZero && !this->isZero) {
        return false;
    }
    if (!right.isZero && this->isZero) {
        return false;
    }
    if (right.n != this->n || right.order != this->order) {
        return false;
    }
    for (int i = 0; i < int(this->leafNodes.size()); i++) {
        if (this->leafNodes[i]->a != right.leafNodes[i]->a) {
            return false;
        }
    }
    return true;
}

bool Poly::operator!=(Poly const& right) { return !((*this) == right); }
