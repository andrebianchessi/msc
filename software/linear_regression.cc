#include "linear_regression.h"

#include <boost/math/special_functions/binomial.hpp>
#include <queue>

#include "maybe.h"
#include "utils.h"

// Returns number of polynomial terms given the number of variables and the
// max degree
int polyTerms(int nVars, int maxOrder) {
    // Source:
    // https://mathoverflow.net/questions/225953/number-of-polynomial-terms-for-certain-degree-and-certain-number-of-variables/225963#225963?newreg=2a0208ceb740461d8eaa21e304b0e341
    return int(
        boost::math::binomial_coefficient<double>(nVars + maxOrder, maxOrder));
}

int LinReg::nTerms() { return polyTerms(this->XSize, this->order); }

Maybe<LinReg> LinReg::NewLinReg(int XSize, int order) {
    Maybe<LinReg> r;
    if (XSize < 0 || order < 0) {
        r.isError = true;
        r.errMsg = "XSize and order must be >=0";
        return r;
    }
    LinReg lr;
    lr.XSize = XSize;
    lr.order = order;

    lr.coefficients = new Node();
    lr.coefficients->l = order;
    std::queue<Node*> q;
    q.push(lr.coefficients);

    for (int i = 0; i < XSize; i++) {
        for (int qI = q.size(); qI > 0; qI--) {
            Node* n = q.front();
            q.pop();

            for (int e = n->l; e >= 0; e--) {
                Node* newNode = new Node();
                newNode->e = e;
                newNode->l = n->l - e;
                if (i == XSize - 1) {
                    newNode->a = 0.0;
                }
                n->children.push_back(newNode);
                q.push(newNode);
            }
        }
    }

    r.val = lr;
    return r;
};

double dfs(Node* n, int index, std::vector<double>* X) {
    if (n->IsLeaf()) {
        return n->a * std::pow(X->at(index), n->e);
    }
    double sumOfChildren = 0;
    for (int c = 0; c < int(n->children.size()); c++) {
        if (index == -1) {
            sumOfChildren += dfs(n->children[c], index + 1, X);
        } else {
            sumOfChildren += std::pow(X->at(index), n->e) *
                             dfs(n->children[c], index + 1, X);
        }
    }
    return sumOfChildren;
}
Maybe<double> LinReg::operator()(std::vector<double>* X) {
    Maybe<double> r;
    if (int(X->size()) != this->XSize) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }
    r.val = dfs(this->coefficients, -1, X);
    return r;
};

double dfsD(Node* n, int di, int index, std::vector<double>* X) {
    if (n->IsLeaf()) {
        if (index != di) {
            return n->a * std::pow(X->at(index), n->e);
        }
        if (n->e == 0) {
            return 0;
        }
        return n->a * n->e * std::pow(X->at(index), n->e - 1);
    }
    double sumOfChildren = 0;
    for (int c = 0; c < int(n->children.size()); c++) {
        if (index == -1) {
            sumOfChildren += dfsD(n->children[c], di, index + 1, X);
        } else {
            if (index != di) {
                sumOfChildren += std::pow(X->at(index), n->e) *
                                 dfsD(n->children[c], di, index + 1, X);
            } else {
                if (n->e != 0) {
                    sumOfChildren += n->e * std::pow(X->at(index), n->e - 1) *
                                     dfsD(n->children[c], di, index + 1, X);
                }
            }
        }
    }
    return sumOfChildren;
}
Maybe<double> LinReg::D(int di, std::vector<double>* X) {
    Maybe<double> r;
    if (int(X->size()) != this->XSize) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }
    if (di < 0 || di >= this->XSize) {
        r.isError = true;
        r.errMsg = "invalid i";
        return r;
    }
    r.val = dfsD(this->coefficients, di, -1, X);
    return r;
}

double dfsD2(Node* n, int di, int index, std::vector<double>* X) {
    if (n->IsLeaf()) {
        if (index != di) {
            return n->a * std::pow(X->at(index), n->e);
        }
        if (n->e <= 1) {
            return 0;
        }
        return n->a * n->e * (n->e - 1) * std::pow(X->at(index), n->e - 2);
    }
    double sumOfChildren = 0;
    for (int c = 0; c < int(n->children.size()); c++) {
        if (index == -1) {
            sumOfChildren += dfsD2(n->children[c], di, index + 1, X);
        } else {
            if (index != di) {
                sumOfChildren += std::pow(X->at(index), n->e) *
                                 dfsD2(n->children[c], di, index + 1, X);
            } else {
                if (n->e >= 2) {
                    sumOfChildren += n->e * (n->e - 1) *
                                     std::pow(X->at(index), n->e - 2) *
                                     dfsD2(n->children[c], di, index + 1, X);
                }
            }
        }
    }
    return sumOfChildren;
}
Maybe<double> LinReg::D2(int di, std::vector<double>* X) {
    Maybe<double> r;
    if (int(X->size()) != this->XSize) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }
    if (di < 0 || di >= this->XSize) {
        r.isError = true;
        r.errMsg = "invalid i";
        return r;
    }
    r.val = dfsD2(this->coefficients, di, -1, X);
    return r;
}

Maybe<int> LinReg::CoefficientAt(std::vector<int> powers) {
    Maybe<int> r;

    if (int(powers.size()) > this->XSize) {
        r.isError = true;
        r.errMsg = "powers must be same length of X";
        return r;
    }

    int pSum = 0;
    for (int i = 0; i < int(powers.size()); i++) {
        if (powers[i] < 0) {
            r.isError = true;
            r.errMsg = "powers must be >=0";
            return r;
        }
        pSum += powers[i];
    }
    if (pSum > this->order) {
        r.isError = true;
        r.errMsg = "sum of powers larger than order of regression";
        return r;
    }

    Node* n = this->coefficients;
    int pIndex = 0;
    int childrenIndex;
    while (pIndex < int(powers.size())) {
        childrenIndex = n->l - powers[pIndex];
        n = n->children[childrenIndex];
        pIndex += 1;
    }
    r.val = n->a;
    return r;
}