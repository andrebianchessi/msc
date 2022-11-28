#include "model.h"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "maybe.h"
#include "utils.h"

class TestModel : public Model {
    // Linear regression model of second order
    // y(x) = ax^2 + bx + c
    // Training data:
    // x | y_real
    // 0 | 5
    // 1 | 25
    // 2 | 71
    // Exact solution: a = 13, b = 7, c = 5
   public:
    double a, b, c;
    TestModel() {
        this->a = 0;
        this->b = 0;
        this->c = 0;
    }

    int nParameters() override { return 3; }

    Maybe<double> operator()(std::vector<double>* X) {
        Maybe<double> r;
        if (X->size() != 1) {
            r.errMsg = "Invalid X";
            r.isError = true;
            return r;
        }
        r.val = this->a * X->at(0) * X->at(0) + this->b * X->at(0) + this->c;
        return r;
    };

   private:
    Maybe<Void> GetParameters(std::vector<double>* target) override {
        Maybe<Void> r;
        if (target->size() != 3) {
            r.errMsg = "Invalid target length";
            r.isError = true;
            return r;
        }
        target->at(0) = this->a;
        target->at(1) = this->b;
        target->at(2) = this->c;
        return r;
    };

    Maybe<Void> SetParameters(std::vector<double>* parameters) override {
        Maybe<Void> r;
        if (parameters->size() != 3) {
            r.errMsg = "Invalid parameters length";
            r.isError = true;
            return r;
        }
        this->a = parameters->at(0);
        this->b = parameters->at(1);
        this->c = parameters->at(2);
        return r;
    };

    FRIEND_TEST(ModelTest, LossTest);
    double Loss() override {
        // Standard quadratic loss
        // L = Sum ( (y(x)-y_real)^2 )
        double l = 0;
        std::vector<double> X;

        // x | y
        // 0 | 5
        // 1 | 25
        // 2 | 71
        X = {0};
        l += pow((*this)(&X).val - (5), 2);
        X = {1};
        l += pow((*this)(&X).val - (25), 2);
        X = {2};
        l += pow((*this)(&X).val - (71), 2);

        return l;
    };

    FRIEND_TEST(ModelTest, LossGradientTest);
    std::vector<double> LossGradient() override {
        // L = Sum ( (y(x)-y_real)^2 )
        // dL/di = Sum ( 2*(y(x)-y_real)* dy(x)/di  )
        std::vector<double> grad = std::vector<double>(3);

        std::vector<double> X = std::vector<double>(1);

        // dL/da = x^2
        X = {0};
        grad[0] += 2 * ((*this)(&X).val - (5)) * X[0] * X[0];
        X = {1};
        grad[0] += 2 * ((*this)(&X).val - (25)) * X[0] * X[0];
        X = {2};
        grad[0] += 2 * ((*this)(&X).val - (71)) * X[0] * X[0];

        // dL/db = x
        X = {0};
        grad[1] += 2 * ((*this)(&X).val - (5)) * X[0];
        X = {1};
        grad[1] += 2 * ((*this)(&X).val - (25)) * X[0];
        X = {2};
        grad[1] += 2 * ((*this)(&X).val - (71)) * X[0];

        // dL/dc = 1
        X = {0};
        grad[2] += 2 * ((*this)(&X).val - (5)) * 1.0;
        X = {1};
        grad[2] += 2 * ((*this)(&X).val - (25)) * 1.0;
        X = {2};
        grad[2] += 2 * ((*this)(&X).val - (71)) * 1.0;

        return grad;
    };
};

TEST(ModelTest, ConstructorTest) {
    TestModel m = TestModel();
    ASSERT_EQ(m.a, 0);
    ASSERT_EQ(m.b, 0);
    ASSERT_EQ(m.c, 0);
}

TEST(ModelTest, OperatorTest) {
    TestModel m = TestModel();
    auto x = std::vector<double>{1.0};
    auto v = m(&x);
    ASSERT_FALSE(v.isError);
    ASSERT_DOUBLE_EQ(v.val, 0.0);

    double a = 99.0;
    double b = 98.0;
    double c = 97.0;
    m.a = a;
    m.b = b;
    m.c = c;
    x = std::vector<double>{7.0};
    v = m(&x);
    ASSERT_FALSE(v.isError);
    ASSERT_DOUBLE_EQ(v.val, a * 7 * 7 + b * 7 + c);

    x = {1.0, 2.0};
    v = m(&x);
    ASSERT_TRUE(v.isError);
}

TEST(ModelTest, LossTest) {
    TestModel m = TestModel();
    double a = 99.0;
    double b = 98.0;
    double c = 97.0;
    m.a = a;
    m.b = b;
    m.c = c;
    auto loss = m.Loss();

    // Training data:
    // x | y
    // 0 | 5
    // 1 | 25
    // 2 | 71
    double expected = 0;

    expected += pow(a * 0 * 0 + b * 0 + c - 5, 2);
    expected += pow(a * 1 * 1 + b * 1 + c - 25, 2);
    expected += pow(a * 2 * 2 + b * 2 + c - 71, 2);
    ASSERT_DOUBLE_EQ(loss, expected);
}

TEST(ModelTest, LossGradientTest) {
    TestModel m = TestModel();
    double a = 99.0;
    double b = 98.0;
    double c = 97.0;
    m.a = a;
    m.b = b;
    m.c = c;

    auto lossGrad = m.LossGradient();
    // Training data:
    // x | y
    // 0 | 5
    // 1 | 25
    // 2 | 71
    double da = 0;
    da += 2 * (a * 0 * 0 + b * 0 + c - 5) * 0 * 0;
    da += 2 * (a * 1 * 1 + b * 1 + c - 25) * 1 * 1;
    da += 2 * (a * 2 * 2 + b * 2 + c - 71) * 2 * 2;
    ASSERT_DOUBLE_EQ(da, lossGrad[0]);

    double db = 0;
    db += 2 * (a * 0 * 0 + b * 0 + c - 5) * 0;
    db += 2 * (a * 1 * 1 + b * 1 + c - 25) * 1;
    db += 2 * (a * 2 * 2 + b * 2 + c - 71) * 2;
    ASSERT_DOUBLE_EQ(db, lossGrad[1]);

    double dc = 0;
    dc += 2 * (a * 0 * 0 + b * 0 + c - 5) * 1;
    dc += 2 * (a * 1 * 1 + b * 1 + c - 25) * 1;
    dc += 2 * (a * 2 * 2 + b * 2 + c - 71) * 1;
    ASSERT_DOUBLE_EQ(dc, lossGrad[2]);
}

TEST(ModelTest, TrainTest) {
    TestModel m = TestModel();
    auto status = m.Train(0.01, 0.01, true);
    ASSERT_TRUE(RelativeAbsError(m.a, 13) < 0.001);
    ASSERT_TRUE(RelativeAbsError(m.b, 7) < 0.001);
    ASSERT_TRUE(RelativeAbsError(m.c, 5) < 0.001);
}