#include <gtest/gtest.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// Check that boost vectors and matrices are correctly imported

using namespace boost::numeric::ublas;

TEST(BoostTest, VectorTest){
  vector<int> v(3);
  v[0] = 1;
  ASSERT_EQ(v[0],1);
  ASSERT_EQ(v[1],0);
  ASSERT_EQ(v[2],0);
  ASSERT_EQ(v.size(),3);
}

TEST(BoostTest, MatrixTest){
  matrix<int> m(2,2);
  m(0,0) = 1;
  m(0,1) = 2;
  m(1,0) = 3;
  m(1,1) = 4;
  ASSERT_EQ(m(0,0),1);
  ASSERT_EQ(m(1,1),4);
  ASSERT_EQ(m.size1(),2);
  ASSERT_EQ(m.size2(),2);
}

TEST(BoostTest, BoundedMatrixTest){
  bounded_matrix<int,2,2> m(2,2);
  m(0,0) = 1;
  m(0,1) = 2;
  m(1,0) = 3;
  m(1,1) = 4;
  ASSERT_EQ(m(0,0),1);
  ASSERT_EQ(m(1,1),4);
  ASSERT_EQ(m.size1(),2);
  ASSERT_EQ(m.size2(),2);
}