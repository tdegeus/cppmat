
#include <catch/catch.hpp>

#define EQ(a,b) REQUIRE_THAT( (a), Catch::WithinAbs((b), 1.e-10) );

#define CPPMAT_NOCONVERT
// #include <cppmat/cppmat.h>
#include "../src/cppmat/cppmat.h"

#include <Eigen/Eigen>

#include "support.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatD;
typedef Eigen::Matrix<double, Eigen::Dynamic,              1, Eigen::ColMajor> ColD;

static const size_t ND = 3;

typedef cppmat::cartesian::tensor4 <double> T4;
typedef cppmat::cartesian::tensor2 <double> T2;
typedef cppmat::cartesian::tensor2s<double> T2s;
typedef cppmat::cartesian::tensor2d<double> T2d;
typedef cppmat::cartesian::vector  <double> V;

// =================================================================================================

TEST_CASE("cppmat::cartesian::vector", "var_cartesian_vector.h")
{

// =================================================================================================
// operations
// =================================================================================================

SECTION( "V.dot(V), length" )
{
  V A = V::Random(ND);

  double B = A.dot(A);

  double b = std::pow(A.length(), 2.);

  EQ( B, b );
}

// -------------------------------------------------------------------------------------------------

SECTION( "V.dot(T2)" )
{
  T2  I = T2 ::I(ND);
  V   A = V  ::Random(ND);

  V   B = A.dot(I);

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    EQ( B(i), A(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "V.dot(T2s)" )
{
  T2s I = T2s::I(ND);
  V   A = V  ::Random(ND);

  V   B = A.dot(I);

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    EQ( B(i), A(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "V.dot(T2d)" )
{
  T2d I = T2d::I(ND);
  V   A = V  ::Random(ND);

  V   B = A.dot(I);

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    EQ( B(i), A(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "V.dyadic(V)" )
{
  V  A = V::Random(ND);
  V  B = V::Random(ND);
  V  C = V::Random(ND);

  V  D = C.dot(A.dyadic(B));
  V  E = C.dot(A) * B;

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    EQ( D(i), E(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "V.cross(V)" )
{
  V  A = V::Random(3);
  V  B = V::Random(3);
  V  C = A.cross(B);

  double d = A.dot(C);
  double e = B.dot(C);

  EQ( d, 0.0 );
  EQ( e, 0.0 );
}

// =================================================================================================

}
