
#include <catch/catch.hpp>

#define EQ(a,b) REQUIRE_THAT( (a), Catch::WithinAbs((b), 1.e-10) );

#define CPPMAT_NOCONVERT
// #include <cppmat/cppmat.h>
#include "../src/cppmat/cppmat.h"

#include <Eigen/Eigen>

#include "support.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatD;
typedef Eigen::Matrix<double, Eigen::Dynamic,              1, Eigen::ColMajor> ColD;

static const size_t ND = 11;

typedef cppmat::cartesian::tensor4 <double> T4;
typedef cppmat::cartesian::tensor2 <double> T2;
typedef cppmat::cartesian::tensor2s<double> T2s;
typedef cppmat::cartesian::tensor2d<double> T2d;
typedef cppmat::cartesian::vector  <double> V;

// =================================================================================================

TEST_CASE("cppmat::cartesian::tensor2d", "var_cartesian_tensor2s.h")
{

// =================================================================================================
// unit tensors
// =================================================================================================

SECTION( "I, T2d.dot(T2)" )
{
  T2d I = T2d::I(ND);
  T2  A = T2 ::Random(ND);

  T2  B = I.dot(A);

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    for ( size_t j = 0 ; j < B.ndim() ; ++j )
      EQ( B(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "I, T2d.dot(T2s)" )
{
  T2d I = T2d::I(ND);
  T2s A = T2s::Random(ND);

  T2  B = I.dot(A);

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    for ( size_t j = 0 ; j < B.ndim() ; ++j )
      EQ( B(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "I, T2d.dot(T2d)" )
{
  T2d I = T2d::I(ND);
  T2d A = T2d::Random(ND);

  T2d B = I.dot(A);

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    for ( size_t j = 0 ; j < B.ndim() ; ++j )
      EQ( B(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "I, T2d.dot(V)" )
{
  T2d I = T2d::I(ND);
  V   A = V  ::Random(ND);

  V   B = I.dot(A);

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    EQ( B(i), A(i) );
}

// =================================================================================================
// tensor products
// =================================================================================================

SECTION( "I, T2d.dot(T4)" )
{
  T2d A = T2d::Random(ND);
  T4  I = T4 ::I(ND);

  T2  B = A.ddot(I);

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    for ( size_t j = 0 ; j < B.ndim() ; ++j )
      EQ( B(i,j), A(i,j) );
}

// =================================================================================================
// tensor products
// =================================================================================================

SECTION( "T2d.ddot(T2), T2d.dot(T2), T2d.trace()" )
{
  T2d A = T2d::Random(ND);
  T2  B = T2 ::Random(ND);

  double C = A.ddot(B);

  double c = A.dot(B).trace();

  EQ( C, c );
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2d.ddot(T2s), T2d.dot(T2s), T2d.trace()" )
{
  T2d A = T2d::Random(ND);
  T2s B = T2s::Random(ND);

  double C = A.ddot(B);

  double c = A.dot(B).trace();

  EQ( C, c );
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2d.ddot(T2d), T2d.dot(T2d), T2d.trace()" )
{
  T2d A = T2d::Random(ND);
  T2d B = T2d::Random(ND);

  double C = A.ddot(B);

  double c = A.dot(B).trace();

  EQ( C, c );
}

// =================================================================================================
// tensor products
// =================================================================================================

SECTION( "T2d.dyadic(T2)" )
{
  T2d A = T2d::I(ND);
  T2  B = T2 ::I(ND);

  T4  C = A.dyadic(B);

  T4  D = T4::II(ND);

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    for ( size_t j = 0 ; j < B.ndim() ; ++j )
      for ( size_t k = 0 ; k < B.ndim() ; ++k )
        for ( size_t l = 0 ; l < B.ndim() ; ++l )
          EQ( C(i,j,k,l), D(i,j,k,l) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2d.dyadic(T2s)" )
{
  T2d A = T2d::I(ND);
  T2s B = T2s::I(ND);

  T4  C = A.dyadic(B);

  T4  D = T4::II(ND);

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    for ( size_t j = 0 ; j < B.ndim() ; ++j )
      for ( size_t k = 0 ; k < B.ndim() ; ++k )
        for ( size_t l = 0 ; l < B.ndim() ; ++l )
          EQ( C(i,j,k,l), D(i,j,k,l) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2d.dyadic(T2d)" )
{
  T2d A = T2d::I(ND);
  T2d B = T2d::I(ND);

  T4  C = A.dyadic(B);

  T4  D = T4::II(ND);

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    for ( size_t j = 0 ; j < B.ndim() ; ++j )
      for ( size_t k = 0 ; k < B.ndim() ; ++k )
        for ( size_t l = 0 ; l < B.ndim() ; ++l )
          EQ( C(i,j,k,l), D(i,j,k,l) );
}

// =================================================================================================
// transpositions
// =================================================================================================

SECTION( "T" )
{
  T2d A = T2d::Random(ND);

  T2d B = A.T();

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    for ( size_t j = 0 ; j < B.ndim() ; ++j )
      EQ( B(i,j), A(j,i) );
}

// =================================================================================================
// determinant
// =================================================================================================

SECTION("T2d.det() -- 3D")
{
  MatD   a = makeDiagonal(MatD::Random(3,3));
  double c = a.determinant();

  T2s    A = T2s::CopyDense(3, a.data(), a.data()+a.size());
  double C = A.det();

  EQ( C, c );
}

// -------------------------------------------------------------------------------------------------

SECTION("T2d.det() -- 2D")
{
  MatD   a = makeDiagonal(MatD::Random(2,2));
  double c = a.determinant();

  T2s    A = T2s::CopyDense(2, a.data(), a.data()+a.size());
  double C = A.det();

  EQ( C, c );
}

// =================================================================================================
// inverse
// =================================================================================================

SECTION("T2d.inv() -- 3D")
{
  MatD a = makeDiagonal(MatD::Random(3,3));
  MatD c = a.inverse();

  T2s A = T2s::CopyDense(3, a.data(), a.data()+a.size());
  T2s C = A.inv();

  for ( size_t i = 0 ; i < A.ndim() ; ++i )
    for ( size_t j = 0 ; j < A.ndim() ; ++j )
      EQ( C(i,j), c(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION("T2d.inv() -- 2D")
{
  MatD a = makeDiagonal(MatD::Random(2,2));
  MatD c = a.inverse();

  T2s A = T2s::CopyDense(2, a.data(), a.data()+a.size());
  T2s C = A.inv();

  for ( size_t i = 0 ; i < A.ndim() ; ++i )
    for ( size_t j = 0 ; j < A.ndim() ; ++j )
      EQ( C(i,j), c(i,j) );
}

// =================================================================================================

}
