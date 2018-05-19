
#include <catch/catch.hpp>

#define EQ(a,b) REQUIRE_THAT( (a), Catch::WithinAbs((b), 1.e-10) );

#define CPPMAT_NOCONVERT
// #include <cppmat/cppmat.h>
#include "../src/cppmat/cppmat.h"

#include <Eigen/Eigen>

#include "support.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatD;
typedef Eigen::Matrix<double, Eigen::Dynamic,              1, Eigen::ColMajor> ColD;

static const size_t M = 11;
static const size_t N = 11;

typedef cppmat::diagonal ::matrix<double> dMat;
typedef cppmat::symmetric::matrix<double> sMat;
typedef cppmat           ::matrix<double>  Mat;

// =================================================================================================

TEST_CASE("cppmat::diagonal::matrix", "matrix.h")
{

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "matrix += matrix" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));
  MatD b = makeDiagonal(MatD::Random(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) += b(i,j);

  A += B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix -= matrix" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));
  MatD b = makeDiagonal(MatD::Random(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) -= b(i,j);

  A -= B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix *= matrix" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));
  MatD b = makeDiagonal(MatD::Random(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) *= b(i,j);

  A *= B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "matrix *= scalar" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));
  double b = a(0,0);

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) *= b;

  A *= b;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix /= scalar" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));
  double b = a(0,0) + 1.0;

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) /= b;

  A /= b;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "matrix + matrix" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));
  MatD b = makeDiagonal(MatD::Random(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) + b(i,j);

  dMat C = A + B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix - matrix" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));
  MatD b = makeDiagonal(MatD::Random(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) - b(i,j);

  dMat C = A - B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix * matrix" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));
  MatD b = makeDiagonal(MatD::Random(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) * b(i,j);

  dMat C = A * B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "matrix * scalar" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));
  double b = a(0,0);

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) * b;

  dMat C = A * b;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix / scalar" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));
  double b = a(0,0) + 1.0;

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) / b;

  dMat C = A / b;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "scalar * matrix" )
{
  MatD b = makeDiagonal(MatD::Random(M,N));
  double a = b(0,0);

  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a * b(i,j);

  dMat C = a * B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// =================================================================================================
// algebra
// =================================================================================================

SECTION( "min" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());

  double C = A.min();

  double c = a.minCoeff();

  EQ( c, C );
}

// -------------------------------------------------------------------------------------------------

SECTION( "max" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());

  double C = A.max();

  double c = a.maxCoeff();

  EQ( c, C );
}

// -------------------------------------------------------------------------------------------------

SECTION( "sum" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());

  double C = A.sum();

  double c = 0.0;

  for ( auto i = 0 ; i < a.rows() ; ++i )
    for ( auto j = 0 ; j < a.cols() ; ++j )
      c += a(i,j);

  EQ( c, C );
}

// -------------------------------------------------------------------------------------------------

SECTION( "mean" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());

  double C = A.mean();

  double c = 0.0;

  for ( auto i = 0 ; i < a.rows() ; ++i )
    for ( auto j = 0 ; j < a.cols() ; ++j )
      c += a(i,j);

  c /= static_cast<double>(M*N);

  EQ( c, C );
}

// -------------------------------------------------------------------------------------------------

SECTION( "average" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));
  MatD b = makeDiagonal(MatD::Random(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  double C = A.average(B);

  double c = 0.0;
  double d = 0.0;

  for ( auto i = 0 ; i < a.rows() ; ++i ) {
    for ( auto j = 0 ; j < a.cols() ; ++j ) {
      c += b(i,j) * a(i,j);
      d += b(i,j);
    }
  }

  c /= d;

  EQ( c, C );
}

// =================================================================================================
// absolute value
// =================================================================================================

SECTION( "abs" )
{
  MatD a = makeDiagonal(MatD::Random(M,N) - MatD::Constant(M,N,.5));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());

  dMat C = A.abs();

  MatD c = a.cwiseAbs();

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// =================================================================================================
// index operators
// =================================================================================================

SECTION( "decompress" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());

  std::vector<size_t> idx = A.decompress(A.compress(1,1));
  std::vector<size_t> jdx = {1,1};

  REQUIRE( idx == jdx );
}

// =================================================================================================

}
