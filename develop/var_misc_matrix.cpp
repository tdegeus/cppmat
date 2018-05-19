
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

typedef cppmat::diagonal::matrix<double> dMat;
typedef cppmat::symmetric::matrix<double> sMat;
typedef cppmat::matrix<double> Mat;

// =================================================================================================

TEST_CASE("cppmat::misc::matrix", "matrix.h")
{

// =================================================================================================
// extra arithmetic operators : cppmat::matrix
// =================================================================================================

SECTION( "matrix *= symmetric::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeSymmetric(MatD::Random(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) *=  b(i,j);

  A *= B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix /= symmetric::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeSymmetric(MatD::Random(M,N) + MatD::Ones(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) /=  b(i,j);

  A /= B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix += symmetric::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeSymmetric(MatD::Random(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) +=  b(i,j);

  A += B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix -= symmetric::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeSymmetric(MatD::Random(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) -=  b(i,j);

  A -= B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix *= diagonal::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeDiagonal(MatD::Random(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) *=  b(i,j);

  A *= B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix += diagonal::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeDiagonal(MatD::Random(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) +=  b(i,j);

  A += B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix -= diagonal::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeDiagonal(MatD::Random(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) -=  b(i,j);

  A -= B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// =================================================================================================
// extra arithmetic operators : cppmat::symmetric::matrix
// =================================================================================================

SECTION( "symmetric::matrix *= diagonal::matrix" )
{
  MatD a = makeSymmetric(MatD::Random(M,N));
  MatD b = makeDiagonal (MatD::Random(M,N));

  sMat A = sMat::CopyDense(M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) *=  b(i,j);

  A *= B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "symmetric::matrix += diagonal::matrix" )
{
  MatD a = makeSymmetric(MatD::Random(M,N));
  MatD b = makeDiagonal (MatD::Random(M,N));

  sMat A = sMat::CopyDense(M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) +=  b(i,j);

  A += B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "symmetric::matrix -= diagonal::matrix" )
{
  MatD a = makeSymmetric(MatD::Random(M,N));
  MatD b = makeDiagonal (MatD::Random(M,N));

  sMat A = sMat::CopyDense(M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) -=  b(i,j);

  A -= B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// =================================================================================================
// extra arithmetic operators : cppmat::diagonal::matrix
// =================================================================================================

SECTION( "diagonal::matrix *= matrix" )
{
  MatD a = makeDiagonal (MatD::Random(M,N));
  MatD b = MatD::Random(M,N);

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
   Mat B =  Mat::Copy     (M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) *=  b(i,j);

  A *= B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "diagonal::matrix /= matrix" )
{
  MatD a = makeDiagonal (MatD::Random(M,N));
  MatD b = MatD::Random(M,N) + MatD::Ones(M,N);

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
   Mat B =  Mat::Copy     (M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) /=  b(i,j);

  A /= B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "diagonal::matrix *= symmetric::matrix" )
{
  MatD a = makeDiagonal (MatD::Random(M,N));
  MatD b = makeSymmetric(MatD::Random(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) *=  b(i,j);

  A *= B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "diagonal::matrix /= symmetric::matrix" )
{
  MatD a = makeDiagonal (MatD::Random(M,N));
  MatD b = makeSymmetric(MatD::Random(M,N) + MatD::Ones(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) /=  b(i,j);

  A /= B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// =================================================================================================
// extra external arithmetic operators -> cppmat::matrix
// =================================================================================================

SECTION( "matrix * symmetric::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeSymmetric(MatD::Random(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) * b(i,j);

  Mat C = A * B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix / symmetric::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeSymmetric(MatD::Random(M,N) + MatD::Ones(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) / b(i,j);

  Mat C = A / B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix + symmetric::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeSymmetric(MatD::Random(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) + b(i,j);

  Mat C = A + B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix - symmetric::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeSymmetric(MatD::Random(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) - b(i,j);

  Mat C = A - B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix + diagonal::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeDiagonal(MatD::Random(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) + b(i,j);

  Mat C = A + B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix - diagonal::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeDiagonal(MatD::Random(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) - b(i,j);

  Mat C = A - B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "symmetric::matrix * matrix" )
{
  MatD a = makeSymmetric(MatD::Random(M,N));
  MatD b = MatD::Random(M,N);

  sMat A = sMat::CopyDense(M, N, a.data(), a.data()+a.size());
   Mat B =  Mat::Copy     (M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) * b(i,j);

  Mat C = A * B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "symmetric::matrix / matrix" )
{
  MatD a = makeSymmetric(MatD::Random(M,N) + MatD::Ones(M,N));
  MatD b = MatD::Random(M,N);

  sMat A = sMat::CopyDense(M, N, a.data(), a.data()+a.size());
   Mat B =  Mat::Copy     (M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) / b(i,j);

  Mat C = A / B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "symmetric::matrix + matrix" )
{
  MatD a = makeSymmetric(MatD::Random(M,N));
  MatD b = MatD::Random(M,N);

  sMat A = sMat::CopyDense(M, N, a.data(), a.data()+a.size());
   Mat B =  Mat::Copy     (M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) + b(i,j);

  Mat C = A + B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "symmetric::matrix - matrix" )
{
  MatD a = makeSymmetric(MatD::Random(M,N));
  MatD b = MatD::Random(M,N);

  sMat A = sMat::CopyDense(M, N, a.data(), a.data()+a.size());
   Mat B =  Mat::Copy     (M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) - b(i,j);

  Mat C = A - B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "diagonal::matrix + matrix" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));
  MatD b = MatD::Random(M,N);

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
   Mat B =  Mat::Copy     (M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) + b(i,j);

  Mat C = A + B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "diagonal::matrix - matrix" )
{
  MatD a = makeDiagonal(MatD::Random(M,N));
  MatD b = MatD::Random(M,N);

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
   Mat B =  Mat::Copy     (M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) - b(i,j);

  Mat C = A - B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// =================================================================================================
// extra external arithmetic operators -> cppmat::symmetric::matrix
// =================================================================================================

SECTION( "symmetric::matrix + diagonal::matrix" )
{
  MatD a = makeSymmetric(MatD::Random(M,N));
  MatD b = makeDiagonal (MatD::Random(M,N));

  sMat A = sMat::CopyDense(M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) + b(i,j);

  sMat C = A + B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "symmetric::matrix - diagonal::matrix" )
{
  MatD a = makeSymmetric(MatD::Random(M,N));
  MatD b = makeDiagonal (MatD::Random(M,N));

  sMat A = sMat::CopyDense(M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) - b(i,j);

  sMat C = A - B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "diagonal::matrix + scalar" )
{
  MatD   a = makeDiagonal(MatD::Random(M,N));
  double b = a(0,0);

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) + b;

  sMat C = A + b;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "diagonal::matrix - scalar" )
{
  MatD   a = makeDiagonal(MatD::Random(M,N));
  double b = a(0,0);

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) - b;

  sMat C = A - b;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "diagonal::matrix + symmetric::matrix" )
{
  MatD a = makeSymmetric(MatD::Random(M,N));
  MatD b = makeDiagonal (MatD::Random(M,N));

  sMat A = sMat::CopyDense(M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = b(i,j) + a(i,j);

  sMat C = B + A;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "diagonal::matrix - symmetric::matrix" )
{
  MatD a = makeSymmetric(MatD::Random(M,N));
  MatD b = makeDiagonal (MatD::Random(M,N));

  sMat A = sMat::CopyDense(M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = b(i,j) - a(i,j);

  sMat C = B - A;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "scalar + diagonal::matrix" )
{
  MatD   a = makeDiagonal(MatD::Random(M,N));
  double b = a(0,0);

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = b + a(i,j);

  sMat C = b + A;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "scalar - diagonal::matrix" )
{
  MatD   a = makeDiagonal(MatD::Random(M,N));
  double b = a(0,0);

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = b - a(i,j);

  sMat C = b - A;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// =================================================================================================
// extra external arithmetic operators -> cppmat::diagonal::matrix
// =================================================================================================

SECTION( "diagonal::matrix * matrix" )
{
  MatD a = makeDiagonal(MatD::Random(M,N) + MatD::Ones(M,N));
  MatD b = MatD::Random(M,N);

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
   Mat B =  Mat::Copy     (M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) * b(i,j);

  dMat C = A * B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "diagonal::matrix / matrix" )
{
  MatD a = makeDiagonal(MatD::Random(M,N) + MatD::Ones(M,N));
  MatD b = MatD::Random(M,N);

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
   Mat B =  Mat::Copy     (M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) / b(i,j);

  dMat C = A / B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "diagonal::matrix * symmetric::matrix" )
{
  MatD a = makeDiagonal (MatD::Random(M,N) + MatD::Ones(M,N));
  MatD b = makeSymmetric(MatD::Random(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) * b(i,j);

  dMat C = A * B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "diagonal::matrix / symmetric::matrix" )
{
  MatD a = makeDiagonal (MatD::Random(M,N) + MatD::Ones(M,N));
  MatD b = makeSymmetric(MatD::Random(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) / b(i,j);

  dMat C = A / B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "diagonal::matrix * matrix" )
{
  MatD a = makeDiagonal(MatD::Random(M,N) + MatD::Ones(M,N));
  MatD b = MatD::Random(M,N);

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
   Mat B =  Mat::Copy     (M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = b(i,j) * a(i,j);

  dMat C = B * A;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "diagonal::matrix * symmetric::matrix" )
{
  MatD a = makeDiagonal (MatD::Random(M,N) + MatD::Ones(M,N));
  MatD b = makeSymmetric(MatD::Random(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = b(i,j) * a(i,j);

  dMat C = B * A;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// =================================================================================================

}
