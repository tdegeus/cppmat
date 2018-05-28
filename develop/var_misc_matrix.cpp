
#include "support.h"

static const size_t M = 11;
static const size_t N = 11;

typedef cppmat::diagonal ::matrix<double> dMat;
typedef cppmat::symmetric::matrix<double> sMat;
typedef cppmat::           matrix<double> Mat;

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

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) *= b(i,j);

  A *= B;

  Equal(A, a);
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix /= symmetric::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeSymmetric(MatD::Random(M,N) + MatD::Ones(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) /= b(i,j);

  A /= B;

  Equal(A, a);
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix += symmetric::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeSymmetric(MatD::Random(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) += b(i,j);

  A += B;

  Equal(A, a);
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix -= symmetric::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeSymmetric(MatD::Random(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) -= b(i,j);

  A -= B;

  Equal(A, a);
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix *= diagonal::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeDiagonal(MatD::Random(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) *= b(i,j);

  A *= B;

  Equal(A, a);
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix += diagonal::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeDiagonal(MatD::Random(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) += b(i,j);

  A += B;

  Equal(A, a);
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix -= diagonal::matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = makeDiagonal(MatD::Random(M,N));

   Mat A =  Mat::Copy     (M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) -= b(i,j);

  A -= B;

  Equal(A, a);
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

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) *= b(i,j);

  A *= B;

  Equal(A, a);
}

// -------------------------------------------------------------------------------------------------

SECTION( "symmetric::matrix += diagonal::matrix" )
{
  MatD a = makeSymmetric(MatD::Random(M,N));
  MatD b = makeDiagonal (MatD::Random(M,N));

  sMat A = sMat::CopyDense(M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) += b(i,j);

  A += B;

  Equal(A, a);
}

// -------------------------------------------------------------------------------------------------

SECTION( "symmetric::matrix -= diagonal::matrix" )
{
  MatD a = makeSymmetric(MatD::Random(M,N));
  MatD b = makeDiagonal (MatD::Random(M,N));

  sMat A = sMat::CopyDense(M, N, a.data(), a.data()+a.size());
  dMat B = dMat::CopyDense(M, N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) -= b(i,j);

  A -= B;

  Equal(A, a);
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

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) *= b(i,j);

  A *= B;

  Equal(A, a);
}

// -------------------------------------------------------------------------------------------------

SECTION( "diagonal::matrix /= matrix" )
{
  MatD a = makeDiagonal (MatD::Random(M,N));
  MatD b = MatD::Random(M,N) + MatD::Ones(M,N);

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
   Mat B =  Mat::Copy     (M, N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) /= b(i,j);

  A /= B;

  Equal(A, a);
}

// -------------------------------------------------------------------------------------------------

SECTION( "diagonal::matrix *= symmetric::matrix" )
{
  MatD a = makeDiagonal (MatD::Random(M,N));
  MatD b = makeSymmetric(MatD::Random(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) *= b(i,j);

  A *= B;

  Equal(A, a);
}

// -------------------------------------------------------------------------------------------------

SECTION( "diagonal::matrix /= symmetric::matrix" )
{
  MatD a = makeDiagonal (MatD::Random(M,N));
  MatD b = makeSymmetric(MatD::Random(M,N) + MatD::Ones(M,N));

  dMat A = dMat::CopyDense(M, N, a.data(), a.data()+a.size());
  sMat B = sMat::CopyDense(M, N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) /= b(i,j);

  A /= B;

  Equal(A, a);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
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

  Equal(C, c);
}

// =================================================================================================

}
