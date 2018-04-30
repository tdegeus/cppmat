
#include <catch/catch.hpp>

#define EQ(a,b) REQUIRE_THAT( (a), Catch::WithinAbs((b), 1.e-10) );

#define CPPMAT_NOCONVERT
#include <cppmat.h>

#include <Eigen/Eigen>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatD;
typedef Eigen::Matrix<double, Eigen::Dynamic,              1, Eigen::ColMajor> ColD;

static const size_t M = 11;
static const size_t N = 9;

typedef cppmat::matrix<double> Mat;
typedef cppmat::vector<double> Vec;

// =================================================================================================

TEST_CASE("cppmat::matrix", "matrix.h")
{

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "matrix += matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());
  Mat B = Mat::Copy(M,N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) += b(i,j);

  A += B;

  REQUIRE( static_cast<size_t>(a.size()) == A.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix -= matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());
  Mat B = Mat::Copy(M,N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) -= b(i,j);

  A -= B;

  REQUIRE( static_cast<size_t>(a.size()) == A.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix *= matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());
  Mat B = Mat::Copy(M,N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) *= b(i,j);

  A *= B;

  REQUIRE( static_cast<size_t>(a.size()) == A.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix /= matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N) + MatD::Ones(M,N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());
  Mat B = Mat::Copy(M,N, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) /= b(i,j);

  A /= B;

  REQUIRE( static_cast<size_t>(a.size()) == A.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "matrix += scalar" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) += b;

  A += b;

  REQUIRE( static_cast<size_t>(a.size()) == A.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix -= scalar" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) -= b;

  A -= b;

  REQUIRE( static_cast<size_t>(a.size()) == A.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix *= scalar" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) *= b;

  A *= b;

  REQUIRE( static_cast<size_t>(a.size()) == A.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix /= scalar" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0) + 1.0;

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) /= b;

  A /= b;

  REQUIRE( static_cast<size_t>(a.size()) == A.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "matrix + matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());
  Mat B = Mat::Copy(M,N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) + b(i,j);

  Mat C = A + B;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix - matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());
  Mat B = Mat::Copy(M,N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) - b(i,j);

  Mat C = A - B;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix * matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());
  Mat B = Mat::Copy(M,N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) * b(i,j);

  Mat C = A * B;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix / matrix" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N) + MatD::Ones(M,N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());
  Mat B = Mat::Copy(M,N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) / b(i,j);

  Mat C = A / B;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "matrix + scalar" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) + b;

  Mat C = A + b;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix - scalar" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) - b;

  Mat C = A - b;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix * scalar" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) * b;

  Mat C = A * b;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "matrix / scalar" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0) + 1.0;

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) / b;

  Mat C = A / b;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "scalar + matrix" )
{
  MatD b = MatD::Random(M,N);
  double a = b(0,0);

  Mat B = Mat::Copy(M,N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a + b(i,j);

  Mat C = a + B;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "scalar - matrix" )
{
  MatD b = MatD::Random(M,N);
  double a = b(0,0);

  Mat B = Mat::Copy(M,N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a - b(i,j);

  Mat C = a - B;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "scalar * matrix" )
{
  MatD b = MatD::Random(M,N);
  double a = b(0,0);

  Mat B = Mat::Copy(M,N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a * b(i,j);

  Mat C = a * B;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "scalar / matrix" )
{
  MatD b = MatD::Random(M,N) + MatD::Ones(M,N);
  double a = b(0,0);

  Mat B = Mat::Copy(M,N, b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a / b(i,j);

  Mat C = a / B;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( c(i,j), C(i,j) );
}

// =================================================================================================
// basic algebra
// =================================================================================================

SECTION( "minCoeff" )
{
  ColD a = ColD::Random(M*N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  Vec C = A.minCoeff(-1);

  Vec c = Vec::Constant(M, A.maxCoeff());

  for ( size_t i = 0 ; i < M ; i++ )
    for ( size_t j = 0 ; j < N ; j++ )
      c(i) = std::min( c(i), A(i,j) );

  REQUIRE( c.size() == C.size() );

  for ( size_t i = 0 ; i < c.size() ; ++i ) EQ( c[i], C[i] );
}

// -------------------------------------------------------------------------------------------------

SECTION( "maxCoeff" )
{
  ColD a = ColD::Random(M*N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  Vec C = A.maxCoeff(-1);

  Vec c = Vec::Constant(M, A.minCoeff());

  for ( size_t i = 0 ; i < M ; i++ )
    for ( size_t j = 0 ; j < N ; j++ )
      c(i) = std::max( c(i), A(i,j) );

  REQUIRE( c.size() == C.size() );

  for ( size_t i = 0 ; i < c.size() ; ++i ) EQ( c[i], C[i] );
}

// -------------------------------------------------------------------------------------------------

SECTION( "sum" )
{
  ColD a = ColD::Random(M*N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  Vec C = A.sum(0);

  Vec c = Vec::Zero(N);

  for ( size_t i = 0 ; i < M ; i++ )
    for ( size_t j = 0 ; j < N ; j++ )
      c(j) += A(i,j);

  REQUIRE( c.size() == C.size() );

  for ( size_t i = 0 ; i < c.size() ; ++i ) EQ( c[i], C[i] );
}

// -------------------------------------------------------------------------------------------------

SECTION( "mean" )
{
  ColD a = ColD::Random(M*N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  Vec C = A.mean(0);

  Vec c = Vec::Zero(N);

  for ( size_t i = 0 ; i < M ; i++ )
    for ( size_t j = 0 ; j < N ; j++ )
      c(j) += A(i,j);

  c /= static_cast<double>(M);

  REQUIRE( c.size() == C.size() );

  for ( size_t i = 0 ; i < c.size() ; ++i ) EQ( c[i], C[i] );
}

// -------------------------------------------------------------------------------------------------

SECTION( "average" )
{
  ColD a = ColD::Random(M*N);
  ColD b = ColD::Random(M*N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());
  Mat B = Mat::Copy(M,N, b.data(), b.data()+b.size());

  Vec C = A.average(B, 0);

  Vec c = Vec::Zero(N);
  Vec d = Vec::Zero(N);

  for ( size_t i = 0 ; i < M ; i++ ) {
    for ( size_t j = 0 ; j < N ; j++ ) {
      c(j) += B(i,j) * A(i,j);
      d(j) += B(i,j);
    }
  }

  c /= d;

  REQUIRE( c.size() == C.size() );

  for ( size_t i = 0 ; i < c.size() ; ++i ) EQ( c[i], C[i] );
}

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "minCoeff" )
{
  ColD a = ColD::Random(M*N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  double C = A.minCoeff();

  double c = a.minCoeff();

  EQ( c, C );
}

// -------------------------------------------------------------------------------------------------

SECTION( "maxCoeff" )
{
  ColD a = ColD::Random(M*N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  double C = A.maxCoeff();

  double c = a.maxCoeff();

  EQ( c, C );
}

// -------------------------------------------------------------------------------------------------

SECTION( "sum" )
{
  ColD a = ColD::Random(M*N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  double C = A.sum();

  double c = 0.0;

  for ( auto i = 0 ; i < a.size() ; ++i ) c += a[i];

  EQ( c, C );
}

// -------------------------------------------------------------------------------------------------

SECTION( "mean" )
{
  ColD a = ColD::Random(M*N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  double C = A.mean();

  double c = 0.0;

  for ( auto i = 0 ; i < a.size() ; ++i ) c += a[i];

  c /= static_cast<double>(M*N);

  EQ( c, C );
}

// -------------------------------------------------------------------------------------------------

SECTION( "average" )
{
  ColD a = ColD::Random(M*N);
  ColD b = ColD::Random(M*N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());
  Mat B = Mat::Copy(M,N, b.data(), b.data()+b.size());

  double C = A.average(B);

  double c = 0.0;
  double d = 0.0;

  for ( auto i = 0 ; i < a.size() ; ++i ) c += b[i] * a[i];
  for ( auto i = 0 ; i < a.size() ; ++i ) d += b[i];

  c /= d;

  EQ( c, C );
}

// =================================================================================================
// absolute value
// =================================================================================================

SECTION( "abs" )
{
  ColD a = ColD::Random(M*N) - ColD::Constant(M*N, .5);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  Mat C = cppmat::abs(A);

  ColD c = ColD(M*N);

  for ( auto i = 0 ; i < a.size() ; ++i ) c[i] = std::abs(a[i]);

  REQUIRE( c.size() == C.size() );

  for ( auto i = 0 ; i < c.size() ; ++i ) EQ( c[i], C[i] );
}

// -------------------------------------------------------------------------------------------------

SECTION( "abs" )
{
  ColD a = ColD::Random(M*N) - ColD::Constant(M*N, .5);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());
  Mat C = A;

  C.abs();

  ColD c = ColD(M*N);

  for ( auto i = 0 ; i < a.size() ; ++i ) c[i] = std::abs(a[i]);

  REQUIRE( c.size() == C.size() );

  for ( auto i = 0 ; i < c.size() ; ++i ) EQ( c[i], C[i] );
}

// =================================================================================================
// index operators
// =================================================================================================

SECTION( "at" )
{
  ColD a = ColD::Random(M*N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  std::vector<size_t> idx = {1,2};

  EQ( A.at(idx.begin(), idx.end()), A(1,2) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "decompress" )
{
  ColD a = ColD::Random(M*N);

  Mat A = Mat::Copy(M,N, a.data(), a.data()+a.size());

  std::vector<size_t> idx = A.decompress(A.compress(1,2));
  std::vector<size_t> jdx = {1,2};

  REQUIRE( idx == jdx );
}

// =================================================================================================

}
