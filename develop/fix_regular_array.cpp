
#include "support.h"

static const size_t M = 11;
static const size_t N = 9;

typedef cppmat::tiny::array<double,2,M,N> Arr;

// =================================================================================================

TEST_CASE("cppmat::tiny::array", "matrix.h")
{

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "array += array" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());
  Arr B = Arr::Copy(b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) += b(i,j);

  A += B;

  Equal(A, a);
}

// -------------------------------------------------------------------------------------------------

SECTION( "array -= array" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());
  Arr B = Arr::Copy(b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) -= b(i,j);

  A -= B;

  Equal(A, a);
}

// -------------------------------------------------------------------------------------------------

SECTION( "array *= array" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());
  Arr B = Arr::Copy(b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) *= b(i,j);

  A *= B;

  Equal(A, a);
}

// -------------------------------------------------------------------------------------------------

SECTION( "array /= array" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N) + MatD::Ones(M,N);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());
  Arr B = Arr::Copy(b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) /= b(i,j);

  A /= B;

  Equal(A, a);
}

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "array += scalar" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) += b;

  A += b;

  Equal(A, a);
}

// -------------------------------------------------------------------------------------------------

SECTION( "array -= scalar" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) -= b;

  A -= b;

  Equal(A, a);
}

// -------------------------------------------------------------------------------------------------

SECTION( "array *= scalar" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) *= b;

  A *= b;

  Equal(A, a);
}

// -------------------------------------------------------------------------------------------------

SECTION( "array /= scalar" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0) + 1.0;

  Arr A = Arr::Copy(a.data(), a.data()+a.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) /= b;

  A /= b;

  Equal(A, a);
}

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "array + array" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());
  Arr B = Arr::Copy(b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) + b(i,j);

  Arr C = A + B;

  Equal(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "array - array" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());
  Arr B = Arr::Copy(b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) - b(i,j);

  Arr C = A - B;

  Equal(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "array * array" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());
  Arr B = Arr::Copy(b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) * b(i,j);

  Arr C = A * B;

  Equal(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "array / array" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N) + MatD::Ones(M,N);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());
  Arr B = Arr::Copy(b.data(), b.data()+b.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) / b(i,j);

  Arr C = A / B;

  Equal(C, c);
}

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "array + scalar" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) + b;

  Arr C = A + b;

  Equal(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "array - scalar" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) - b;

  Arr C = A - b;

  Equal(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "array * scalar" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) * b;

  Arr C = A * b;

  Equal(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "array / scalar" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0) + 1.0;

  Arr A = Arr::Copy(a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = a(i,j) / b;

  Arr C = A / b;

  Equal(C, c);
}

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "scalar + array" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = b + a(i,j);

  Arr C = b + A;

  Equal(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "scalar - array" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = b - a(i,j);

  Arr C = b - A;

  Equal(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "scalar * array" )
{
  MatD a = MatD::Random(M,N);
  double b = a(0,0);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = b * a(i,j);

  Arr C = b * A;

  Equal(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "scalar / array" )
{
  MatD a = MatD::Random(M,N) + MatD::Ones(M,N);
  double b = a(0,0);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = b / a(i,j);

  Arr C = b / A;

  Equal(C, c);
}

// =================================================================================================
// algebra
// =================================================================================================

SECTION( "min" )
{
  MatD a = MatD::Random(M,N);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());

  double c = a.minCoeff();

  double C = A.min();

  EQ(c, C);
}

// -------------------------------------------------------------------------------------------------

SECTION( "max" )
{
  MatD a = MatD::Random(M,N);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());

  double c = a.maxCoeff();

  double C = A.max();

  EQ(c, C);
}

// -------------------------------------------------------------------------------------------------

SECTION( "sum" )
{
  MatD a = MatD::Random(M,N);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());

  double c = 0.0;

  for ( auto i = 0 ; i < a.rows() ; ++i )
    for ( auto j = 0 ; j < a.cols() ; ++j )
      c += a(i,j);

  double C = A.sum();

  EQ(c, C);
}

// -------------------------------------------------------------------------------------------------

SECTION( "mean" )
{
  MatD a = MatD::Random(M,N);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());

  double c = 0.0;

  for ( auto i = 0 ; i < a.rows() ; ++i )
    for ( auto j = 0 ; j < a.cols() ; ++j )
      c += a(i,j);

  c /= static_cast<double>(M*N);

  double C = A.mean();

  EQ(c, C);
}

// -------------------------------------------------------------------------------------------------

SECTION( "average" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());
  Arr B = Arr::Copy(b.data(), b.data()+b.size());

  double c = 0.0;
  double d = 0.0;

  for ( auto i = 0 ; i < a.rows() ; ++i ) {
    for ( auto j = 0 ; j < a.cols() ; ++j ) {
      c += b(i,j) * a(i,j);
      d += b(i,j);
    }
  }

  c /= d;

  double C = A.average(B);

  EQ(c, C);
}

// =================================================================================================
// absolute value
// =================================================================================================

SECTION( "abs" )
{
  MatD a = MatD::Random(M,N) - MatD::Constant(M,N,.5);

  Arr A = Arr::Copy(a.data(), a.data()+a.size());

  MatD c = a.cwiseAbs();

  Arr C = A.abs();

  Equal(C, c);
}

// =================================================================================================
// index operators
// =================================================================================================

SECTION( "at" )
{
  Arr A = Arr::Random();

  std::vector<size_t> idx = {1,2};

  EQ( A.at(idx.begin(), idx.end()), A(1,2) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "decompress" )
{
  Arr A = Arr::Random();

  std::vector<size_t> idx = A.decompress(A.compress(1,2));
  std::vector<size_t> jdx = {1,2};

  REQUIRE( idx == jdx );
}

// =================================================================================================

}
