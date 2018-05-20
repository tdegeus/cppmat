
#include "support.h"

static const size_t M = 11;
static const size_t N = 9;
static const size_t O = 6;
static const size_t P = 5;

typedef cppmat::array<double> Arr;

// =================================================================================================

TEST_CASE("cppmat::array", "matrix.h")
{

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "array += array" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N);

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
  Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
  Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
  Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
  Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
  Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
  Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
  Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
  Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

  MatD c = MatD::Zero(M,N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      c(i,j) = b / a(i,j);

  Arr C = b / A;

  Equal(C, c);
}

// =================================================================================================
// algebra - partial
// =================================================================================================

SECTION( "min" )
{
  Arr A = Arr::Random({M,N,O,P});

  Arr C = A.min({-1,-2});

  Arr c = Arr::Constant({M,N}, A.max());

  for ( size_t i = 0 ; i < M ; i++ )
    for ( size_t j = 0 ; j < N ; j++ )
      for ( size_t k = 0 ; k < O ; k++ )
        for ( size_t l = 0 ; l < P ; l++ )
          c(i,j) = std::min( c(i,j), A(i,j,k,l) );

  Equal(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "max" )
{
  Arr A = Arr::Random({M,N,O,P});

  Arr C = A.max({-1,-2});

  Arr c = Arr::Constant({M,N}, A.min());

  for ( size_t i = 0 ; i < M ; i++ )
    for ( size_t j = 0 ; j < N ; j++ )
      for ( size_t k = 0 ; k < O ; k++ )
        for ( size_t l = 0 ; l < P ; l++ )
          c(i,j) = std::max( c(i,j), A(i,j,k,l) );

  Equal(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "sum" )
{
  Arr A = Arr::Random({M,N,O,P});

  Arr C = A.sum({0,1});

  Arr c = Arr::Zero({O,P});

  for ( size_t i = 0 ; i < M ; i++ )
    for ( size_t j = 0 ; j < N ; j++ )
      for ( size_t k = 0 ; k < O ; k++ )
        for ( size_t l = 0 ; l < P ; l++ )
          c(k,l) += A(i,j,k,l);

  Equal(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "mean" )
{
  Arr A = Arr::Random({M,N,O,P});

  Arr C = A.mean({0,1});

  Arr c = Arr::Zero({O,P});

  for ( size_t i = 0 ; i < M ; i++ )
    for ( size_t j = 0 ; j < N ; j++ )
      for ( size_t k = 0 ; k < O ; k++ )
        for ( size_t l = 0 ; l < P ; l++ )
          c(k,l) += A(i,j,k,l);

  c /= static_cast<double>(M*N);

  Equal(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "average" )
{
  Arr A = Arr::Random({M,N,O,P});
  Arr B = Arr::Random({M,N,O,P});

  Arr C = A.average(B,{-2,0});

  Arr c = Arr::Zero({N,P});
  Arr d = Arr::Zero({N,P});

  for ( size_t i = 0 ; i < M ; i++ ) {
    for ( size_t j = 0 ; j < N ; j++ ) {
      for ( size_t k = 0 ; k < O ; k++ ) {
        for ( size_t l = 0 ; l < P ; l++ ) {
          c(j,l) += B(i,j,k,l) * A(i,j,k,l);
          d(j,l) += B(i,j,k,l);
        }
      }
    }
  }

  c /= d;

  Equal(C, c);
}

// =================================================================================================
// algebra
// =================================================================================================

SECTION( "min" )
{
  MatD a = MatD::Random(M,N);

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

  double c = a.minCoeff();

  double C = A.min();

  EQ(c, C);
}

// -------------------------------------------------------------------------------------------------

SECTION( "max" )
{
  MatD a = MatD::Random(M,N);

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

  double c = a.maxCoeff();

  double C = A.max();

  EQ(c, C);
}

// -------------------------------------------------------------------------------------------------

SECTION( "sum" )
{
  MatD a = MatD::Random(M,N);

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
  Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

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

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

  MatD c = a.cwiseAbs();

  Arr C = A.abs();

  Equal(C, c);
}

// =================================================================================================
// index operators
// =================================================================================================

SECTION( "at" )
{
  Arr A = Arr::Random({M,N,O,P});

  std::vector<size_t> idx = {1,2,3,4};

  EQ( A.at(idx.begin(), idx.end()), A(1,2,3,4) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "decompress" )
{
  ColD a = ColD::Random(M*N*O*P);

  Arr A = Arr::Copy({M,N,O,P}, a.data(), a.data()+a.size());

  std::vector<size_t> idx = A.decompress(A.compress(1,2,3,4));
  std::vector<size_t> jdx = {1,2,3,4};

  REQUIRE( idx == jdx );
}

// =================================================================================================

}
