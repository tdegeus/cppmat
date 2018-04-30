
#include <catch/catch.hpp>

#define EQ(a,b) REQUIRE_THAT( (a), Catch::WithinAbs((b), 1.e-10) );

#define CPPMAT_NOCONVERT
#include <cppmat.h>

#include <Eigen/Eigen>

typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor> ColD;

static const size_t N = 9;

typedef cppmat::periodic::vector<double> Vec;

// =================================================================================================

TEST_CASE("cppmat::periodic::vector", "vector.h")
{

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "vector += vector" )
{
  ColD a = ColD::Random(N);
  ColD b = ColD::Random(N);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());
  Vec B = Vec::Copy(b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < N ; ++i )
    a(i) += b(i);

  A += B;

  REQUIRE( static_cast<size_t>(a.size()) == A.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( a(i), A(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "vector -= vector" )
{
  ColD a = ColD::Random(N);
  ColD b = ColD::Random(N);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());
  Vec B = Vec::Copy(b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < N ; ++i )
    a(i) -= b(i);

  A -= B;

  REQUIRE( static_cast<size_t>(a.size()) == A.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( a(i), A(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "vector *= vector" )
{
  ColD a = ColD::Random(N);
  ColD b = ColD::Random(N);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());
  Vec B = Vec::Copy(b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < N ; ++i )
    a(i) *= b(i);

  A *= B;

  REQUIRE( static_cast<size_t>(a.size()) == A.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( a(i), A(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "vector /= vector" )
{
  ColD a = ColD::Random(N);
  ColD b = ColD::Random(N) + ColD::Ones(N);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());
  Vec B = Vec::Copy(b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < N ; ++i )
    a(i) /= b(i);

  A /= B;

  REQUIRE( static_cast<size_t>(a.size()) == A.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( a(i), A(i) );
}

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "vector += scalar" )
{
  ColD a = ColD::Random(N);
  double b = a(0,0);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());

  for ( size_t i = 0 ; i < N ; ++i )
    a(i) += b;

  A += b;

  REQUIRE( static_cast<size_t>(a.size()) == A.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( a(i), A(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "vector -= scalar" )
{
  ColD a = ColD::Random(N);
  double b = a(0,0);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());

  for ( size_t i = 0 ; i < N ; ++i )
    a(i) -= b;

  A -= b;

  REQUIRE( static_cast<size_t>(a.size()) == A.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( a(i), A(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "vector *= scalar" )
{
  ColD a = ColD::Random(N);
  double b = a(0,0);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());

  for ( size_t i = 0 ; i < N ; ++i )
    a(i) *= b;

  A *= b;

  REQUIRE( static_cast<size_t>(a.size()) == A.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( a(i), A(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "vector /= scalar" )
{
  ColD a = ColD::Random(N);
  double b = a(0,0) + 1.0;

  Vec A = Vec::Copy(a.data(), a.data()+a.size());

  for ( size_t i = 0 ; i < N ; ++i )
    a(i) /= b;

  A /= b;

  REQUIRE( static_cast<size_t>(a.size()) == A.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( a(i), A(i) );
}

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "vector + vector" )
{
  ColD a = ColD::Random(N);
  ColD b = ColD::Random(N);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());
  Vec B = Vec::Copy(b.data(), b.data()+b.size());

  ColD c = ColD::Zero(N);

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a(i) + b(i);

  Vec C = A + B;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( c(i), C(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "vector - vector" )
{
  ColD a = ColD::Random(N);
  ColD b = ColD::Random(N);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());
  Vec B = Vec::Copy(b.data(), b.data()+b.size());

  ColD c = ColD::Zero(N);

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a(i) - b(i);

  Vec C = A - B;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( c(i), C(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "vector * vector" )
{
  ColD a = ColD::Random(N);
  ColD b = ColD::Random(N);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());
  Vec B = Vec::Copy(b.data(), b.data()+b.size());

  ColD c = ColD::Zero(N);

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a(i) * b(i);

  Vec C = A * B;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( c(i), C(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "vector / vector" )
{
  ColD a = ColD::Random(N);
  ColD b = ColD::Random(N) + ColD::Ones(N);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());
  Vec B = Vec::Copy(b.data(), b.data()+b.size());

  ColD c = ColD::Zero(N);

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a(i) / b(i);

  Vec C = A / B;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( c(i), C(i) );
}

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "vector + scalar" )
{
  ColD a = ColD::Random(N);
  double b = a(0,0);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());

  ColD c = ColD::Zero(N);

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a(i) + b;

  Vec C = A + b;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( c(i), C(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "vector - scalar" )
{
  ColD a = ColD::Random(N);
  double b = a(0,0);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());

  ColD c = ColD::Zero(N);

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a(i) - b;

  Vec C = A - b;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( c(i), C(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "vector * scalar" )
{
  ColD a = ColD::Random(N);
  double b = a(0,0);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());

  ColD c = ColD::Zero(N);

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a(i) * b;

  Vec C = A * b;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( c(i), C(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "vector / scalar" )
{
  ColD a = ColD::Random(N);
  double b = a(0,0) + 1.0;

  Vec A = Vec::Copy(a.data(), a.data()+a.size());

  ColD c = ColD::Zero(N);

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a(i) / b;

  Vec C = A / b;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( c(i), C(i) );
}

// =================================================================================================
// arithmetic
// =================================================================================================

SECTION( "scalar + vector" )
{
  ColD b = ColD::Random(N);
  double a = b(0,0);

  Vec B = Vec::Copy(b.data(), b.data()+b.size());

  ColD c = ColD::Zero(N);

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a + b(i);

  Vec C = a + B;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( c(i), C(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "scalar - vector" )
{
  ColD b = ColD::Random(N);
  double a = b(0,0);

  Vec B = Vec::Copy(b.data(), b.data()+b.size());

  ColD c = ColD::Zero(N);

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a - b(i);

  Vec C = a - B;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( c(i), C(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "scalar * vector" )
{
  ColD b = ColD::Random(N);
  double a = b(0,0);

  Vec B = Vec::Copy(b.data(), b.data()+b.size());

  ColD c = ColD::Zero(N);

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a * b(i);

  Vec C = a * B;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( c(i), C(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "scalar / vector" )
{
  ColD b = ColD::Random(N) + ColD::Ones(N);
  double a = b(0,0);

  Vec B = Vec::Copy(b.data(), b.data()+b.size());

  ColD c = ColD::Zero(N);

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a / b(i);

  Vec C = a / B;

  REQUIRE( static_cast<size_t>(c.size()) == C.size() );

  for ( size_t i = 0 ; i < N ; ++i )
    EQ( c(i), C(i) );
}

// =================================================================================================
// basic algebra
// =================================================================================================

SECTION( "minCoeff" )
{
  ColD a = ColD::Random(N);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());

  double C = A.minCoeff();

  double c = a.minCoeff();

  EQ( c, C );
}

// -------------------------------------------------------------------------------------------------

SECTION( "maxCoeff" )
{
  ColD a = ColD::Random(N);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());

  double C = A.maxCoeff();

  double c = a.maxCoeff();

  EQ( c, C );
}

// -------------------------------------------------------------------------------------------------

SECTION( "sum" )
{
  ColD a = ColD::Random(N);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());

  double C = A.sum();

  double c = 0.0;

  for ( auto i = 0 ; i < a.size() ; ++i ) c += a[i];

  EQ( c, C );
}

// -------------------------------------------------------------------------------------------------

SECTION( "mean" )
{
  ColD a = ColD::Random(N);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());

  double C = A.mean();

  double c = 0.0;

  for ( auto i = 0 ; i < a.size() ; ++i ) c += a[i];

  c /= static_cast<double>(N);

  EQ( c, C );
}

// -------------------------------------------------------------------------------------------------

SECTION( "average" )
{
  ColD a = ColD::Random(N);
  ColD b = ColD::Random(N);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());
  Vec B = Vec::Copy(b.data(), b.data()+b.size());

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
  ColD a = ColD::Random(N) - ColD::Constant(N, .5);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());

  Vec C = cppmat::abs(A);

  ColD c = ColD(N);

  for ( auto i = 0 ; i < a.size() ; ++i ) c[i] = std::abs(a[i]);

  REQUIRE( c.size() == C.size() );

  for ( auto i = 0 ; i < c.size() ; ++i ) EQ( c[i], C[i] );
}

// -------------------------------------------------------------------------------------------------

SECTION( "abs" )
{
  ColD a = ColD::Random(N) - ColD::Constant(N, .5);

  Vec A = Vec::Copy(a.data(), a.data()+a.size());
  Vec C = A;

  C.abs();

  ColD c = ColD(N);

  for ( auto i = 0 ; i < a.size() ; ++i ) c[i] = std::abs(a[i]);

  REQUIRE( c.size() == C.size() );

  for ( auto i = 0 ; i < c.size() ; ++i ) EQ( c[i], C[i] );
}

// =================================================================================================

}
