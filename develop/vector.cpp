/* =================================================================================================

Compile using:

$ clang++ `pkg-config --cflags Eigen3 cppmat` -std=c++14 -pedantic -Wall -o test *.cpp
================================================================================================= */

#include <catch/catch.hpp>

#define CPPMAT_NOCONVERT
#include <cppmat/cppmat.h>

#include <Eigen/Eigen>

#define N 9

// =================================================================================================

TEST_CASE("cppmat::vector", "tiny_vector.h")
{

// =================================================================================================

SECTION( "vector + vector" )
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(N);
  Eigen::VectorXd b = Eigen::VectorXd::Random(N);
  Eigen::VectorXd c = a + b;

  // compute using cppmat

  cppmat::vector<double> A(N);
  cppmat::vector<double> B(N);

  for ( size_t i = 0 ; i < N ; ++i ) {
    A(i) = a(i);
    B(i) = b(i);
  }

  cppmat::vector<double> C = A + B;

  A += B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( C(i) - c(i) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( A(i) - c(i) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "vector - vector" )
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(N);
  Eigen::VectorXd b = Eigen::VectorXd::Random(N);
  Eigen::VectorXd c = a - b;

  // compute using cppmat

  cppmat::vector<double> A(N);
  cppmat::vector<double> B(N);

  for ( size_t i = 0 ; i < N ; ++i ) {
    A(i) = a(i);
    B(i) = b(i);
  }

  cppmat::vector<double> C = A - B;

  A -= B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( C(i) - c(i) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( A(i) - c(i) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "vector * vector" )
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(N);
  Eigen::VectorXd b = Eigen::VectorXd::Random(N);
  Eigen::VectorXd c(N);

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a(i) * b(i);

  // compute using cppmat

  cppmat::vector<double> A(N);
  cppmat::vector<double> B(N);

  for ( size_t i = 0 ; i < N ; ++i ) {
    A(i) = a(i);
    B(i) = b(i);
  }

  cppmat::vector<double> C = A * B;

  A *= B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( C(i) - c(i) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( A(i) - c(i) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "vector / vector" )
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(N);
  Eigen::VectorXd b = Eigen::VectorXd::Random(N);
  Eigen::VectorXd c(N);

  for ( size_t i = 0 ; i < N ; ++i )
    if ( b(i) == 0.0 )
        b(i) = 0.1;

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a(i) / b(i);

  // compute using cppmat

  cppmat::vector<double> A(N);
  cppmat::vector<double> B(N);

  for ( size_t i = 0 ; i < N ; ++i ) {
    A(i) = a(i);
    B(i) = b(i);
  }

  cppmat::vector<double> C = A / B;

  A /= B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( C(i) - c(i) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( A(i) - c(i) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION( "vector + scalar" )
{
  // compute using Eigen

  Eigen::VectorXd a  = Eigen::VectorXd::Random(N);
  Eigen::VectorXd bb = Eigen::VectorXd::Random(N);
  Eigen::VectorXd c(N);
  double b = bb(0,0);
  double B = b;

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a(i) + b;

  // compute using cppmat

  cppmat::vector<double> A(N);

  for ( size_t i = 0 ; i < N ; ++i )
    A(i) = a(i);

  cppmat::vector<double> C = A + B;

  A += B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( C(i) - c(i) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( A(i) - c(i) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "vector - scalar" )
{
  // compute using Eigen

  Eigen::VectorXd a  = Eigen::VectorXd::Random(N);
  Eigen::VectorXd bb = Eigen::VectorXd::Random(N);
  Eigen::VectorXd c(N);
  double b = bb(0,0);
  double B = b;

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a(i) - b;

  // compute using cppmat

  cppmat::vector<double> A(N);

  for ( size_t i = 0 ; i < N ; ++i )
    A(i) = a(i);

  cppmat::vector<double> C = A - B;

  A -= B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( C(i) - c(i) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( A(i) - c(i) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "vector * scalar" )
{
  // compute using Eigen

  Eigen::VectorXd a  = Eigen::VectorXd::Random(N);
  Eigen::VectorXd bb = Eigen::VectorXd::Random(N);
  Eigen::VectorXd c(N);
  double b = bb(0,0);
  double B = b;

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a(i) * b;

  // compute using cppmat

  cppmat::vector<double> A(N);

  for ( size_t i = 0 ; i < N ; ++i )
    A(i) = a(i);

  cppmat::vector<double> C = A * B;

  A *= B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( C(i) - c(i) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( A(i) - c(i) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "vector / scalar" )
{
  // compute using Eigen

  Eigen::VectorXd a  = Eigen::VectorXd::Random(N);
  Eigen::VectorXd bb = Eigen::VectorXd::Random(N);
  Eigen::VectorXd c(N);
  double b;

  for ( size_t i = 0 ; i < N ; ++i )
    if ( bb(i) != 0. )
        b = bb(i);

  double B = b;

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a(i) / b;

  // compute using cppmat

  cppmat::vector<double> A(N);

  for ( size_t i = 0 ; i < N ; ++i )
    A(i) = a(i);

  cppmat::vector<double> C = A / B;

  A /= B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( C(i) - c(i) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( A(i) - c(i) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION( "scalar + vector" )
{
  // compute using Eigen

  Eigen::VectorXd aa = Eigen::VectorXd::Random(N);
  Eigen::VectorXd b  = Eigen::VectorXd::Random(N);
  Eigen::VectorXd c(N);
  double a = aa(0,0);
  double A = a;

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a + b(i);

  // compute using cppmat

  cppmat::vector<double> B(N);

  for ( size_t i = 0 ; i < N ; ++i )
    B(i) = b(i);

  cppmat::vector<double> C = A + B;

  B += A;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( C(i) - c(i) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( B(i) - c(i) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "scalar - vector" )
{
  // compute using Eigen

  Eigen::VectorXd aa = Eigen::VectorXd::Random(N);
  Eigen::VectorXd b  = Eigen::VectorXd::Random(N);
  Eigen::VectorXd c(N);
  double a = aa(0,0);
  double A = a;

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a - b(i);

  // compute using cppmat

  cppmat::vector<double> B(N);

  for ( size_t i = 0 ; i < N ; ++i )
    B(i) = b(i);

  cppmat::vector<double> C = A - B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( C(i) - c(i) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "scalar * vector" )
{
  // compute using Eigen

  Eigen::VectorXd aa = Eigen::VectorXd::Random(N);
  Eigen::VectorXd b  = Eigen::VectorXd::Random(N);
  Eigen::VectorXd c(N);
  double a = aa(0,0);
  double A = a;

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a * b(i);

  // compute using cppmat

  cppmat::vector<double> B(N);

  for ( size_t i = 0 ; i < N ; ++i )
    B(i) = b(i);

  cppmat::vector<double> C = A * B;

  B *= A;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( C(i) - c(i) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( B(i) - c(i) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "scalar / vector" )
{
  // compute using Eigen

  Eigen::VectorXd aa = Eigen::VectorXd::Random(N);
  Eigen::VectorXd b  = Eigen::VectorXd::Random(N);
  Eigen::VectorXd c(N);
  double a = aa(0,0);
  double A = a;

  for ( size_t i = 0 ; i < N ; ++i )
    c(i) = a / b(i);

  // compute using cppmat

  cppmat::vector<double> B(N);

  for ( size_t i = 0 ; i < N ; ++i )
    B(i) = b(i);

  cppmat::vector<double> C = A / B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < N ; ++i )
    n += std::abs( C(i) - c(i) );

  REQUIRE( n < 1.e-12 );
}

}
