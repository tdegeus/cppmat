/* =================================================================================================

Compile using:

$ clang++ `pkg-config --cflags Eigen3 cppmat` -std=c++14 -Wpedantic -Wall -o test *.cpp
================================================================================================= */

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include <catch/catch.hpp>

#include <cppmat/matrix.h>
#include <Eigen/Eigen>

// =================================================================================================

TEST_CASE("cppmat::matrix", "matrix.h")
{

// =================================================================================================

SECTION( "matrix + matrix" )
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd c = a + b;

  // compute using cppmat

  cppmat::matrix<double> A({10,10});
  cppmat::matrix<double> B({10,10});

  for ( size_t i = 0 ; i < 10 ; ++i ) {
    for ( size_t j = 0 ; j < 10 ; ++j ) {
      A(i,j) = a(i,j);
      B(i,j) = b(i,j);
    }
  }

  cppmat::matrix<double> C = A + B;

  A += B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( A(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "matrix - matrix" )
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd c = a - b;

  // compute using cppmat

  cppmat::matrix<double> A({10,10});
  cppmat::matrix<double> B({10,10});

  for ( size_t i = 0 ; i < 10 ; ++i ) {
    for ( size_t j = 0 ; j < 10 ; ++j ) {
      A(i,j) = a(i,j);
      B(i,j) = b(i,j);
    }
  }

  cppmat::matrix<double> C = A - B;

  A -= B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( A(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "matrix * matrix" )
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd c(10,10);

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      c(i,j) = a(i,j) * b(i,j);

  // compute using cppmat

  cppmat::matrix<double> A({10,10});
  cppmat::matrix<double> B({10,10});

  for ( size_t i = 0 ; i < 10 ; ++i ) {
    for ( size_t j = 0 ; j < 10 ; ++j ) {
      A(i,j) = a(i,j);
      B(i,j) = b(i,j);
    }
  }

  cppmat::matrix<double> C = A * B;

  A *= B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( A(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "matrix / matrix" )
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd c(10,10);

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      if ( b(i,j) == 0.0 )
        b(i,j) = 0.1;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      c(i,j) = a(i,j) / b(i,j);

  // compute using cppmat

  cppmat::matrix<double> A({10,10});
  cppmat::matrix<double> B({10,10});

  for ( size_t i = 0 ; i < 10 ; ++i ) {
    for ( size_t j = 0 ; j < 10 ; ++j ) {
      A(i,j) = a(i,j);
      B(i,j) = b(i,j);
    }
  }

  cppmat::matrix<double> C = A / B;

  A /= B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( A(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION( "matrix + scalar" )
{
  // compute using Eigen

  Eigen::MatrixXd a  = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd bb = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd c(10,10);
  double b = bb(0,0);
  double B = b;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      c(i,j) = a(i,j) + b;

  // compute using cppmat

  cppmat::matrix<double> A({10,10});

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      A(i,j) = a(i,j);

  cppmat::matrix<double> C = A + B;

  A += B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( A(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "matrix - scalar" )
{
  // compute using Eigen

  Eigen::MatrixXd a  = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd bb = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd c(10,10);
  double b = bb(0,0);
  double B = b;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      c(i,j) = a(i,j) - b;

  // compute using cppmat

  cppmat::matrix<double> A({10,10});

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      A(i,j) = a(i,j);

  cppmat::matrix<double> C = A - B;

  A -= B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( A(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "matrix * scalar" )
{
  // compute using Eigen

  Eigen::MatrixXd a  = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd bb = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd c(10,10);
  double b = bb(0,0);
  double B = b;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      c(i,j) = a(i,j) * b;

  // compute using cppmat

  cppmat::matrix<double> A({10,10});

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      A(i,j) = a(i,j);

  cppmat::matrix<double> C = A * B;

  A *= B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( A(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "matrix / scalar" )
{
  // compute using Eigen

  Eigen::MatrixXd a  = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd bb = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd c(10,10);
  double b;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      if ( bb(i,j) != 0. )
        b = bb(i,j);

  double B = b;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      c(i,j) = a(i,j) / b;

  // compute using cppmat

  cppmat::matrix<double> A({10,10});

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      A(i,j) = a(i,j);

  cppmat::matrix<double> C = A / B;

  A /= B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( A(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION( "scalar + matrix" )
{
  // compute using Eigen

  Eigen::MatrixXd aa = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd b  = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd c(10,10);
  double a = aa(0,0);
  double A = a;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      c(i,j) = a + b(i,j);

  // compute using cppmat

  cppmat::matrix<double> B({10,10});

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      B(i,j) = b(i,j);

  cppmat::matrix<double> C = A + B;

  B += A;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( B(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "scalar - matrix" )
{
  // compute using Eigen

  Eigen::MatrixXd aa = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd b  = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd c(10,10);
  double a = aa(0,0);
  double A = a;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      c(i,j) = a - b(i,j);

  // compute using cppmat

  cppmat::matrix<double> B({10,10});

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      B(i,j) = b(i,j);

  cppmat::matrix<double> C = A - B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "scalar * matrix" )
{
  // compute using Eigen

  Eigen::MatrixXd aa = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd b  = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd c(10,10);
  double a = aa(0,0);
  double A = a;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      c(i,j) = a * b(i,j);

  // compute using cppmat

  cppmat::matrix<double> B({10,10});

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      B(i,j) = b(i,j);

  cppmat::matrix<double> C = A * B;

  B *= A;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( B(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );
}

// -----------------------------------------------------------------------------------------------

SECTION( "scalar / matrix" )
{
  // compute using Eigen

  Eigen::MatrixXd aa = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd b  = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd c(10,10);
  double a = aa(0,0);
  double A = a;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      c(i,j) = a / b(i,j);

  // compute using cppmat

  cppmat::matrix<double> B({10,10});

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      B(i,j) = b(i,j);

  cppmat::matrix<double> C = A / B;

  // verify

  double n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  REQUIRE( n < 1.e-12 );
}

}
