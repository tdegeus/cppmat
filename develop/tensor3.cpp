/* =================================================================================================

Compile using:

$ clang++ `pkg-config --cflags Eigen3 cppmat` -std=c++14 -pedantic -Wall -o test *.cpp
================================================================================================= */

#include <catch/catch.hpp>

#define CPPMAT_NOCONVERT
#include <cppmat/cppmat.h>

#include <Eigen/Eigen>

// =================================================================================================

TEST_CASE("cppmat::cartesian3d", "tensor3.h")
{

using     T4  = cppmat::cartesian3d::tensor4<double>;
using     T2  = cppmat::cartesian3d::tensor2<double>;
using     T2d = cppmat::cartesian3d::tensor2d<double>;
using     T2s = cppmat::cartesian3d::tensor2s<double>;
using     V   = cppmat::cartesian3d::vector<double>;
namespace cm  = cppmat::cartesian3d;

size_t nd = 3;
double n;

// =================================================================================================

SECTION("tensor4 arithmetic")
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::VectorXd b = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::VectorXd c = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::VectorXd d = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::VectorXd e = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::VectorXd f(nd*nd*nd*nd);

  for ( size_t i=0 ; i<nd*nd*nd*nd ; ++i )
    f(i) = ( ( 10.*a(i) ) * ( b(i)/3. ) ) / ( 0.5*c(i) ) - 5./d(i) + 2.*e(i);

  // compute using cppmat

  T4 A(nd),B(nd),C(nd),D(nd),E(nd),F(nd);

  for ( size_t i=0; i<nd; ++i ) {
    for ( size_t j=0; j<nd; ++j ) {
      for ( size_t k=0; k<nd; ++k ) {
        for ( size_t l=0; l<nd; ++l ) {
          A(i,j,k,l) = a(i*nd*nd*nd+j*nd*nd+k*nd+l);
          B(i,j,k,l) = b(i*nd*nd*nd+j*nd*nd+k*nd+l);
          C(i,j,k,l) = c(i*nd*nd*nd+j*nd*nd+k*nd+l);
          D(i,j,k,l) = d(i*nd*nd*nd+j*nd*nd+k*nd+l);
          E(i,j,k,l) = e(i*nd*nd*nd+j*nd*nd+k*nd+l);
        }
      }
    }
  }

  F = ( ( A*10. ) * ( B/3. ) ) / ( 0.5*C ) - 5./D + 2.*E;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( f(i*nd*nd*nd+j*nd*nd+k*nd+l)-F(i,j,k,l) );

  REQUIRE( n < 1.e-12 );

  // compute using cppmat

  F.zeros();

  A *= 10.;
  F += A;

  B /= 3.;
  F *= B;

  C *= 0.5;
  F /= C;

  D  = 5./D;
  F -= D;

  E *= 2.;
  F += E;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( f(i*nd*nd*nd+j*nd*nd+k*nd+l)-F(i,j,k,l) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor4.ddot( tensor4 )")
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::VectorXd b = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::VectorXd c(nd*nd*nd*nd);

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          for ( size_t m=0; m<nd; ++m )
            for ( size_t n=0; n<nd; ++n )
              c(i*nd*nd*nd+j*nd*nd+m*nd+n) += a(i*nd*nd*nd+j*nd*nd+k*nd+l)*b(l*nd*nd*nd+k*nd*nd+m*nd+n);

  // compute using cppmat

  T4 A(nd),B(nd);

  for ( size_t i=0; i<nd; ++i ) {
    for ( size_t j=0; j<nd; ++j ) {
      for ( size_t k=0; k<nd; ++k ) {
        for ( size_t l=0; l<nd; ++l ) {
          A(i,j,k,l) = a(i*nd*nd*nd+j*nd*nd+k*nd+l);
          B(i,j,k,l) = b(i*nd*nd*nd+j*nd*nd+k*nd+l);
        }
      }
    }
  }

  T4 C = A.ddot( B );
  T4 D = cm::ddot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-C(i,j,k,l) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor4.ddot( tensor4.T() )")
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::VectorXd b = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::VectorXd c(nd*nd*nd*nd),d(nd*nd*nd*nd);

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          d(l*nd*nd*nd+k*nd*nd+j*nd+i) = b(i*nd*nd*nd+j*nd*nd+k*nd+l);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          for ( size_t m=0; m<nd; ++m )
            for ( size_t n=0; n<nd; ++n )
              c(i*nd*nd*nd+j*nd*nd+m*nd+n) += a(i*nd*nd*nd+j*nd*nd+k*nd+l)*d(l*nd*nd*nd+k*nd*nd+m*nd+n);

  // compute using cppmat

  T4 A(nd),B(nd);

  for ( size_t i=0; i<nd; ++i ) {
    for ( size_t j=0; j<nd; ++j ) {
      for ( size_t k=0; k<nd; ++k ) {
        for ( size_t l=0; l<nd; ++l ) {
          A(i,j,k,l) = a(i*nd*nd*nd+j*nd*nd+k*nd+l);
          B(i,j,k,l) = b(i*nd*nd*nd+j*nd*nd+k*nd+l);
        }
      }
    }
  }

  T4 C = A.ddot( B.T() );
  T4 D = cm::ddot( A , cm::transpose(B) );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-C(i,j,k,l) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor4.ddot( tensor4.LT() )")
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::VectorXd b = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::VectorXd c(nd*nd*nd*nd),d(nd*nd*nd*nd);

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          d(j*nd*nd*nd+i*nd*nd+k*nd+l) = b(i*nd*nd*nd+j*nd*nd+k*nd+l);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          for ( size_t m=0; m<nd; ++m )
            for ( size_t n=0; n<nd; ++n )
              c(i*nd*nd*nd+j*nd*nd+m*nd+n) += a(i*nd*nd*nd+j*nd*nd+k*nd+l)*d(l*nd*nd*nd+k*nd*nd+m*nd+n);

  // compute using cppmat

  T4 A(nd),B(nd);

  for ( size_t i=0; i<nd; ++i ) {
    for ( size_t j=0; j<nd; ++j ) {
      for ( size_t k=0; k<nd; ++k ) {
        for ( size_t l=0; l<nd; ++l ) {
          A(i,j,k,l) = a(i*nd*nd*nd+j*nd*nd+k*nd+l);
          B(i,j,k,l) = b(i*nd*nd*nd+j*nd*nd+k*nd+l);
        }
      }
    }
  }

  T4 C = A.ddot( B.LT() );
  T4 D = cm::ddot( A , cm::transposeL(B) );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-C(i,j,k,l) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor4.ddot( tensor4.RT() )")
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::VectorXd b = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::VectorXd c(nd*nd*nd*nd),d(nd*nd*nd*nd);

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          d(i*nd*nd*nd+j*nd*nd+l*nd+k) = b(i*nd*nd*nd+j*nd*nd+k*nd+l);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          for ( size_t m=0; m<nd; ++m )
            for ( size_t n=0; n<nd; ++n )
              c(i*nd*nd*nd+j*nd*nd+m*nd+n) += a(i*nd*nd*nd+j*nd*nd+k*nd+l)*d(l*nd*nd*nd+k*nd*nd+m*nd+n);

  // compute using cppmat

  T4 A(nd),B(nd);

  for ( size_t i=0; i<nd; ++i ) {
    for ( size_t j=0; j<nd; ++j ) {
      for ( size_t k=0; k<nd; ++k ) {
        for ( size_t l=0; l<nd; ++l ) {
          A(i,j,k,l) = a(i*nd*nd*nd+j*nd*nd+k*nd+l);
          B(i,j,k,l) = b(i*nd*nd*nd+j*nd*nd+k*nd+l);
        }
      }
    }
  }

  T4 C = A.ddot( B.RT() );
  T4 D = cm::ddot( A , cm::transposeR(B) );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-C(i,j,k,l) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor4.ddot( tensor2 )")
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          c(i,j) += a(i*nd*nd*nd+j*nd*nd+k*nd+l)*b(l,k);

  // compute using cppmat

  T4 A(nd);
  cm::tensor2<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          A(i,j,k,l) = a(i*nd*nd*nd+j*nd*nd+k*nd+l);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2<double> C = A.ddot( B );
  cm::tensor2<double> D = cm::ddot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor4.ddot( tensor2s )")
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          c(i,j) += a(i*nd*nd*nd+j*nd*nd+k*nd+l)*b(l,k);

  // compute using cppmat

  T4 A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          A(i,j,k,l) = a(i*nd*nd*nd+j*nd*nd+k*nd+l);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2<double> C = A.ddot( B );
  cm::tensor2<double> D = cm::ddot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor4.ddot( tensor2d )")
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          c(i,j) += a(i*nd*nd*nd+j*nd*nd+k*nd+l)*b(l,k);

  // compute using cppmat

  T4 A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          A(i,j,k,l) = a(i*nd*nd*nd+j*nd*nd+k*nd+l);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2<double> C = A.ddot( B );
  cm::tensor2<double> D = cm::ddot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor2.ddot( tensor4 )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd b = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::MatrixXd c(nd,nd);

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          c(k,l) += a(i,j)*b(j*nd*nd*nd+i*nd*nd+k*nd+l);

  // compute using cppmat

  cm::tensor2<double> A(nd);
  T4 B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          B(i,j,k,l) = b(i*nd*nd*nd+j*nd*nd+k*nd+l);

  cm::tensor2<double> C = A.ddot( B );
  cm::tensor2<double> D = cm::ddot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor2.ddot( tensor2 )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  double c = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c += a(i,j)*b(j,i);

  // compute using cppmat

  cm::tensor2<double> A(nd), B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  double C = A.ddot( B );
  double D = cm::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = std::abs( D - c );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2.ddot( tensor2s )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  double c = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c += a(i,j)*b(j,i);

  // compute using cppmat

  T2  A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  double C = A.ddot( B );
  double D = cm::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = std::abs( D - c );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2.ddot( tensor2d )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  double c = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c += a(i,j)*b(j,i);

  // compute using cppmat

  T2  A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  double C = A.ddot( B );
  double D = cm::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = std::abs( D - c );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2.dot( tensor2 )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        c(i,k) += a(i,j)*b(j,k);

  // compute using cppmat

  T2  A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  T2  C = A.dot( B );
  T2  D = cm::dot ( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor2.dot( tensor2s )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        c(i,k) += a(i,j)*b(j,k);

  // compute using cppmat

  T2  A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  T2  C = A.dot( B );
  T2  D = cm::dot ( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor2.dot( tensor2d )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        c(i,k) += a(i,j)*b(j,k);

  // compute using cppmat

  T2  A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  T2  C = A.dot( B );
  T2  D = cm::dot ( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor2.dot( vector )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd b = Eigen::VectorXd::Random(nd);
  Eigen::VectorXd c(nd);

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i) += a(i,j)*b(j);

  // compute using cppmat

  T2  A(nd);
  V B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i) = b(i);

  cm::vector <double> C = A.dot( B );
  cm::vector <double> D = cm::dot( A, B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-C(i) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-D(i) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2.dyadic( tensor2 )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd c(nd*nd*nd*nd);

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          c(i*nd*nd*nd+j*nd*nd+k*nd+l) += a(i,j)*b(k,l);

  // compute using cppmat

  T2  A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  T4 C = A.dyadic( B );
  T4 D = cm::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2.dyadic( tensor2s )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd c(nd*nd*nd*nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          c(i*nd*nd*nd+j*nd*nd+k*nd+l) += a(i,j)*b(k,l);

  // compute using cppmat

  T2  A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  T4 C = A.dyadic( B );
  T4 D = cm::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2.dyadic( tensor2d )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd c(nd*nd*nd*nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          c(i*nd*nd*nd+j*nd*nd+k*nd+l) += a(i,j)*b(k,l);

  // compute using cppmat

  T2  A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  T4 C = A.dyadic( B );
  T4 D = cm::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2 arithmetic")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd d = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd e = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd f(nd,nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      f(i,j) = ( 7.*a(i,j) + 10./b(i,j) ) / ( c(i,j)+2. ) * ( d(i,j)-1. ) - 2.*e(i,j);

  // compute using cppmat

  T2  A(nd);
  T2  B(nd);
  T2  C(nd);
  T2  D(nd);
  T2  E(nd);

  for ( size_t i=0; i<nd; ++i ) {
    for ( size_t j=0; j<nd; ++j ) {
      A(i,j) = a(i,j);
      B(i,j) = b(i,j);
      C(i,j) = c(i,j);
      D(i,j) = d(i,j);
      E(i,j) = e(i,j);
    }
  }

  T2  F = ( 7.*A + 10./B ) / ( C+2. ) * ( D-1. ) - 2.*E;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( f(i,j)-F(i,j) );

  REQUIRE( n < 1.e-12 );

  // compute using cppmat

  F.zeros();

  A *= 7.;
  F += A;

  B  = 10./B;
  F += B;

  C += 2.;
  F /= C;

  D -= 1.;
  F *= D;

  E *= 2.;
  F -= E;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( f(i,j)-F(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2.dot( tensor2.T() )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd d(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      d(j,i) = b(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        c(i,k) += a(i,j)*d(j,k);

  // compute using cppmat

  T2  A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  T2  C = A.dot( B.T() );
  T2  D = cm::dot ( A , cm::transpose( B ) );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor2.trace()")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  double c = a.trace();

  // compute using cppmat

  T2  A(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  double C = A.trace();
  double D = cm::trace ( A );

  // check the result

  n = std::abs( C - c );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = std::abs( D - c );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2.det() -- 3D")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3,3);
  double c = a.determinant();

  // compute using cppmat

  T2  A(3);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      A(i,j) = a(i,j);

  double C = A.det();
  double D = cm::det ( A );

  // check the result

  n = std::abs( C - c );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = std::abs( D - c );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2.inv() -- 3D")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3,3);
  Eigen::MatrixXd c = a.inverse();

  // compute using cppmat

  T2  A(3);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      A(i,j) = a(i,j);

  T2  C = A.inv();
  T2  D = cm::inv( A );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor2s.ddot( tensor4 )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd b = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          c(k,l) += a(i,j)*b(j*nd*nd*nd+i*nd*nd+k*nd+l);

  // compute using cppmat

  T2s A(nd);
  T4 B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          B(i,j,k,l) = b(i*nd*nd*nd+j*nd*nd+k*nd+l);

  T2  C = A.ddot( B );
  T2  D = cm::ddot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor2s.ddot( tensor2 )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  double c = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c += a(i,j)*b(j,i);

  // compute using cppmat

  T2s A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  double C = A.ddot( B );
  double D = cm::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = std::abs( D - c );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s.ddot( tensor2s )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  double c = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c += a(i,j)*b(j,i);

  // compute using cppmat

  T2s A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  double C = A.ddot( B );
  double D = cm::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = std::abs( D - c );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s.ddot( tensor2d )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  double c = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c += a(i,j)*b(j,i);

  // compute using cppmat

  T2s A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  double C = A.ddot( B );
  double D = cm::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = std::abs( D - c );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s.dot( tensor2 )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        c(i,k) += a(i,j)*b(j,k);

  // compute using cppmat

  T2s A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  T2  C = A.dot( B );
  T2  D = cm::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor2s.dot( tensor2s )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        c(i,k) += a(i,j)*b(j,k);

  // compute using cppmat

  T2s A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  T2  C = A.dot( B );
  T2  D = cm::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor2s.dot( tensor2d )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        c(i,k) += a(i,j)*b(j,k);

  // compute using cppmat

  T2s A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  T2  C = A.dot( B );
  T2  D = cm::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor2s.dot( vector )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd b = Eigen::VectorXd::Random(nd);
  Eigen::VectorXd c(nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i) += a(i,j)*b(j);

  // compute using cppmat

  T2s A(nd);
  V B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i) = b(i);

  cm::vector <double> C = A.dot( B );
  cm::vector <double> D = cm::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-C(i) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-D(i) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s.dyadic( tensor2 )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd c(nd*nd*nd*nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          c(i*nd*nd*nd+j*nd*nd+k*nd+l) += a(i,j)*b(k,l);

  // compute using cppmat

  T2s A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  T4 C = A.dyadic( B );
  T4 D = cm::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s.dyadic( tensor2s )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd c(nd*nd*nd*nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          c(i*nd*nd*nd+j*nd*nd+k*nd+l) += a(i,j)*b(k,l);

  // compute using cppmat

  T2s A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  T4 C = A.dyadic( B );
  T4 D = cm::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s.dyadic( tensor2d )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd c(nd*nd*nd*nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          c(i*nd*nd*nd+j*nd*nd+k*nd+l) += a(i,j)*b(k,l);

  // compute using cppmat

  T2s A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  T4 C = A.dyadic( B );
  T4 D = cm::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s arithmetic")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd d = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd e = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd f(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      a(j,i) = a(i,j);
      b(j,i) = b(i,j);
      c(j,i) = c(i,j);
      d(j,i) = d(i,j);
      e(j,i) = e(i,j);
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      f(i,j) = ( 7.*a(i,j) + 10./b(i,j) ) / ( c(i,j)+2. ) * ( d(i,j)-1. ) - 2.*e(i,j);

  // compute using cppmat

  T2s A(nd);
  T2s B(nd);
  T2s C(nd);
  T2s D(nd);
  T2s E(nd);

  for ( size_t i=0; i<nd; ++i ) {
    for ( size_t j=i; j<nd; ++j ) {
      A(i,j) = a(i,j);
      B(i,j) = b(i,j);
      C(i,j) = c(i,j);
      D(i,j) = d(i,j);
      E(i,j) = e(i,j);
    }
  }

  T2s F = ( 7.*A + 10./B ) / ( C+2. ) * ( D-1. ) - 2.*E;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( f(i,j)-F(i,j) );

  REQUIRE( n < 1.e-12 );

  // compute using cppmat

  F.zeros();

  A *= 7.;
  F += A;

  B  = 10./B;
  F += B;

  C += 2.;
  F /= C;

  D -= 1.;
  F *= D;

  E *= 2.;
  F -= E;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( f(i,j)-F(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s.dot( tensor2s.T() )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd d(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      d(j,i) = b(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        c(i,k) += a(i,j)*d(j,k);

  // compute using cppmat

  T2s A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  T2  C = A.dot( B.T() );
  T2  D = cm::dot( A , cm::transpose( B ) );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor2s.trace()")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  double c = a.trace();

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  // compute using cppmat

  T2s A(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  double C = A.trace();
  double D = cm::trace( A );

  // check the result

  n = std::abs( C - c );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = std::abs( D - c );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s.det() -- 3D")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3,3);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = i+1 ; j < 3 ; ++j )
      a(j,i) = a(i,j);

  double c = a.determinant();

  // compute using cppmat

  T2s A(3);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=i; j<3; ++j )
      A(i,j) = a(i,j);

  double C = A.det();
  double D = cm::det( A );

  // check the result

  n = std::abs( C - c );

  REQUIRE( n < 1.e-10 );

  // check the result

  n = std::abs( D - c );

  REQUIRE( n < 1.e-10 );

}

// =================================================================================================

SECTION("tensor2s.inv() -- 3D")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3,3);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = i+1 ; j < 3 ; ++j )
      a(j,i) = a(i,j);

  Eigen::MatrixXd c = a.inverse();

  // compute using cppmat

  T2s A(3);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=i; j<3; ++j )
      A(i,j) = a(i,j);

  T2s C = A.inv();
  T2s D = cm::inv( A );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-10 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor2d.ddot( tensor4 )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd b = Eigen::VectorXd::Random(nd*nd*nd*nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          c(k,l) += a(i,j)*b(j*nd*nd*nd+i*nd*nd+k*nd+l);

  // compute using cppmat

  T2d A(nd);
  T4 B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          B(i,j,k,l) = b(i*nd*nd*nd+j*nd*nd+k*nd+l);

  T2  C = A.ddot( B );
  T2  D = cm::ddot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor2d.ddot( tensor2 )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  double c = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c += a(i,j)*b(j,i);

  // compute using cppmat

  T2d A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  double C = A.ddot( B );
  double D = cm::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = std::abs( D - c );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d.ddot( tensor2s )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  double c = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c += a(i,j)*b(j,i);

  // compute using cppmat

  T2d A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  double C = A.ddot( B );
  double D = cm::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = std::abs( D - c );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d.ddot( tensor2d )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  double c = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c += a(i,j)*b(j,i);

  // compute using cppmat

  T2d A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  double C = A.ddot( B );
  double D = cm::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = std::abs( D - c );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d.dot( tensor2 )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        c(i,k) += a(i,j)*b(j,k);

  // compute using cppmat

  T2d A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  T2  C = A.dot( B );
  T2  D = cm::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor2d.dot( tensor2s )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        c(i,k) += a(i,j)*b(j,k);

  // compute using cppmat

  T2d A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  T2  C = A.dot( B );
  T2  D = cm::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor2d.dot( tensor2d )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        c(i,k) += a(i,j)*b(j,k);

  // compute using cppmat

  T2d A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  T2d C = A.dot( B );
  T2d D = cm::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d.dot( vector )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd b = Eigen::VectorXd::Random(nd);
  Eigen::VectorXd c(nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i) += a(i,j)*b(j);

  // compute using cppmat

  T2d A(nd);
  V B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    B(i) = b(i);

  cm::vector <double> C = A.dot( B );
  cm::vector <double> D = cm::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-C(i) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-D(i) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d.dyadic( tensor2 )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd c(nd*nd*nd*nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          c(i*nd*nd*nd+j*nd*nd+k*nd+l) += a(i,j)*b(k,l);

  // compute using cppmat

  T2d A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  T4 C = A.dyadic( B );
  T4 D = cm::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d.dyadic( tensor2s )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd c(nd*nd*nd*nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          c(i*nd*nd*nd+j*nd*nd+k*nd+l) += a(i,j)*b(k,l);

  // compute using cppmat

  T2d A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  T4 C = A.dyadic( B );
  T4 D = cm::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d.dyadic( tensor2d )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd c(nd*nd*nd*nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          c(i*nd*nd*nd+j*nd*nd+k*nd+l) += a(i,j)*b(k,l);

  // compute using cppmat

  T2d A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  T4 C = A.dyadic( B );
  T4 D = cm::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d arithmetic")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd d = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd e = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd f(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = 0 ; j < nd ; ++j ) {
      if ( i != j ) {
        a(i,j) = 0.0;
        b(i,j) = 0.0;
        c(i,j) = 0.0;
        d(i,j) = 0.0;
        e(i,j) = 0.0;
      }
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      f(i,j) = ( 7.*a(i,j) + 10.*b(i,j) ) * ( c(i,j) )  * ( d(i,j) ) - 2.*e(i,j);

  // compute using cppmat

  T2d A(nd);
  T2d B(nd);
  T2d C(nd);
  T2d D(nd);
  T2d E(nd);

  for ( size_t i=0; i<nd; ++i ) {
    A(i,i) = a(i,i);
    B(i,i) = b(i,i);
    C(i,i) = c(i,i);
    D(i,i) = d(i,i);
    E(i,i) = e(i,i);
  }

  T2d F = ( 7.*A + 10.*B ) * ( C )  * ( D ) - 2.*E;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( f(i,j)-F(i,j) );

  REQUIRE( n < 1.e-12 );

  // compute using cppmat

  F.zeros();

  A *= 7.;
  F += A;

  B  = 10.*B;
  F += B;

  F *= C;

  F *= D;

  E *= 2.;
  F -= E;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( f(i,j)-F(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d.dot( tensor2s.T() )")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd d(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      d(j,i) = b(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        c(i,k) += a(i,j)*d(j,k);

  // compute using cppmat

  T2d A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  T2  C = A.dot( B.T() );
  T2  D = cm::dot( A , cm::transpose( B ) );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("tensor2d.trace()")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  double c = a.trace();

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  // compute using cppmat

  T2d A(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  double C = A.trace();
  double D = cm::trace( A );

  // check the result

  n = std::abs( C - c );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = std::abs( D - c );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d.det() -- 3D")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3,3);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  double c = a.determinant();

  // compute using cppmat

  T2d A(3);

  for ( size_t i=0; i<3; ++i )
    A(i,i) = a(i,i);

  double C = A.det();
  double D = cm::det( A );

  // check the result

  n = std::abs( C - c );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = std::abs( D - c );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d.inv() -- 3D")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3,3);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  Eigen::MatrixXd c = a.inverse();

  // compute using cppmat

  T2d A(3);

  for ( size_t i=0; i<3; ++i )
    A(i,i) = a(i,i);

  T2d C = A.inv();
  T2d D = cm::inv( A );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-10 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("vector.dot( vector )")
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(nd);
  Eigen::VectorXd b = Eigen::VectorXd::Random(nd);
  double c = 0.;

  for ( size_t i=0; i<nd; ++i )
    c += a(i)*b(i);

  // compute using cppmat

  V A(nd);
  V B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i) = a(i);

  for ( size_t i=0; i<nd; ++i )
    B(i) = b(i);

  double C = A.dot( B );
  double D = cm::dot( A , B );

  // check the result

  n = std::abs ( C - c );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = std::abs ( D - c );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("vector.dot( tensor2 )")
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd c(nd);

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(j) += a(i)*b(i,j);

  // compute using cppmat

  V A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i) = a(i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  V C = A.dot( B );
  V D = cm::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-C(i) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-D(i) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("vector.dot( tensor2s )")
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd c(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(j) += a(i)*b(i,j);

  // compute using cppmat

  V A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i) = a(i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  V C = A.dot( B );
  V D = cm::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-C(i) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-D(i) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("vector.dot( tensor2d )")
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::VectorXd c(nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(j) += a(i)*b(i,j);

  // compute using cppmat

  V A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i) = a(i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  V C = A.dot( B );
  V D = cm::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-C(i) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-D(i) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("vector.dyadic( vector )")
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(nd);
  Eigen::VectorXd b = Eigen::VectorXd::Random(nd);
  Eigen::MatrixXd c(nd,nd);

  c.setZero();

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) += a(i)*b(j);

  // compute using cppmat

  V A(nd);
  V B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i) = a(i);

  for ( size_t i=0; i<nd; ++i )
    B(i) = b(i);

  T2  C = A.dyadic( B );
  T2  D = cm::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j) - C(i,j) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j) - D(i,j) );

  REQUIRE( n < 1.e-12 );

  // check for symmetry

  REQUIRE ( ! C.issymmetric() );
  REQUIRE ( ! C.isdiagonal () );

}

// =================================================================================================

SECTION("vector.cross( vector )")
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(3);
  Eigen::VectorXd b = Eigen::VectorXd::Random(3);

  Eigen::Vector3d aa(a(0),a(1),a(2));
  Eigen::Vector3d bb(b(0),b(1),b(2));

  Eigen::Vector3d c = aa.cross(bb);

  // compute using cppmat

  V A(3);
  V B(3);

  for ( size_t i=0; i<3; ++i )
    A(i) = a(i);

  for ( size_t i=0; i<3; ++i )
    B(i) = b(i);

  cm::vector <double> C = A.cross( B );
  cm::vector <double> D = cm::cross( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<3; ++i )
    n += std::abs( c(i) - C(i) );

  REQUIRE( n < 1.e-12 );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<3; ++i )
    n += std::abs( c(i) - D(i) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("vector arithmetic")
{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(nd);
  Eigen::VectorXd b = Eigen::VectorXd::Random(nd);
  Eigen::VectorXd c = Eigen::VectorXd::Random(nd);
  Eigen::VectorXd d = Eigen::VectorXd::Random(nd);
  Eigen::VectorXd e = Eigen::VectorXd::Random(nd);
  Eigen::VectorXd f(nd);

  for ( size_t i=0; i<nd; ++i )
    f(i) = ( 7.*a(i) + 10./b(i) ) / ( c(i)+2. ) * ( d(i)-1. ) - 2.*e(i);

  // compute using cppmat

  cm::vector<double> A(nd);
  cm::vector<double> B(nd);
  cm::vector<double> C(nd);
  cm::vector<double> D(nd);
  cm::vector<double> E(nd);

  for ( size_t i=0; i<nd; ++i ) {
    A(i) = a(i);
    B(i) = b(i);
    C(i) = c(i);
    D(i) = d(i);
    E(i) = e(i);
  }

  cm::vector<double> F = ( 7.*A + 10./B ) / ( C+2. ) * ( D-1. ) - 2.*E;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( f(i)-F(i) );

  REQUIRE( n < 1.e-12 );

  // compute using cppmat

  F.zeros();

  A *= 7.;
  F += A;

  B  = 10./B;
  F += B;

  C += 2.;
  F /= C;

  D -= 1.;
  F *= D;

  E *= 2.;
  F -= E;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( f(i)-F(i) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2/tensor2d/tensor2d arithmetic")
{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd d = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd e = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd f(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = 0 ; j < nd ; ++j ) {
      if ( i != j ) {
        a(i,j) = 0.0;
        e(i,j) = 0.0;
      }
    }
  }

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      if ( i != j ) {
        c(j,i) = c(i,j);
        d(j,i) = d(i,j);
      }
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      f(i,j) = ( 7.*a(i,j) + 10./b(i,j) ) / ( c(i,j)+2. ) * ( d(i,j)-1. ) - 2.*e(i,j);

  // compute using cppmat

  T2d A(nd);
  T2  B(nd);
  T2s C(nd);
  T2s D(nd);
  T2d E(nd);

  for ( size_t i=0; i<nd; ++i ) {
    for ( size_t j=0; j<nd; ++j ) {
      A(i,j) = a(i,j);
      B(i,j) = b(i,j);
      C(i,j) = c(i,j);
      D(i,j) = d(i,j);
      E(i,j) = e(i,j);
    }
  }

  T2  F = ( 7.*A + 10./B ) / ( C+2. ) * ( D-1. ) - 2.*E;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( f(i,j)-F(i,j) );

  REQUIRE( n < 1.e-12 );

  // compute using cppmat

  F.zeros();

  A *= 7.;
  F += A;

  B  = 10./B;
  F += B;

  C += 2.;
  F /= C;

  D -= 1.;
  F *= D;

  E *= 2.;
  F -= E;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( f(i,j)-F(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2 * tensor2")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) * b(i,j);

  T2  A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2<double> C = A * B;

  A *= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2 / tensor2")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) / b(i,j);

  T2  A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2<double> C = A / B;

  A /= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2 + tensor2")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) + b(i,j);

  T2  A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2<double> C = A + B;

  A += B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2 - tensor2")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) - b(i,j);

  T2  A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2<double> C = A - B;

  A -= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2 * tensor2s")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) * b(i,j);

  T2  A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2<double> C = A * B;

  A *= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2 / tensor2s")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) / b(i,j);

  T2  A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2<double> C = A / B;

  A /= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2 + tensor2s")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) + b(i,j);

  T2  A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2<double> C = A + B;

  A += B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2 - tensor2s")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) - b(i,j);

  T2  A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2<double> C = A - B;

  A -= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2 * tensor2d")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) * b(i,j);

  T2  A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  T2d C = A * B;

  A *= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2 + tensor2d")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) + b(i,j);

  T2  A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  cm::tensor2<double> C = A + B;

  A += B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2 - tensor2d")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) - b(i,j);

  T2  A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  cm::tensor2<double> C = A - B;

  A -= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s * tensor2")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) * b(i,j);

  T2s A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2<double> C = A * B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s / tensor2")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) / b(i,j);

  T2s A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2<double> C = A / B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s + tensor2")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) + b(i,j);

  T2s A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2<double> C = A + B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s - tensor2")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) - b(i,j);

  T2s A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2<double> C = A - B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s * tensor2s")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) * b(i,j);

  T2s A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  T2s C = A * B;

  A *= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s / tensor2s")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) / b(i,j);

  T2s A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  T2s C = A / B;

  A /= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s + tensor2s")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) + b(i,j);

  T2s A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  T2s C = A + B;

  A += B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s - tensor2s")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) - b(i,j);

  T2s A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  T2s C = A - B;

  A -= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s * tensor2d")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) * b(i,j);

  T2s A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  T2d C = A * B;

  A *= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s + tensor2d")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) + b(i,j);

  T2s A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  T2s C = A + B;

  A += B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2s - tensor2d")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) - b(i,j);

  T2s A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  T2s C = A - B;

  A -= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d * tensor2")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) * b(i,j);

  T2d A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  T2d C = A * B;

  A *= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d / tensor2")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) / b(i,j);

  T2d A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  T2d C = A / B;

  A /= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d + tensor2")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) + b(i,j);

  T2d A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2<double> C = A + B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d - tensor2")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) - b(i,j);

  T2d A(nd);
  T2  B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2<double> C = A - B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d * tensor2s")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) * b(i,j);

  T2d A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2d<double> C = A * B;

  A *= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d / tensor2s")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) / b(i,j);

  T2d A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  T2d C = A / B;

  A /= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d + tensor2s")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) + b(i,j);

  T2d A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  T2s C = A + B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d - tensor2s")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i+1 ; j < nd ; ++j ) {
      b(j,i) = b(i,j);
    }
  }

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) - b(i,j);

  T2d A(nd);
  T2s B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cm::tensor2s<double> C = A - B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d + tensor2d", "cppmat::cartesian3d")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) + b(i,j);

  T2d A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  cm::tensor2d<double> C = A + B;

  A += B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d - tensor2d", "cppmat::cartesian3d")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) - b(i,j);

  T2d A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  cm::tensor2d<double> C = A - B;

  A -= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

// =================================================================================================

SECTION("tensor2d * tensor2d", "cppmat::cartesian3d")
{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        b(i,j) = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) * b(i,j);

  T2d A(nd);
  T2d B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  T2d C = A * B;

  A *= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  REQUIRE( n < 1.e-12 );

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  REQUIRE( n < 1.e-12 );
}

}
