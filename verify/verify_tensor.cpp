/* =================================================================================================

Compile using:

$ clang++ `pkg-config --cflags Eigen3 cppmat` -std=c++14 -Wpedantic -Wall -o test verify_tensor.cpp
================================================================================================= */

#include <cppmat/tensor.h>
#include <Eigen/Eigen>

int main()
{

  double n;

  size_t nd = 6;


// =================================================================================================
// tensor4 arithmetic
// =================================================================================================

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

  cppmat::tensor4 <double> A(nd),B(nd),C(nd),D(nd),E(nd),F(nd);

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

  if ( n > 1.e-12 ) std::cout << "test   1   : failed   " << std::endl;
  else              std::cout << "test   1   : completed" << std::endl;

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

  if ( n > 1.e-12 ) std::cout << "test   2   : failed   " << std::endl;
  else              std::cout << "test   2   : completed" << std::endl;

}

// =================================================================================================
// tensor4.ddot( tensor4 )
// =================================================================================================

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

  cppmat::tensor4 <double> A(nd),B(nd);

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

  cppmat::tensor4 <double> C = A.ddot( B );
  cppmat::tensor4 <double> D = cppmat::ddot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-C(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test   3(a): failed   " << std::endl;
  else              std::cout << "test   3(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test   3(b): failed   " << std::endl;
  else              std::cout << "test   3(b): completed" << std::endl;
}

// =================================================================================================
// tensor4.ddot( tensor4.T() )
// =================================================================================================

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

  cppmat::tensor4 <double> A(nd),B(nd);

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

  cppmat::tensor4 <double> C = A.ddot( B.T() );
  cppmat::tensor4 <double> D = cppmat::ddot( A , cppmat::transpose(B) );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-C(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test   4(a): failed   " << std::endl;
  else              std::cout << "test   4(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test   4(b): failed   " << std::endl;
  else              std::cout << "test   4(b): completed" << std::endl;
}

// =================================================================================================
// tensor4.ddot( tensor4.LT() )
// =================================================================================================

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

  cppmat::tensor4 <double> A(nd),B(nd);

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

  cppmat::tensor4 <double> C = A.ddot( B.LT() );
  cppmat::tensor4 <double> D = cppmat::ddot( A , cppmat::transposeL(B) );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-C(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test   5(a): failed   " << std::endl;
  else              std::cout << "test   5(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test   5(b): failed   " << std::endl;
  else              std::cout << "test   5(b): completed" << std::endl;
}

// =================================================================================================
// tensor4.ddot( tensor4.RT() )
// =================================================================================================

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

  cppmat::tensor4 <double> A(nd),B(nd);

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

  cppmat::tensor4 <double> C = A.ddot( B.RT() );
  cppmat::tensor4 <double> D = cppmat::ddot( A , cppmat::transposeR(B) );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-C(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test   6(a): failed   " << std::endl;
  else              std::cout << "test   6(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test   6(b): failed   " << std::endl;
  else              std::cout << "test   6(b): completed" << std::endl;
}

// =================================================================================================
// tensor4.ddot( tensor2 )
// =================================================================================================

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

  cppmat::tensor4 <double> A(nd);
  cppmat::tensor2<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          A(i,j,k,l) = a(i*nd*nd*nd+j*nd*nd+k*nd+l);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A.ddot( B );
  cppmat::tensor2<double> D = cppmat::ddot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   7(a): failed   " << std::endl;
  else              std::cout << "test   7(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   7(b): failed   " << std::endl;
  else              std::cout << "test   7(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor4.ddot( tensor2s )
// =================================================================================================

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

  cppmat::tensor4 <double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          A(i,j,k,l) = a(i*nd*nd*nd+j*nd*nd+k*nd+l);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A.ddot( B );
  cppmat::tensor2<double> D = cppmat::ddot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   8(a): failed   " << std::endl;
  else              std::cout << "test   8(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   8(b): failed   " << std::endl;
  else              std::cout << "test   8(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor4.ddot( tensor2d )
// =================================================================================================

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

  cppmat::tensor4 <double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          A(i,j,k,l) = a(i*nd*nd*nd+j*nd*nd+k*nd+l);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A.ddot( B );
  cppmat::tensor2<double> D = cppmat::ddot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   9(a): failed   " << std::endl;
  else              std::cout << "test   9(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   9(b): failed   " << std::endl;
  else              std::cout << "test   9(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2.ddot( tensor4 )
// =================================================================================================

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

  cppmat::tensor2<double> A(nd);
  cppmat::tensor4 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          B(i,j,k,l) = b(i*nd*nd*nd+j*nd*nd+k*nd+l);

  cppmat::tensor2<double> C = A.ddot( B );
  cppmat::tensor2<double> D = cppmat::ddot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  10(a): failed   " << std::endl;
  else              std::cout << "test  10(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  10(b): failed   " << std::endl;
  else              std::cout << "test  10(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2.ddot( tensor2 )
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  double c = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c += a(i,j)*b(j,i);

  // compute using cppmat

  cppmat::tensor2<double> A(nd), B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  double C = A.ddot( B );
  double D = cppmat::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-12 ) std::cout << "test  11(a): failed   " << std::endl;
  else              std::cout << "test  11(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-12 ) std::cout << "test  11(b): failed   " << std::endl;
  else              std::cout << "test  11(b): completed" << std::endl;

}

// =================================================================================================
// tensor2.ddot( tensor2s )
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  double C = A.ddot( B );
  double D = cppmat::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-12 ) std::cout << "test  12(a): failed   " << std::endl;
  else              std::cout << "test  12(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-12 ) std::cout << "test  12(b): failed   " << std::endl;
  else              std::cout << "test  12(b): completed" << std::endl;

}

// =================================================================================================
// tensor2.ddot( tensor2d )
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  double C = A.ddot( B );
  double D = cppmat::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-12 ) std::cout << "test  13(a): failed   " << std::endl;
  else              std::cout << "test  13(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-12 ) std::cout << "test  13(b): failed   " << std::endl;
  else              std::cout << "test  13(b): completed" << std::endl;

}

// =================================================================================================
// tensor2.dot( tensor2 )
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2 <double> C = A.dot( B );
  cppmat::tensor2 <double> D = cppmat::dot ( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  14(a): failed   " << std::endl;
  else              std::cout << "test  14(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  14(b): failed   " << std::endl;
  else              std::cout << "test  14(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2.dot( tensor2s )
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2 <double> C = A.dot( B );
  cppmat::tensor2 <double> D = cppmat::dot ( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  15(a): failed   " << std::endl;
  else              std::cout << "test  15(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  15(b): failed   " << std::endl;
  else              std::cout << "test  15(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2.dot( tensor2d )
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2 <double> C = A.dot( B );
  cppmat::tensor2 <double> D = cppmat::dot ( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  16(a): failed   " << std::endl;
  else              std::cout << "test  16(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  16(b): failed   " << std::endl;
  else              std::cout << "test  16(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2.dot( vector )
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::vector  <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i) = b(i);

  cppmat::vector <double> C = A.dot( B );
  cppmat::vector <double> D = cppmat::dot( A, B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-C(i) );

  if ( n > 1.e-12 ) std::cout << "test  17(a): failed   " << std::endl;
  else              std::cout << "test  17(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-D(i) );

  if ( n > 1.e-12 ) std::cout << "test  17(b): failed   " << std::endl;
  else              std::cout << "test  17(b): completed" << std::endl;

}

// =================================================================================================
// tensor2.dyadic( tensor2 )
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor4 <double> C = A.dyadic( B );
  cppmat::tensor4 <double> D = cppmat::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  18(a): failed   " << std::endl;
  else              std::cout << "test  18(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  18(b): failed   " << std::endl;
  else              std::cout << "test  18(b): completed" << std::endl;

}

// =================================================================================================
// tensor2.dyadic( tensor2s )
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor4 <double> C = A.dyadic( B );
  cppmat::tensor4 <double> D = cppmat::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  19(a): failed   " << std::endl;
  else              std::cout << "test  19(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  19(b): failed   " << std::endl;
  else              std::cout << "test  19(b): completed" << std::endl;

}

// =================================================================================================
// tensor2.dyadic( tensor2d )
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  cppmat::tensor4 <double> C = A.dyadic( B );
  cppmat::tensor4 <double> D = cppmat::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  20(a): failed   " << std::endl;
  else              std::cout << "test  20(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  20(b): failed   " << std::endl;
  else              std::cout << "test  20(b): completed" << std::endl;

}

// =================================================================================================
// tensor2 - arithmetic
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2 <double> B(nd);
  cppmat::tensor2 <double> C(nd);
  cppmat::tensor2 <double> D(nd);
  cppmat::tensor2 <double> E(nd);

  for ( size_t i=0; i<nd; ++i ) {
    for ( size_t j=0; j<nd; ++j ) {
      A(i,j) = a(i,j);
      B(i,j) = b(i,j);
      C(i,j) = c(i,j);
      D(i,j) = d(i,j);
      E(i,j) = e(i,j);
    }
  }

  cppmat::tensor2 <double> F = ( 7.*A + 10./B ) / ( C+2. ) * ( D-1. ) - 2.*E;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( f(i,j)-F(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  21   : failed   " << std::endl;
  else              std::cout << "test  21   : completed" << std::endl;

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

  if ( n > 1.e-12 ) std::cout << "test  22   : failed   " << std::endl;
  else              std::cout << "test  22   : completed" << std::endl;

}

// =================================================================================================
// tensor2.dot( tensor2.T() )
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2 <double> C = A.dot( B.T() );
  cppmat::tensor2 <double> D = cppmat::dot ( A , cppmat::transpose( B ) );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  23(a): failed   " << std::endl;
  else              std::cout << "test  23(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  23(b): failed   " << std::endl;
  else              std::cout << "test  23(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2.trace()
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  double c = a.trace();

  // compute using cppmat

  cppmat::tensor2 <double> A(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  double C = A.trace();
  double D = cppmat::trace ( A );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-12 ) std::cout << "test  24(a): failed   " << std::endl;
  else              std::cout << "test  24(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-12 ) std::cout << "test  24(b): failed   " << std::endl;
  else              std::cout << "test  24(b): completed" << std::endl;

}

// =================================================================================================
// tensor2.det() -- 3D
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3,3);
  double c = a.determinant();

  // compute using cppmat

  cppmat::tensor2 <double> A(3);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      A(i,j) = a(i,j);

  double C = A.det();
  double D = cppmat::det ( A );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-12 ) std::cout << "test  25(a): failed   " << std::endl;
  else              std::cout << "test  25(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-12 ) std::cout << "test  25(b): failed   " << std::endl;
  else              std::cout << "test  25(b): completed" << std::endl;

}

// =================================================================================================
// tensor2.det() -- 2D
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(2,2);
  double c = a.determinant();

  // compute using cppmat

  cppmat::tensor2 <double> A(2);

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      A(i,j) = a(i,j);

  double C = A.det();
  double D = cppmat::det ( A );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-12 ) std::cout << "test  26(a): failed   " << std::endl;
  else              std::cout << "test  26(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-12 ) std::cout << "test  26(b): failed   " << std::endl;
  else              std::cout << "test  26(b): completed" << std::endl;

}

// =================================================================================================
// tensor2.inv() -- 3D
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3,3);
  Eigen::MatrixXd c = a.inverse();

  // compute using cppmat

  cppmat::tensor2 <double> A(3);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      A(i,j) = a(i,j);

  cppmat::tensor2 <double> C = A.inv();
  cppmat::tensor2 <double> D = cppmat::inv( A );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  27(a): failed   " << std::endl;
  else              std::cout << "test  27(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  27(b): failed   " << std::endl;
  else              std::cout << "test  27(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2.inv() -- 2D
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(2,2);
  Eigen::MatrixXd c = a.inverse();

  // compute using cppmat

  cppmat::tensor2 <double> A(2);

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      A(i,j) = a(i,j);

  cppmat::tensor2 <double> C = A.inv();
  cppmat::tensor2 <double> D = cppmat::inv( A );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  28(a): failed   " << std::endl;
  else              std::cout << "test  28(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  28(b): failed   " << std::endl;
  else              std::cout << "test  28(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2s.ddot( tensor4 )
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor4 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          B(i,j,k,l) = b(i*nd*nd*nd+j*nd*nd+k*nd+l);

  cppmat::tensor2 <double> C = A.ddot( B );
  cppmat::tensor2 <double> D = cppmat::ddot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  29(a): failed   " << std::endl;
  else              std::cout << "test  29(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  29(b): failed   " << std::endl;
  else              std::cout << "test  29(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2s.ddot( tensor2 )
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  double C = A.ddot( B );
  double D = cppmat::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-12 ) std::cout << "test  30(a): failed   " << std::endl;
  else              std::cout << "test  30(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-12 ) std::cout << "test  30(b): failed   " << std::endl;
  else              std::cout << "test  30(b): completed" << std::endl;

}

// =================================================================================================
// tensor2s.ddot( tensor2s )
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  double C = A.ddot( B );
  double D = cppmat::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-12 ) std::cout << "test  31(a): failed   " << std::endl;
  else              std::cout << "test  31(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-12 ) std::cout << "test  31(b): failed   " << std::endl;
  else              std::cout << "test  31(b): completed" << std::endl;

}

// =================================================================================================
// tensor2s.ddot( tensor2d )
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  double C = A.ddot( B );
  double D = cppmat::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-12 ) std::cout << "test  32(a): failed   " << std::endl;
  else              std::cout << "test  32(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-12 ) std::cout << "test  32(b): failed   " << std::endl;
  else              std::cout << "test  32(b): completed" << std::endl;

}

// =================================================================================================
// tensor2s.dot( tensor2 )
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2 <double> C = A.dot( B );
  cppmat::tensor2 <double> D = cppmat::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  33(a): failed   " << std::endl;
  else              std::cout << "test  33(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  33(b): failed   " << std::endl;
  else              std::cout << "test  33(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2s.dot( tensor2s )
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2 <double> C = A.dot( B );
  cppmat::tensor2 <double> D = cppmat::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  34(a): failed   " << std::endl;
  else              std::cout << "test  34(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  34(b): failed   " << std::endl;
  else              std::cout << "test  34(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2s.dot( tensor2d )
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2 <double> C = A.dot( B );
  cppmat::tensor2 <double> D = cppmat::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  35(a): failed   " << std::endl;
  else              std::cout << "test  35(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  35(b): failed   " << std::endl;
  else              std::cout << "test  35(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2s.dot( vector )
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::vector  <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i) = b(i);

  cppmat::vector <double> C = A.dot( B );
  cppmat::vector <double> D = cppmat::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-C(i) );

  if ( n > 1.e-12 ) std::cout << "test  36(a): failed   " << std::endl;
  else              std::cout << "test  36(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-D(i) );

  if ( n > 1.e-12 ) std::cout << "test  36(b): failed   " << std::endl;
  else              std::cout << "test  36(b): completed" << std::endl;

}

// =================================================================================================
// tensor2s.dyadic( tensor2 )
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor4 <double> C = A.dyadic( B );
  cppmat::tensor4 <double> D = cppmat::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  37(a): failed   " << std::endl;
  else              std::cout << "test  37(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  37(b): failed   " << std::endl;
  else              std::cout << "test  37(b): completed" << std::endl;

}

// =================================================================================================
// tensor2s.dyadic( tensor2s )
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor4 <double> C = A.dyadic( B );
  cppmat::tensor4 <double> D = cppmat::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  38(a): failed   " << std::endl;
  else              std::cout << "test  38(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  38(b): failed   " << std::endl;
  else              std::cout << "test  38(b): completed" << std::endl;

}

// =================================================================================================
// tensor2s.dyadic( tensor2d )
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  cppmat::tensor4 <double> C = A.dyadic( B );
  cppmat::tensor4 <double> D = cppmat::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  39(a): failed   " << std::endl;
  else              std::cout << "test  39(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  39(b): failed   " << std::endl;
  else              std::cout << "test  39(b): completed" << std::endl;

}

// =================================================================================================
// tensor2s - arithmetic
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2s<double> B(nd);
  cppmat::tensor2s<double> C(nd);
  cppmat::tensor2s<double> D(nd);
  cppmat::tensor2s<double> E(nd);

  for ( size_t i=0; i<nd; ++i ) {
    for ( size_t j=i; j<nd; ++j ) {
      A(i,j) = a(i,j);
      B(i,j) = b(i,j);
      C(i,j) = c(i,j);
      D(i,j) = d(i,j);
      E(i,j) = e(i,j);
    }
  }

  cppmat::tensor2s<double> F = ( 7.*A + 10./B ) / ( C+2. ) * ( D-1. ) - 2.*E;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( f(i,j)-F(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  40   : failed   " << std::endl;
  else              std::cout << "test  40   : completed" << std::endl;

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

  if ( n > 1.e-12 ) std::cout << "test  41   : failed   " << std::endl;
  else              std::cout << "test  41   : completed" << std::endl;

}

// =================================================================================================
// tensor2s.dot( tensor2s.T() )
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2 <double> C = A.dot( B.T() );
  cppmat::tensor2 <double> D = cppmat::dot( A , cppmat::transpose( B ) );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  42(b): failed   " << std::endl;
  else              std::cout << "test  42(b): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  42(b): failed   " << std::endl;
  else              std::cout << "test  42(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2s.trace()
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  double c = a.trace();

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = i+1 ; j < nd ; ++j )
      a(j,i) = a(i,j);

  // compute using cppmat

  cppmat::tensor2s<double> A(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  double C = A.trace();
  double D = cppmat::trace( A );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-12 ) std::cout << "test  43(a): failed   " << std::endl;
  else              std::cout << "test  43(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-12 ) std::cout << "test  43(b): failed   " << std::endl;
  else              std::cout << "test  43(b): completed" << std::endl;

}

// =================================================================================================
// tensor2s.det() -- 3D
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3,3);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = i+1 ; j < 3 ; ++j )
      a(j,i) = a(i,j);

  double c = a.determinant();

  // compute using cppmat

  cppmat::tensor2s<double> A(3);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=i; j<3; ++j )
      A(i,j) = a(i,j);

  double C = A.det();
  double D = cppmat::det( A );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-10 ) std::cout << "test  44(a): failed   " << std::endl;
  else              std::cout << "test  44(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-10 ) std::cout << "test  44(b): failed   " << std::endl;
  else              std::cout << "test  44(b): completed" << std::endl;

}

// =================================================================================================
// tensor2s.det() -- 2D
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(2,2);

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = i+1 ; j < 2 ; ++j )
      a(j,i) = a(i,j);

  double c = a.determinant();

  // compute using cppmat

  cppmat::tensor2 <double> A(2);

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      A(i,j) = a(i,j);

  double C = A.det();
  double D = cppmat::det( A );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-10 ) std::cout << "test  45(a): failed   " << std::endl;
  else              std::cout << "test  45(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-10 ) std::cout << "test  45(b): failed   " << std::endl;
  else              std::cout << "test  45(b): completed" << std::endl;

}

// =================================================================================================
// tensor2s.inv() -- 3D
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3,3);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = i+1 ; j < 3 ; ++j )
      a(j,i) = a(i,j);

  Eigen::MatrixXd c = a.inverse();

  // compute using cppmat

  cppmat::tensor2s<double> A(3);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=i; j<3; ++j )
      A(i,j) = a(i,j);

  cppmat::tensor2s<double> C = A.inv();
  cppmat::tensor2s<double> D = cppmat::inv( A );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-10 ) std::cout << "test  46(a): failed   " << std::endl;
  else              std::cout << "test  46(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  46(b): failed   " << std::endl;
  else              std::cout << "test  46(b): completed" << std::endl;

  // check for symmetry

  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2s.inv() -- 2D
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(2,2);

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = i+1 ; j < 2 ; ++j )
      a(j,i) = a(i,j);

  Eigen::MatrixXd c = a.inverse();

  // compute using cppmat

  cppmat::tensor2s<double> A(2);

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      A(i,j) = a(i,j);

  cppmat::tensor2s<double> C = A.inv();
  cppmat::tensor2s<double> D = cppmat::inv( A );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-10 ) std::cout << "test  47(a): failed   " << std::endl;
  else              std::cout << "test  47(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  47(b): failed   " << std::endl;
  else              std::cout << "test  47(b): completed" << std::endl;

  // check for symmetry

  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2d.ddot( tensor4 )
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor4 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          B(i,j,k,l) = b(i*nd*nd*nd+j*nd*nd+k*nd+l);

  cppmat::tensor2 <double> C = A.ddot( B );
  cppmat::tensor2 <double> D = cppmat::ddot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  48(b): failed   " << std::endl;
  else              std::cout << "test  48(b): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  48(b): failed   " << std::endl;
  else              std::cout << "test  48(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2d.ddot( tensor2 )
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  double C = A.ddot( B );
  double D = cppmat::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-12 ) std::cout << "test  49(a): failed   " << std::endl;
  else              std::cout << "test  49(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-12 ) std::cout << "test  49(b): failed   " << std::endl;
  else              std::cout << "test  49(b): completed" << std::endl;

}

// =================================================================================================
// tensor2d.ddot( tensor2s )
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  double C = A.ddot( B );
  double D = cppmat::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-12 ) std::cout << "test  50(a): failed   " << std::endl;
  else              std::cout << "test  50(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-12 ) std::cout << "test  50(b): failed   " << std::endl;
  else              std::cout << "test  50(b): completed" << std::endl;

}

// =================================================================================================
// tensor2d.ddot( tensor2d )
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  double C = A.ddot( B );
  double D = cppmat::ddot( A , B );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-12 ) std::cout << "test  51(a): failed   " << std::endl;
  else              std::cout << "test  51(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-12 ) std::cout << "test  51(b): failed   " << std::endl;
  else              std::cout << "test  51(b): completed" << std::endl;

}

// =================================================================================================
// tensor2d.dot( tensor2 )
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2 <double> C = A.dot( B );
  cppmat::tensor2 <double> D = cppmat::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  52(b): failed   " << std::endl;
  else              std::cout << "test  52(b): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  52(b): failed   " << std::endl;
  else              std::cout << "test  52(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2d.dot( tensor2s )
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2 <double> C = A.dot( B );
  cppmat::tensor2 <double> D = cppmat::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  53(b): failed   " << std::endl;
  else              std::cout << "test  53(b): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  53(b): failed   " << std::endl;
  else              std::cout << "test  53(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2d.dot( tensor2d )
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2d<double> C = A.dot( B );
  cppmat::tensor2d<double> D = cppmat::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  54(b): failed   " << std::endl;
  else              std::cout << "test  54(b): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  54(b): failed   " << std::endl;
  else              std::cout << "test  54(b): completed" << std::endl;

}

// =================================================================================================
// tensor2d.dot( vector )
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::vector  <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    B(i) = b(i);

  cppmat::vector <double> C = A.dot( B );
  cppmat::vector <double> D = cppmat::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-C(i) );

  if ( n > 1.e-12 ) std::cout << "test  55(a): failed   " << std::endl;
  else              std::cout << "test  55(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-D(i) );

  if ( n > 1.e-12 ) std::cout << "test  55(b): failed   " << std::endl;
  else              std::cout << "test  55(b): completed" << std::endl;

}

// =================================================================================================
// tensor2d.dyadic( tensor2 )
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor4 <double> C = A.dyadic( B );
  cppmat::tensor4 <double> D = cppmat::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  56(a): failed   " << std::endl;
  else              std::cout << "test  56(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  56(b): failed   " << std::endl;
  else              std::cout << "test  56(b): completed" << std::endl;

}

// =================================================================================================
// tensor2d.dyadic( tensor2s )
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor4 <double> C = A.dyadic( B );
  cppmat::tensor4 <double> D = cppmat::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  57(a): failed   " << std::endl;
  else              std::cout << "test  57(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  57(b): failed   " << std::endl;
  else              std::cout << "test  57(b): completed" << std::endl;

}

// =================================================================================================
// tensor2d.dyadic( tensor2d )
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  cppmat::tensor4 <double> C = A.dyadic( B );
  cppmat::tensor4 <double> D = cppmat::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l) - C(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  58(a): failed   " << std::endl;
  else              std::cout << "test  58(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          n += std::abs( c(i*nd*nd*nd+j*nd*nd+k*nd+l)-D(i,j,k,l) );

  if ( n > 1.e-12 ) std::cout << "test  58(b): failed   " << std::endl;
  else              std::cout << "test  58(b): completed" << std::endl;

}

// =================================================================================================
// tensor2d - arithmetic
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2d<double> B(nd);
  cppmat::tensor2d<double> C(nd);
  cppmat::tensor2d<double> D(nd);
  cppmat::tensor2d<double> E(nd);

  for ( size_t i=0; i<nd; ++i ) {
    A(i,i) = a(i,i);
    B(i,i) = b(i,i);
    C(i,i) = c(i,i);
    D(i,i) = d(i,i);
    E(i,i) = e(i,i);
  }

  cppmat::tensor2d<double> F = ( 7.*A + 10.*B ) * ( C )  * ( D ) - 2.*E;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( f(i,j)-F(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  59   : failed   " << std::endl;
  else              std::cout << "test  59   : completed" << std::endl;

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

  if ( n > 1.e-12 ) std::cout << "test  60   : failed   " << std::endl;
  else              std::cout << "test  60   : completed" << std::endl;

}

// =================================================================================================
// tensor2d.dot( tensor2s.T() )
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2 <double> C = A.dot( B.T() );
  cppmat::tensor2 <double> D = cppmat::dot( A , cppmat::transpose( B ) );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  61(b): failed   " << std::endl;
  else              std::cout << "test  61(b): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  61(b): failed   " << std::endl;
  else              std::cout << "test  61(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// tensor2d.trace()
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  double c = a.trace();

  for ( size_t i = 0 ; i < nd ; ++i )
    for ( size_t j = 0 ; j < nd ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  // compute using cppmat

  cppmat::tensor2d<double> A(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  double C = A.trace();
  double D = cppmat::trace( A );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-12 ) std::cout << "test  62(a): failed   " << std::endl;
  else              std::cout << "test  62(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-12 ) std::cout << "test  62(b): failed   " << std::endl;
  else              std::cout << "test  62(b): completed" << std::endl;

}

// =================================================================================================
// tensor2d.det() -- 3D
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3,3);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  double c = a.determinant();

  // compute using cppmat

  cppmat::tensor2d<double> A(3);

  for ( size_t i=0; i<3; ++i )
    A(i,i) = a(i,i);

  double C = A.det();
  double D = cppmat::det( A );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-12 ) std::cout << "test  63(a): failed   " << std::endl;
  else              std::cout << "test  63(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-12 ) std::cout << "test  63(b): failed   " << std::endl;
  else              std::cout << "test  63(b): completed" << std::endl;

}

// =================================================================================================
// tensor2d.det() -- 2D
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(2,2);

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  double c = a.determinant();

  // compute using cppmat

  cppmat::tensor2 <double> A(2);

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      A(i,j) = a(i,j);

  double C = A.det();
  double D = cppmat::det( A );

  // check the result

  n = std::abs( C - c );

  if ( n > 1.e-12 ) std::cout << "test  64(a): failed   " << std::endl;
  else              std::cout << "test  64(a): completed" << std::endl;

  // check the result

  n = std::abs( D - c );

  if ( n > 1.e-12 ) std::cout << "test  64(b): failed   " << std::endl;
  else              std::cout << "test  64(b): completed" << std::endl;

}

// =================================================================================================
// tensor2d.inv() -- 3D
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3,3);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  Eigen::MatrixXd c = a.inverse();

  // compute using cppmat

  cppmat::tensor2d<double> A(3);

  for ( size_t i=0; i<3; ++i )
    A(i,i) = a(i,i);

  cppmat::tensor2d<double> C = A.inv();
  cppmat::tensor2d<double> D = cppmat::inv( A );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-10 ) std::cout << "test  65(a): failed   " << std::endl;
  else              std::cout << "test  65(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  65(b): failed   " << std::endl;
  else              std::cout << "test  65(b): completed" << std::endl;

}

// =================================================================================================
// tensor2d.inv() -- 2D
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(2,2);

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      if ( i != j )
        a(i,j) = 0.0;

  Eigen::MatrixXd c = a.inverse();

  // compute using cppmat

  cppmat::tensor2d<double> A(2);

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      A(i,j) = a(i,j);

  cppmat::tensor2d<double> C = A.inv();
  cppmat::tensor2d<double> D = cppmat::inv( A );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-10 ) std::cout << "test  66(a): failed   " << std::endl;
  else              std::cout << "test  66(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      n += std::abs( c(i,j)-D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  66(b): failed   " << std::endl;
  else              std::cout << "test  66(b): completed" << std::endl;

}

// =================================================================================================
// vector.dot( vector )
// =================================================================================================

{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(nd);
  Eigen::VectorXd b = Eigen::VectorXd::Random(nd);
  double c = 0.;

  for ( size_t i=0; i<nd; ++i )
    c += a(i)*b(i);

  // compute using cppmat

  cppmat::vector  <double> A(nd);
  cppmat::vector  <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i) = a(i);

  for ( size_t i=0; i<nd; ++i )
    B(i) = b(i);

  double C = A.dot( B );
  double D = cppmat::dot( A , B );

  // check the result

  n = std::abs ( C - c );

  if ( n > 1.e-12 ) std::cout << "test  67(a): failed   " << std::endl;
  else              std::cout << "test  67(a): completed" << std::endl;

  // check the result

  n = std::abs ( D - c );

  if ( n > 1.e-12 ) std::cout << "test  68(a): failed   " << std::endl;
  else              std::cout << "test  68(a): completed" << std::endl;

}

// =================================================================================================
// vector.dot( tensor2 )
// =================================================================================================

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

  cppmat::vector  <double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i) = a(i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::vector  <double> C = A.dot( B );
  cppmat::vector  <double> D = cppmat::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-C(i) );

  if ( n > 1.e-12 ) std::cout << "test  68(a): failed   " << std::endl;
  else              std::cout << "test  68(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-D(i) );

  if ( n > 1.e-12 ) std::cout << "test  68(b): failed   " << std::endl;
  else              std::cout << "test  68(b): completed" << std::endl;

}

// =================================================================================================
// vector.dot( tensor2s )
// =================================================================================================

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

  cppmat::vector  <double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i) = a(i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::vector  <double> C = A.dot( B );
  cppmat::vector  <double> D = cppmat::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-C(i) );

  if ( n > 1.e-12 ) std::cout << "test  69(a): failed   " << std::endl;
  else              std::cout << "test  69(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-D(i) );

  if ( n > 1.e-12 ) std::cout << "test  69(b): failed   " << std::endl;
  else              std::cout << "test  69(b): completed" << std::endl;

}

// =================================================================================================
// vector.dot( tensor2d )
// =================================================================================================

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

  cppmat::vector  <double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i) = a(i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::vector  <double> C = A.dot( B );
  cppmat::vector  <double> D = cppmat::dot( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-C(i) );

  if ( n > 1.e-12 ) std::cout << "test  70(a): failed   " << std::endl;
  else              std::cout << "test  70(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( c(i)-D(i) );

  if ( n > 1.e-12 ) std::cout << "test  70(b): failed   " << std::endl;
  else              std::cout << "test  70(b): completed" << std::endl;

}

// =================================================================================================
// vector.dyadic( vector )
// =================================================================================================

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

  cppmat::vector  <double> A(nd);
  cppmat::vector  <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i) = a(i);

  for ( size_t i=0; i<nd; ++i )
    B(i) = b(i);

  cppmat::tensor2 <double> C = A.dyadic( B );
  cppmat::tensor2 <double> D = cppmat::dyadic( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j) - C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  71(a): failed   " << std::endl;
  else              std::cout << "test  71(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j) - D(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  71(b): failed   " << std::endl;
  else              std::cout << "test  71(b): completed" << std::endl;

  // check for symmetry

  if ( C.issymmetric() ) std::cout << "Result unexpected symmetric" << std::endl;
  if ( C.isdiagonal () ) std::cout << "Result unexpected diagonal"  << std::endl;

}

// =================================================================================================
// vector.cross( vector )
// =================================================================================================

{
  // compute using Eigen

  Eigen::VectorXd a = Eigen::VectorXd::Random(3);
  Eigen::VectorXd b = Eigen::VectorXd::Random(3);

  Eigen::Vector3d aa(a(0),a(1),a(2));
  Eigen::Vector3d bb(b(0),b(1),b(2));

  Eigen::Vector3d c = aa.cross(bb);

  // compute using cppmat

  cppmat::vector  <double> A(3);
  cppmat::vector  <double> B(3);

  for ( size_t i=0; i<3; ++i )
    A(i) = a(i);

  for ( size_t i=0; i<3; ++i )
    B(i) = b(i);

  cppmat::vector <double> C = A.cross( B );
  cppmat::vector <double> D = cppmat::cross( A , B );

  // check the result

  n = 0.0;

  for ( size_t i=0; i<3; ++i )
    n += std::abs( c(i) - C(i) );

  if ( n > 1.e-12 ) std::cout << "test  72(a): failed   " << std::endl;
  else              std::cout << "test  72(a): completed" << std::endl;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<3; ++i )
    n += std::abs( c(i) - D(i) );

  if ( n > 1.e-12 ) std::cout << "test  72(b): failed   " << std::endl;
  else              std::cout << "test  72(b): completed" << std::endl;

}

// =================================================================================================
// tensor2d - arithmetic
// =================================================================================================

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

  cppmat::vector<double> A(nd);
  cppmat::vector<double> B(nd);
  cppmat::vector<double> C(nd);
  cppmat::vector<double> D(nd);
  cppmat::vector<double> E(nd);

  for ( size_t i=0; i<nd; ++i ) {
    A(i) = a(i);
    B(i) = b(i);
    C(i) = c(i);
    D(i) = d(i);
    E(i) = e(i);
  }

  cppmat::vector<double> F = ( 7.*A + 10./B ) / ( C+2. ) * ( D-1. ) - 2.*E;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    n += std::abs( f(i)-F(i) );

  if ( n > 1.e-12 ) std::cout << "test  73   : failed   " << std::endl;
  else              std::cout << "test  73   : completed" << std::endl;

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

  if ( n > 1.e-12 ) std::cout << "test  74   : failed   " << std::endl;
  else              std::cout << "test  74   : completed" << std::endl;

}

// =================================================================================================
// tensor2 - arithmetic
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2 <double> B(nd);
  cppmat::tensor2s<double> C(nd);
  cppmat::tensor2s<double> D(nd);
  cppmat::tensor2d<double> E(nd);

  for ( size_t i=0; i<nd; ++i ) {
    for ( size_t j=0; j<nd; ++j ) {
      A(i,j) = a(i,j);
      B(i,j) = b(i,j);
      C(i,j) = c(i,j);
      D(i,j) = d(i,j);
      E(i,j) = e(i,j);
    }
  }

  cppmat::tensor2 <double> F = ( 7.*A + 10./B ) / ( C+2. ) * ( D-1. ) - 2.*E;

  // check the result

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( f(i,j)-F(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  75   : failed   " << std::endl;
  else              std::cout << "test  75   : completed" << std::endl;

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

  if ( n > 1.e-12 ) std::cout << "test  76   : failed   " << std::endl;
  else              std::cout << "test  76   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2 * tensor2
// =================================================================================================

{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) * b(i,j);

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A * B;

  A *= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  77   : failed   " << std::endl;
  else              std::cout << "test  77   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  78   : failed   " << std::endl;
  else              std::cout << "test  78   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2 / tensor2
// =================================================================================================

{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) / b(i,j);

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A / B;

  A /= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  79   : failed   " << std::endl;
  else              std::cout << "test  79   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  80   : failed   " << std::endl;
  else              std::cout << "test  80   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2 + tensor2
// =================================================================================================

{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) + b(i,j);

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A + B;

  A += B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  81   : failed   " << std::endl;
  else              std::cout << "test  81   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  82   : failed   " << std::endl;
  else              std::cout << "test  82   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2 - tensor2
// =================================================================================================

{
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(nd,nd);
  Eigen::MatrixXd c(nd,nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      c(i,j) = a(i,j) - b(i,j);

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A - B;

  A -= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  83   : failed   " << std::endl;
  else              std::cout << "test  83   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  84   : failed   " << std::endl;
  else              std::cout << "test  84   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2 * tensor2s
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A * B;

  A *= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  85   : failed   " << std::endl;
  else              std::cout << "test  85   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  86   : failed   " << std::endl;
  else              std::cout << "test  86   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2 / tensor2s
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A / B;

  A /= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  87   : failed   " << std::endl;
  else              std::cout << "test  87   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  88   : failed   " << std::endl;
  else              std::cout << "test  88   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2 + tensor2s
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A + B;

  A += B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  89   : failed   " << std::endl;
  else              std::cout << "test  89   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  90   : failed   " << std::endl;
  else              std::cout << "test  90   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2 - tensor2s
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A - B;

  A -= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  91   : failed   " << std::endl;
  else              std::cout << "test  91   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  92   : failed   " << std::endl;
  else              std::cout << "test  92   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2 * tensor2d
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  cppmat::tensor2d<double> C = A * B;

  A *= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  93   : failed   " << std::endl;
  else              std::cout << "test  93   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  94   : failed   " << std::endl;
  else              std::cout << "test  94   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2 + tensor2d
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  cppmat::tensor2<double> C = A + B;

  A += B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  95   : failed   " << std::endl;
  else              std::cout << "test  95   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  96   : failed   " << std::endl;
  else              std::cout << "test  96   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2 - tensor2d
// =================================================================================================

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

  cppmat::tensor2 <double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  cppmat::tensor2<double> C = A - B;

  A -= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  97   : failed   " << std::endl;
  else              std::cout << "test  97   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  98   : failed   " << std::endl;
  else              std::cout << "test  98   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2s * tensor2
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A * B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  99   : failed   " << std::endl;
  else              std::cout << "test  99   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2s / tensor2
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A / B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 100   : failed   " << std::endl;
  else              std::cout << "test 100   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2s + tensor2
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A + B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 101   : failed   " << std::endl;
  else              std::cout << "test 101   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2s - tensor2
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A - B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 102   : failed   " << std::endl;
  else              std::cout << "test 102   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2s * tensor2s
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2s<double> C = A * B;

  A *= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 103   : failed   " << std::endl;
  else              std::cout << "test 103   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 104   : failed   " << std::endl;
  else              std::cout << "test 104   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2s / tensor2s
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2s<double> C = A / B;

  A /= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 105   : failed   " << std::endl;
  else              std::cout << "test 105   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 106   : failed   " << std::endl;
  else              std::cout << "test 106   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2s + tensor2s
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2s<double> C = A + B;

  A += B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 107   : failed   " << std::endl;
  else              std::cout << "test 107   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 108   : failed   " << std::endl;
  else              std::cout << "test 108   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2s - tensor2s
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2s<double> C = A - B;

  A -= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 109   : failed   " << std::endl;
  else              std::cout << "test 109   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 110   : failed   " << std::endl;
  else              std::cout << "test 110   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2s * tensor2d
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  cppmat::tensor2d<double> C = A * B;

  A *= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 111   : failed   " << std::endl;
  else              std::cout << "test 111   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 112   : failed   " << std::endl;
  else              std::cout << "test 112   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2s + tensor2d
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  cppmat::tensor2s<double> C = A + B;

  A += B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 113   : failed   " << std::endl;
  else              std::cout << "test 113   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 114   : failed   " << std::endl;
  else              std::cout << "test 114   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2s - tensor2d
// =================================================================================================

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

  cppmat::tensor2s<double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      A(i,j) = a(i,j);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  cppmat::tensor2s<double> C = A - B;

  A -= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 116   : failed   " << std::endl;
  else              std::cout << "test 116   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 117   : failed   " << std::endl;
  else              std::cout << "test 117   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2d * tensor2
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2d<double> C = A * B;

  A *= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 118   : failed   " << std::endl;
  else              std::cout << "test 118   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 119   : failed   " << std::endl;
  else              std::cout << "test 119   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2d / tensor2
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2d<double> C = A / B;

  A /= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 120   : failed   " << std::endl;
  else              std::cout << "test 120   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 121   : failed   " << std::endl;
  else              std::cout << "test 121   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2d + tensor2
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A + B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 122   : failed   " << std::endl;
  else              std::cout << "test 122   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2d - tensor2
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2 <double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A - B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 123   : failed   " << std::endl;
  else              std::cout << "test 123   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2d * tensor2s
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A * B;

  A *= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 124   : failed   " << std::endl;
  else              std::cout << "test 124   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 125   : failed   " << std::endl;
  else              std::cout << "test 125   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2d / tensor2s
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2d<double> C = A / B;

  A /= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 126   : failed   " << std::endl;
  else              std::cout << "test 126   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 127   : failed   " << std::endl;
  else              std::cout << "test 127   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2d + tensor2s
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2s<double> C = A + B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 128   : failed   " << std::endl;
  else              std::cout << "test 128   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2d - tensor2s
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2s<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=i; j<nd; ++j )
      B(i,j) = b(i,j);

  cppmat::tensor2<double> C = A - B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 129   : failed   " << std::endl;
  else              std::cout << "test 129   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2d * tensor2d
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  cppmat::tensor2d<double> C = A * B;

  A *= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 130   : failed   " << std::endl;
  else              std::cout << "test 130   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 131   : failed   " << std::endl;
  else              std::cout << "test 131   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2d + tensor2d
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  cppmat::tensor2<double> C = A + B;

  A += B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 132   : failed   " << std::endl;
  else              std::cout << "test 132   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 133   : failed   " << std::endl;
  else              std::cout << "test 133   : completed" << std::endl;

}

// =================================================================================================
// arithmetic -- tensor2d - tensor2d
// =================================================================================================

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

  cppmat::tensor2d<double> A(nd);
  cppmat::tensor2d<double> B(nd);

  for ( size_t i=0; i<nd; ++i )
    A(i,i) = a(i,i);

  for ( size_t i=0; i<nd; ++i )
    B(i,i) = b(i,i);

  cppmat::tensor2<double> C = A - B;

  A -= B;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-C(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 134   : failed   " << std::endl;
  else              std::cout << "test 134   : completed" << std::endl;

  n = 0.0;

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      n += std::abs( c(i,j)-A(i,j) );

  if ( n > 1.e-12 ) std::cout << "test 135   : failed   " << std::endl;
  else              std::cout << "test 135   : completed" << std::endl;

}

  return 0;
}
