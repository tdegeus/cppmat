/* =================================================================================================

Compile using:

$ clang++ `pkg-config --cflags Eigen3 cppmat` -std=c++14 -Wpedantic -Wall -o test verify_matrix.cpp
================================================================================================= */

#include <cppmat/tensor.h>
#include <cppmat/matrix.h>
#include <Eigen/Eigen>

int main()
{

  double n;

// =================================================================================================
// matrix + matrix (2D)
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd c(10,10);

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      c(i,j) = a(i,j) + b(i,j);

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

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   1(a): failed   " << std::endl;
  else              std::cout << "test   1(a): completed" << std::endl;

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( A(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   1(b): failed   " << std::endl;
  else              std::cout << "test   1(b): completed" << std::endl;
}

// =================================================================================================
// matrix - matrix (2D)
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd c(10,10);

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      c(i,j) = a(i,j) - b(i,j);

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

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   2(a): failed   " << std::endl;
  else              std::cout << "test   2(a): completed" << std::endl;

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( A(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   2(b): failed   " << std::endl;
  else              std::cout << "test   2(b): completed" << std::endl;
}

// =================================================================================================
// matrix * matrix (2D)
// =================================================================================================

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

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   3(a): failed   " << std::endl;
  else              std::cout << "test   3(a): completed" << std::endl;

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( A(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   3(b): failed   " << std::endl;
  else              std::cout << "test   3(b): completed" << std::endl;
}

// =================================================================================================
// matrix / scalar (2D)
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a  = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd bb = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd c(10,10);
  double b = bb(0,0);
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

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   4(a): failed   " << std::endl;
  else              std::cout << "test   4(a): completed" << std::endl;

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( A(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   4(b): failed   " << std::endl;
  else              std::cout << "test   4(b): completed" << std::endl;
}

// =================================================================================================
// matrix + scalar (2D)
// =================================================================================================

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

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   5(a): failed   " << std::endl;
  else              std::cout << "test   5(a): completed" << std::endl;

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( A(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   5(b): failed   " << std::endl;
  else              std::cout << "test   5(b): completed" << std::endl;
}

// =================================================================================================
// matrix - scalar (2D)
// =================================================================================================

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

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   6(a): failed   " << std::endl;
  else              std::cout << "test   6(a): completed" << std::endl;

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( A(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   6(b): failed   " << std::endl;
  else              std::cout << "test   6(b): completed" << std::endl;
}

// =================================================================================================
// matrix * scalar (2D)
// =================================================================================================

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

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   7(a): failed   " << std::endl;
  else              std::cout << "test   7(a): completed" << std::endl;

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( A(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   7(b): failed   " << std::endl;
  else              std::cout << "test   7(b): completed" << std::endl;
}

// =================================================================================================
// matrix / scalar (2D)
// =================================================================================================

{
  // compute using Eigen

  Eigen::MatrixXd a  = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd bb = Eigen::MatrixXd::Random(10,10);
  Eigen::MatrixXd c(10,10);
  double b = bb(0,0);
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

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   8(a): failed   " << std::endl;
  else              std::cout << "test   8(a): completed" << std::endl;

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( A(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   8(b): failed   " << std::endl;
  else              std::cout << "test   8(b): completed" << std::endl;
}

// =================================================================================================
// matrix + scalar (2D)
// =================================================================================================

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

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test   9   : failed   " << std::endl;
  else              std::cout << "test   9   : completed" << std::endl;

}

// =================================================================================================
// matrix - scalar (2D)
// =================================================================================================

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

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  10   : failed   " << std::endl;
  else              std::cout << "test  10   : completed" << std::endl;

}

// =================================================================================================
// matrix * scalar (2D)
// =================================================================================================

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

  // verify

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  11   : failed   " << std::endl;
  else              std::cout << "test  11   : completed" << std::endl;

}

// =================================================================================================
// matrix / scalar (2D)
// =================================================================================================

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

  n = 0.0;

  for ( size_t i = 0 ; i < 10 ; ++i )
    for ( size_t j = 0 ; j < 10 ; ++j )
      n += std::abs( C(i,j) - c(i,j) );

  if ( n > 1.e-12 ) std::cout << "test  12   : failed   " << std::endl;
  else              std::cout << "test  12   : completed" << std::endl;

}

  return 0;
}
