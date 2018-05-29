
#include "support.h"

static const size_t ND = 11;

typedef cppmat::tiny::cartesian::tensor4 <double,ND> T4;
typedef cppmat::tiny::cartesian::tensor2 <double,ND> T2;
typedef cppmat::tiny::cartesian::tensor2s<double,ND> T2s;
typedef cppmat::tiny::cartesian::tensor2d<double,ND> T2d;
typedef cppmat::tiny::cartesian::vector  <double,ND> V;

// =================================================================================================

TEST_CASE("cppmat::tiny::cartesian::tensor2", "var_cartesian_tensor2.h")
{

// =================================================================================================
// unit tensors
// =================================================================================================

SECTION( "I, T2.dot(T2)" )
{
  T2  I = T2::I();
  T2  A = T2::Random();

  T2  B = I.dot(A);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "I, T2.dot(T2s)" )
{
  T2  I = T2 ::I();
  T2s A = T2s::Random();

  T2  B = I.dot(A);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "I, T2.dot(T2d)" )
{
  T2  I = T2 ::I();
  T2d A = T2d::Random();

  T2  B = I.dot(A);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "I, T2.dot(V)" )
{
  T2  I = T2::I();
  V   A = V ::Random();

  V   B = I.dot(A);

  Equal(A, B);
}

// =================================================================================================
// tensor products
// =================================================================================================

SECTION( "I, T2.dot(T4)" )
{
  T2  A = T2::Random();
  T4  I = T4::I();

  T2  B = A.ddot(I);

  Equal(A, B);
}

// =================================================================================================
// tensor products
// =================================================================================================

SECTION( "T2.ddot(T2), T2.dot(T2), T2.trace()" )
{
  T2 A = T2::Random();
  T2 B = T2::Random();

  double C = A.ddot(B);

  double c = A.dot(B).trace();

  EQ(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2.ddot(T2s), T2.dot(T2s), T2.trace()" )
{
  T2  A = T2 ::Random();
  T2s B = T2s::Random();

  double C = A.ddot(B);

  double c = A.dot(B).trace();

  EQ(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2.ddot(T2d), T2.dot(T2d), T2.trace()" )
{
  T2  A = T2 ::Random();
  T2d B = T2d::Random();

  double C = A.ddot(B);

  double c = A.dot(B).trace();

  EQ(C, c);
}

// =================================================================================================
// tensor products
// =================================================================================================

SECTION( "T2.dyadic(T2)" )
{
  T2  A = T2::I();
  T2  B = T2::I();

  T4  C = A.dyadic(B);

  T4  D = T4::II();

  Equal(C, D);
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2.dyadic(T2s)" )
{
  T2  A = T2 ::I();
  T2s B = T2s::I();

  T4  C = A.dyadic(B);

  T4  D = T4::II();

  Equal(C, D);
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2.dyadic(T2d)" )
{
  T2  A = T2 ::I();
  T2d B = T2d::I();

  T4  C = A.dyadic(B);

  T4  D = T4::II();

  Equal(C, D);
}

// =================================================================================================
// transpositions
// =================================================================================================

SECTION( "T" )
{
  T2 A = T2::Random();

  T2 B = A.T();

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    for ( size_t j = 0 ; j < B.ndim() ; ++j )
      EQ( B(i,j), A(j,i) );
}

// =================================================================================================

}