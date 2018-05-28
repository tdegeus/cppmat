
#include "support.h"

static const size_t ND = 3;

typedef cppmat::tiny::cartesian::tensor4 <double,ND> T4;
typedef cppmat::tiny::cartesian::tensor2 <double,ND> T2;
typedef cppmat::tiny::cartesian::tensor2s<double,ND> T2s;
typedef cppmat::tiny::cartesian::tensor2d<double,ND> T2d;
typedef cppmat::tiny::cartesian::vector  <double,ND> V;

// =================================================================================================

TEST_CASE("cppmat::tiny::cartesian::vector<3>", "var_cartesian_vector.h")
{

// =================================================================================================
// operations
// =================================================================================================

SECTION( "V.dot(V), length" )
{
  V A = V::Random();

  double B = A.dot(A);

  double b = std::pow(A.length(), 2.);

  EQ(B, b);
}

// -------------------------------------------------------------------------------------------------

SECTION( "V.dot(T2)" )
{
  T2  I = T2 ::I();
  V   A = V  ::Random();

  V   B = A.dot(I);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "V.dot(T2s)" )
{
  T2s I = T2s::I();
  V   A = V  ::Random();

  V   B = A.dot(I);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "V.dot(T2d)" )
{
  T2d I = T2d::I();
  V   A = V  ::Random();

  V   B = A.dot(I);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "V.dyadic(V)" )
{
  V  A = V::Random();
  V  B = V::Random();
  V  C = V::Random();

  V  D = C.dot(A.dyadic(B));
  V  E = C.dot(A) * B;

  Equal(D, E);
}

// -------------------------------------------------------------------------------------------------

SECTION( "V.cross(V)" )
{
  V A = V::Random(3);
  V B = V::Random(3);
  V C = A.cross(B);

  double d = A.dot(C);
  double e = B.dot(C);

  EQ( d, 0.0 );
  EQ( e, 0.0 );
}

// =================================================================================================

}
