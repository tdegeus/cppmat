
#include "support.h"

static const size_t ND = 11;

typedef cppmat::cartesian::tensor4 <double> T4;
typedef cppmat::cartesian::tensor2 <double> T2;
typedef cppmat::cartesian::tensor2s<double> T2s;
typedef cppmat::cartesian::tensor2d<double> T2d;
typedef cppmat::cartesian::vector  <double> V;

// =================================================================================================

TEST_CASE("cppmat::cartesian::vector", "var_cartesian_vector.h")
{

// =================================================================================================
// operations
// =================================================================================================

SECTION( "V.dot(V), length" )
{
  V A = V::Random(ND);

  double B = A.dot(A);

  double b = std::pow(A.length(), 2.);

  EQ(B, b);
}

// -------------------------------------------------------------------------------------------------

SECTION( "V.dot(T2)" )
{
  T2  I = T2 ::I(ND);
  V   A = V  ::Random(ND);

  V   B = A.dot(I);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "V.dot(T2s)" )
{
  T2s I = T2s::I(ND);
  V   A = V  ::Random(ND);

  V   B = A.dot(I);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "V.dot(T2d)" )
{
  T2d I = T2d::I(ND);
  V   A = V  ::Random(ND);

  V   B = A.dot(I);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "V.dyadic(V)" )
{
  V  A = V::Random(ND);
  V  B = V::Random(ND);
  V  C = V::Random(ND);

  V  D = C.dot(A.dyadic(B));
  V  E = C.dot(A) * B;

  Equal(D, E);
}

// -------------------------------------------------------------------------------------------------

SECTION( "V.cross(V)" )
{
  V  A = V::Random(3);
  V  B = V::Random(3);
  V  C = A.cross(B);

  double d = A.dot(C);
  double e = B.dot(C);

  EQ( d, 0.0 );
  EQ( e, 0.0 );
}

// =================================================================================================

}
