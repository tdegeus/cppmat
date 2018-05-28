
#include "support.h"

static const size_t ND = 11;

typedef cppmat::cartesian::tensor4 <double> T4;
typedef cppmat::cartesian::tensor2 <double> T2;
typedef cppmat::cartesian::tensor2s<double> T2s;
typedef cppmat::cartesian::tensor2d<double> T2d;
typedef cppmat::cartesian::vector  <double> V;

// =================================================================================================

TEST_CASE("cppmat::cartesian::tensor2s", "var_cartesian_tensor2s.h")
{

// =================================================================================================
// unit tensors
// =================================================================================================

SECTION( "I, T2s.dot(T2)" )
{
  T2s I = T2s::I(ND);
  T2  A = T2 ::Random(ND);

  T2  B = I.dot(A);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "I, T2s.dot(T2s)" )
{
  T2s I = T2s::I(ND);
  T2s A = T2s::Random(ND);

  T2  B = I.dot(A);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "I, T2s.dot(T2d)" )
{
  T2s I = T2s::I(ND);
  T2d A = T2d::Random(ND);

  T2  B = I.dot(A);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "I, T2s.dot(V)" )
{
  T2s I = T2s::I(ND);
  V   A = V  ::Random(ND);

  V   B = I.dot(A);

  Equal(A, B);
}

// =================================================================================================
// tensor products
// =================================================================================================

SECTION( "I, T2s.dot(T4)" )
{
  T2s A = T2s::Random(ND);
  T4  I = T4 ::I(ND);

  T2  B = A.ddot(I);

  Equal(A, B);
}

// =================================================================================================
// tensor products
// =================================================================================================

SECTION( "T2s.ddot(T2), T2s.dot(T2), T2s.trace()" )
{
  T2s A = T2s::Random(ND);
  T2  B = T2 ::Random(ND);

  double C = A.ddot(B);

  double c = A.dot(B).trace();

  EQ(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2s.ddot(T2s), T2s.dot(T2s), T2s.trace()" )
{
  T2s A = T2s::Random(ND);
  T2s B = T2s::Random(ND);

  double C = A.ddot(B);

  double c = A.dot(B).trace();

  EQ(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2s.ddot(T2d), T2s.dot(T2d), T2s.trace()" )
{
  T2s A = T2s::Random(ND);
  T2d B = T2d::Random(ND);

  double C = A.ddot(B);

  double c = A.dot(B).trace();

  EQ(C, c);
}

// =================================================================================================
// tensor products
// =================================================================================================

SECTION( "T2s.dyadic(T2)" )
{
  T2s A = T2s::I(ND);
  T2  B = T2 ::I(ND);

  T4  C = A.dyadic(B);

  T4  D = T4::II(ND);

  Equal(C, D);
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2s.dyadic(T2s)" )
{
  T2s A = T2s::I(ND);
  T2s B = T2s::I(ND);

  T4  C = A.dyadic(B);

  T4  D = T4::II(ND);

  Equal(C, D);
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2s.dyadic(T2d)" )
{
  T2s A = T2s::I(ND);
  T2d B = T2d::I(ND);

  T4  C = A.dyadic(B);

  T4  D = T4::II(ND);

  Equal(C, D);
}

// =================================================================================================
// transpositions
// =================================================================================================

SECTION( "T" )
{
  T2s A = T2s::Random(ND);

  T2s B = A.T();

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    for ( size_t j = 0 ; j < B.ndim() ; ++j )
      EQ( B(i,j), A(j,i) );
}

// =================================================================================================
// determinant
// =================================================================================================

SECTION("T2s.det() -- 3D")
{
  MatD   a = makeSymmetric(MatD::Random(3,3));
  double c = a.determinant();

  T2s    A = T2s::CopyDense(3, a.data(), a.data()+a.size());
  double C = A.det();

  EQ(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION("T2s.det() -- 2D")
{
  MatD   a = makeSymmetric(MatD::Random(2,2));
  double c = a.determinant();

  T2s    A = T2s::CopyDense(2, a.data(), a.data()+a.size());
  double C = A.det();

  EQ(C, c);
}

// =================================================================================================
// inverse
// =================================================================================================

SECTION("T2s.inv() -- 3D")
{
  MatD a = makeSymmetric(MatD::Random(3,3));
  MatD c = a.inverse();

  T2s A = T2s::CopyDense(3, a.data(), a.data()+a.size());
  T2s C = A.inv();

  Equal(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION("T2s.inv() -- 2D")
{
  MatD a = makeSymmetric(MatD::Random(2,2));
  MatD c = a.inverse();

  T2s A = T2s::CopyDense(2, a.data(), a.data()+a.size());
  T2s C = A.inv();

  Equal(C, c);
}

// =================================================================================================

}
