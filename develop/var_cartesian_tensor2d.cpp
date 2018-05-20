
#include "support.h"

static const size_t ND = 11;

typedef cppmat::cartesian::tensor4 <double> T4;
typedef cppmat::cartesian::tensor2 <double> T2;
typedef cppmat::cartesian::tensor2s<double> T2s;
typedef cppmat::cartesian::tensor2d<double> T2d;
typedef cppmat::cartesian::vector  <double> V;

// =================================================================================================

TEST_CASE("cppmat::cartesian::tensor2d", "var_cartesian_tensor2s.h")
{

// =================================================================================================
// unit tensors
// =================================================================================================

SECTION( "I, T2d.dot(T2)" )
{
  T2d I = T2d::I(ND);
  T2  A = T2 ::Random(ND);

  T2  B = I.dot(A);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "I, T2d.dot(T2s)" )
{
  T2d I = T2d::I(ND);
  T2s A = T2s::Random(ND);

  T2  B = I.dot(A);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "I, T2d.dot(T2d)" )
{
  T2d I = T2d::I(ND);
  T2d A = T2d::Random(ND);

  T2d B = I.dot(A);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "I, T2d.dot(V)" )
{
  T2d I = T2d::I(ND);
  V   A = V  ::Random(ND);

  V   B = I.dot(A);

  Equal(A, B);
}

// =================================================================================================
// tensor products
// =================================================================================================

SECTION( "I, T2d.dot(T4)" )
{
  T2d A = T2d::Random(ND);
  T4  I = T4 ::I(ND);

  T2  B = A.ddot(I);

  Equal(A, B);
}

// =================================================================================================
// tensor products
// =================================================================================================

SECTION( "T2d.ddot(T2), T2d.dot(T2), T2d.trace()" )
{
  T2d A = T2d::Random(ND);
  T2  B = T2 ::Random(ND);

  double C = A.ddot(B);

  double c = A.dot(B).trace();

  EQ(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2d.ddot(T2s), T2d.dot(T2s), T2d.trace()" )
{
  T2d A = T2d::Random(ND);
  T2s B = T2s::Random(ND);

  double C = A.ddot(B);

  double c = A.dot(B).trace();

  EQ(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2d.ddot(T2d), T2d.dot(T2d), T2d.trace()" )
{
  T2d A = T2d::Random(ND);
  T2d B = T2d::Random(ND);

  double C = A.ddot(B);

  double c = A.dot(B).trace();

  EQ(C, c);
}

// =================================================================================================
// tensor products
// =================================================================================================

SECTION( "T2d.dyadic(T2)" )
{
  T2d A = T2d::I(ND);
  T2  B = T2 ::I(ND);

  T4  C = A.dyadic(B);

  T4  D = T4::II(ND);

  Equal(C, D);
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2d.dyadic(T2s)" )
{
  T2d A = T2d::I(ND);
  T2s B = T2s::I(ND);

  T4  C = A.dyadic(B);

  T4  D = T4::II(ND);

  Equal(C, D);
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2d.dyadic(T2d)" )
{
  T2d A = T2d::I(ND);
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
  T2d A = T2d::Random(ND);

  T2d B = A.T();

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    for ( size_t j = 0 ; j < B.ndim() ; ++j )
      EQ( B(i,j), A(j,i) );
}

// =================================================================================================
// determinant
// =================================================================================================

SECTION("T2d.det() -- 3D")
{
  MatD   a = makeDiagonal(MatD::Random(3,3));
  double c = a.determinant();

  T2s    A = T2s::CopyDense(3, a.data(), a.data()+a.size());
  double C = A.det();

  EQ(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION("T2d.det() -- 2D")
{
  MatD   a = makeDiagonal(MatD::Random(2,2));
  double c = a.determinant();

  T2s    A = T2s::CopyDense(2, a.data(), a.data()+a.size());
  double C = A.det();

  EQ(C, c);
}

// =================================================================================================
// inverse
// =================================================================================================

SECTION("T2d.inv() -- 3D")
{
  MatD a = makeDiagonal(MatD::Random(3,3));
  MatD c = a.inverse();

  T2s A = T2s::CopyDense(3, a.data(), a.data()+a.size());
  T2s C = A.inv();

  Equal(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION("T2d.inv() -- 2D")
{
  MatD a = makeDiagonal(MatD::Random(2,2));
  MatD c = a.inverse();

  T2s A = T2s::CopyDense(2, a.data(), a.data()+a.size());
  T2s C = A.inv();

  Equal(C, c);
}

// =================================================================================================

}
