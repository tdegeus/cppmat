
#include "support.h"

static const size_t ND = 11;

typedef cppmat::tiny::cartesian::tensor4 <double,ND> T4;
typedef cppmat::tiny::cartesian::tensor2 <double,ND> T2;
typedef cppmat::tiny::cartesian::tensor2s<double,ND> T2s;
typedef cppmat::tiny::cartesian::tensor2d<double,ND> T2d;
typedef cppmat::tiny::cartesian::vector  <double,ND> V;

// =================================================================================================

TEST_CASE("cppmat::tiny::cartesian::tensor2s", "var_cartesian_tensor2s.h")
{

// =================================================================================================
// unit tensors
// =================================================================================================

SECTION( "I, T2s.dot(T2)" )
{
  T2s I = T2s::I();
  T2  A = T2 ::Random();

  T2  B = I.dot(A);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "I, T2s.dot(T2s)" )
{
  T2s I = T2s::I();
  T2s A = T2s::Random();

  T2  B = I.dot(A);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "I, T2s.dot(T2d)" )
{
  T2s I = T2s::I();
  T2d A = T2d::Random();

  T2  B = I.dot(A);

  Equal(A, B);
}

// -------------------------------------------------------------------------------------------------

SECTION( "I, T2s.dot(V)" )
{
  T2s I = T2s::I();
  V   A = V  ::Random();

  V   B = I.dot(A);

  Equal(A, B);
}

// =================================================================================================
// tensor products
// =================================================================================================

SECTION( "I, T2s.dot(T4)" )
{
  T2s A = T2s::Random();
  T4  I = T4 ::I();

  T2  B = A.ddot(I);

  Equal(A, B);
}

// =================================================================================================
// tensor products
// =================================================================================================

SECTION( "T2s.ddot(T2), T2s.dot(T2), T2s.trace()" )
{
  T2s A = T2s::Random();
  T2  B = T2 ::Random();

  double C = A.ddot(B);

  double c = A.dot(B).trace();

  EQ(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2s.ddot(T2s), T2s.dot(T2s), T2s.trace()" )
{
  T2s A = T2s::Random();
  T2s B = T2s::Random();

  double C = A.ddot(B);

  double c = A.dot(B).trace();

  EQ(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2s.ddot(T2d), T2s.dot(T2d), T2s.trace()" )
{
  T2s A = T2s::Random();
  T2d B = T2d::Random();

  double C = A.ddot(B);

  double c = A.dot(B).trace();

  EQ(C, c);
}

// =================================================================================================
// tensor products
// =================================================================================================

SECTION( "T2s.dyadic(T2)" )
{
  T2s A = T2s::I();
  T2  B = T2 ::I();

  T4  C = A.dyadic(B);

  T4  D = T4::II();

  Equal(C, D);
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2s.dyadic(T2s)" )
{
  T2s A = T2s::I();
  T2s B = T2s::I();

  T4  C = A.dyadic(B);

  T4  D = T4::II();

  Equal(C, D);
}

// -------------------------------------------------------------------------------------------------

SECTION( "T2s.dyadic(T2d)" )
{
  T2s A = T2s::I();
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
  T2s A = T2s::Random();

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

  cppmat::tiny::cartesian::tensor2s<double,3> A = cppmat::tiny::cartesian::tensor2s<double,3>::CopyDense(a.data(), a.data()+a.size());
  double C = A.det();

  EQ(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION("T2s.det() -- 2D")
{
  MatD   a = makeSymmetric(MatD::Random(2,2));
  double c = a.determinant();

  cppmat::tiny::cartesian::tensor2s<double,2> A = cppmat::tiny::cartesian::tensor2s<double,2>::CopyDense(a.data(), a.data()+a.size());
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

  cppmat::tiny::cartesian::tensor2s<double,3> A = cppmat::tiny::cartesian::tensor2s<double,3>::CopyDense(a.data(), a.data()+a.size());
  cppmat::tiny::cartesian::tensor2s<double,3> C = A.inv();

  Equal(C, c);
}

// -------------------------------------------------------------------------------------------------

SECTION("T2s.inv() -- 2D")
{
  MatD a = makeSymmetric(MatD::Random(2,2));
  MatD c = a.inverse();

  cppmat::tiny::cartesian::tensor2s<double,2> A = cppmat::tiny::cartesian::tensor2s<double,2>::CopyDense(a.data(), a.data()+a.size());
  cppmat::tiny::cartesian::tensor2s<double,2> C = A.inv();

  Equal(C, c);
}

// =================================================================================================

}
