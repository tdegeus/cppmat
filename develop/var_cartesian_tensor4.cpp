
#include "support.h"

static const size_t ND = 11;

typedef cppmat::cartesian::tensor4 <double> T4;
typedef cppmat::cartesian::tensor2 <double> T2;
typedef cppmat::cartesian::tensor2s<double> T2s;
typedef cppmat::cartesian::tensor2d<double> T2d;
typedef cppmat::cartesian::vector  <double> V;

// =================================================================================================

TEST_CASE("cppmat::cartesian::tensor4", "var_cartesian_tensor4.h")
{

// =================================================================================================
// transpositions
// =================================================================================================

SECTION( "T" )
{
  T4 A = T4::Random(ND);

  T4 B = A.T();

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    for ( size_t j = 0 ; j < B.ndim() ; ++j )
      for ( size_t k = 0 ; k < B.ndim() ; ++k )
        for ( size_t l = 0 ; l < B.ndim() ; ++l )
          EQ( B(i,j,k,l), A(l,k,j,i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "RT" )
{
  T4 A = T4::Random(ND);

  T4 B = A.RT();

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    for ( size_t j = 0 ; j < B.ndim() ; ++j )
      for ( size_t k = 0 ; k < B.ndim() ; ++k )
        for ( size_t l = 0 ; l < B.ndim() ; ++l )
          EQ( B(i,j,k,l), A(i,j,l,k) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "LT" )
{
  T4 A = T4::Random(ND);

  T4 B = A.LT();

  for ( size_t i = 0 ; i < B.ndim() ; ++i )
    for ( size_t j = 0 ; j < B.ndim() ; ++j )
      for ( size_t k = 0 ; k < B.ndim() ; ++k )
        for ( size_t l = 0 ; l < B.ndim() ; ++l )
          EQ( B(i,j,k,l), A(j,i,k,l) );
}

// =================================================================================================
// unit tensors
// =================================================================================================

SECTION( "I, T4.ddot(T4)" )
{
  T4 I = T4::I(ND);
  T4 A = T4::Random(ND);

  T4 C = I.ddot(A);

  Equal(C, A);
}

// -------------------------------------------------------------------------------------------------

SECTION( "I, T4.ddot(T2)" )
{
  T4 I = T4::I(ND);
  T2 A = T2::Random(ND);

  T2 C = I.ddot(A);

  Equal(C, A);
}

// -------------------------------------------------------------------------------------------------

SECTION( "Irt, T4.ddot(T2)" )
{
  T4 I = T4::Irt(ND);
  T2 A = T2::Random(ND);

  T2 C = I.ddot(A);

  Equal(C, A.T());
}

// -------------------------------------------------------------------------------------------------

SECTION( "Is, T4.ddot(T2)" )
{
  T4 I = T4::Is(ND);
  T2 A = T2::Random(ND);

  T2 C = I.ddot(A);

  Equal(C, (A+A.T())/2.);
}

// -------------------------------------------------------------------------------------------------

SECTION( "Id, T4.ddot(T2)" )
{
  T4 I = T4::Id(ND);
  T2 A = T2::Random(ND);

  T2 C = I.ddot(A);

  Equal(C, A-A.trace()/static_cast<double>(ND)*T2::I(ND));
}

// -------------------------------------------------------------------------------------------------

SECTION( "Isd, T4.ddot(T2)" )
{
  T4 I = T4::Isd(ND);
  T2 A = T2::Random(ND);

  T2 C = I.ddot(A);

  T2 B = (A+A.T())/2.;

  Equal(C, B-B.trace()/static_cast<double>(ND)*T2::I(ND));
}

// -------------------------------------------------------------------------------------------------

SECTION( "II, T4.ddot(T2)" )
{
  T4 I = T4::II(ND);
  T2 A = T2::Random(ND);

  T2 C = I.ddot(A);

  Equal(C, A.trace()*T2::I(ND));
}

// =================================================================================================
// unit tensors
// =================================================================================================

SECTION( "I, T4.ddot(T2s)" )
{
  T4  I = T4 ::I(ND);
  T2s A = T2s::Random(ND);

  T2  C = I.ddot(A);

  Equal(C, A);
}

// -------------------------------------------------------------------------------------------------

SECTION( "Irt, T4.ddot(T2s)" )
{
  T4  I = T4 ::Irt(ND);
  T2s A = T2s::Random(ND);

  T2  C = I.ddot(A);

  Equal(C, A.T());
}

// -------------------------------------------------------------------------------------------------

SECTION( "Is, T4.ddot(T2s)" )
{
  T4  I = T4 ::Is(ND);
  T2s A = T2s::Random(ND);

  T2  C = I.ddot(A);

  Equal(C, (A+A.T())/2.);
}

// -------------------------------------------------------------------------------------------------

SECTION( "Id, T4.ddot(T2s)" )
{
  T4  I = T4 ::Id(ND);
  T2s A = T2s::Random(ND);

  T2 C = I.ddot(A);

  Equal(C, A-A.trace()/static_cast<double>(ND)*T2s::I(ND));
}

// -------------------------------------------------------------------------------------------------

SECTION( "Isd, T4.ddot(T2s)" )
{
  T4  I = T4 ::Isd(ND);
  T2s A = T2s::Random(ND);

  T2  C = I.ddot(A);

  T2s B = (A+A.T())/2.;

  Equal(C, B-B.trace()/static_cast<double>(ND)*T2s::I(ND));
}

// -------------------------------------------------------------------------------------------------

SECTION( "II, T4.ddot(T2s)" )
{
  T4  I = T4 ::II(ND);
  T2s A = T2s::Random(ND);

  T2 C = I.ddot(A);

  Equal(C, A.trace()*T2::I(ND));
}

// =================================================================================================
// unit tensors
// =================================================================================================

SECTION( "I, T4.ddot(T2d)" )
{
  T4  I = T4 ::I(ND);
  T2d A = T2d::Random(ND);

  T2  C = I.ddot(A);

  Equal(C, A);
}

// -------------------------------------------------------------------------------------------------

SECTION( "Irt, T4.ddot(T2d)" )
{
  T4  I = T4 ::Irt(ND);
  T2d A = T2d::Random(ND);

  T2  C = I.ddot(A);

  Equal(C, A.T());
}

// -------------------------------------------------------------------------------------------------

SECTION( "Is, T4.ddot(T2d)" )
{
  T4  I = T4 ::Is(ND);
  T2d A = T2d::Random(ND);

  T2  C = I.ddot(A);

  Equal(C, (A+A.T())/2.);
}

// -------------------------------------------------------------------------------------------------

SECTION( "Id, T4.ddot(T2d)" )
{
  T4  I = T4 ::Id(ND);
  T2d A = T2d::Random(ND);

  T2 C = I.ddot(A);

  Equal(C, A-A.trace()/static_cast<double>(ND)*T2d::I(ND));
}

// -------------------------------------------------------------------------------------------------

SECTION( "Isd, T4.ddot(T2d)" )
{
  T4  I = T4 ::Isd(ND);
  T2d A = T2d::Random(ND);

  T2  C = I.ddot(A);

  T2d B = (A+A.T())/2.;

  Equal(C, B-B.trace()/static_cast<double>(ND)*T2d::I(ND));
}

// -------------------------------------------------------------------------------------------------

SECTION( "II, T4.ddot(T2d)" )
{
  T4  I = T4 ::II(ND);
  T2d A = T2d::Random(ND);

  T2 C = I.ddot(A);

  Equal(C, A.trace()*T2::I(ND));
}

// =================================================================================================

}
