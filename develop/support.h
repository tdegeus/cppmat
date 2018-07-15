
#ifndef SUPPORT_H
#define SUPPORT_H

// -------------------------------------------------------------------------------------------------

#include <catch2/catch.hpp>

// support function to simplify the checks
#define EQ(a,b) REQUIRE_THAT( (a), Catch::WithinAbs((b), 1.e-10) );

// -------------------------------------------------------------------------------------------------

// #define CPPMAT_NOCONVERT
// #include <cppmat/cppmat.h>
#include "../src/cppmat/cppmat.h"

// -------------------------------------------------------------------------------------------------

#include <Eigen/Eigen>

// alias Eigen types
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatD;
typedef Eigen::Matrix<double, Eigen::Dynamic,              1, Eigen::ColMajor> ColD;

// =================================================================================================

inline MatD makeSymmetric(const MatD &A)
{
  return ( A + A.transpose() ) / 2.;
}

// =================================================================================================

inline MatD makeDiagonal(const MatD &A)
{
  MatD B = A;

  for ( auto i = 0 ; i < B.rows() ; ++i )
    for ( auto j = 0 ; j < B.cols() ; ++j )
      if ( i != j )
        B(i,j) = 0.;

  return B;
}

// =================================================================================================

inline void Equal(const cppmat::array<double> &A, const cppmat::array<double> &B)
{
  REQUIRE( A.size()  == B.size()  );
  REQUIRE( A.shape() == B.shape() );

  for ( size_t i = 0 ; i < A.size() ; ++i )
    EQ( A[i], B[i] );
}

// -------------------------------------------------------------------------------------------------

inline void Equal(const cppmat::array<double> &A, const cppmat::symmetric::matrix<double> &B)
{
  REQUIRE( A.shape() == B.shape() );

  for ( size_t i = 0 ; i < A.shape(0) ; ++i )
    for ( size_t j = 0 ; j < A.shape(1) ; ++j )
      EQ( A(i,j), B(i,j) );
}

// -------------------------------------------------------------------------------------------------

inline void Equal(const cppmat::symmetric::matrix<double> &B, const cppmat::array<double> &A)
{
  REQUIRE( A.shape() == B.shape() );

  for ( size_t i = 0 ; i < A.shape(0) ; ++i )
    for ( size_t j = 0 ; j < A.shape(1) ; ++j )
      EQ( A(i,j), B(i,j) );
}

// -------------------------------------------------------------------------------------------------

inline void Equal(const cppmat::array<double> &A, const cppmat::diagonal::matrix<double> &B)
{
  REQUIRE( A.shape() == B.shape() );

  for ( size_t i = 0 ; i < A.shape(0) ; ++i )
    for ( size_t j = 0 ; j < A.shape(1) ; ++j )
      EQ( A(i,j), B(i,j) );
}

// -------------------------------------------------------------------------------------------------

inline void Equal(const cppmat::diagonal::matrix<double> &B, const cppmat::array<double> &A)
{
  REQUIRE( A.shape() == B.shape() );

  for ( size_t i = 0 ; i < A.shape(0) ; ++i )
    for ( size_t j = 0 ; j < A.shape(1) ; ++j )
      EQ( A(i,j), B(i,j) );
}

// -------------------------------------------------------------------------------------------------

inline void Equal(const cppmat::diagonal::matrix<double> &B, const cppmat::diagonal::matrix<double> &A)
{
  REQUIRE( A.shape() == B.shape() );

  for ( size_t i = 0 ; i < A.shape(0) ; ++i )
    for ( size_t j = 0 ; j < A.shape(1) ; ++j )
      EQ( A(i,j), B(i,j) );
}

// -------------------------------------------------------------------------------------------------

inline void Equal(const cppmat::array<double> &A, const MatD &B)
{
  REQUIRE( A.size() == static_cast<size_t>(B.size()) );

  for ( auto i = 0 ; i < B.rows() ; ++i )
    for ( auto j = 0 ; j < B.cols() ; ++j )
      EQ( A(i,j), B(i,j) );
}

// -------------------------------------------------------------------------------------------------

inline void Equal(const cppmat::symmetric::matrix<double> &A, const MatD &B)
{
  for ( auto i = 0 ; i < B.rows() ; ++i )
    for ( auto j = 0 ; j < B.cols() ; ++j )
      EQ( A(i,j), B(i,j) );
}

// -------------------------------------------------------------------------------------------------

inline void Equal(const cppmat::diagonal::matrix<double> &A, const MatD &B)
{
  for ( auto i = 0 ; i < B.rows() ; ++i )
    for ( auto j = 0 ; j < B.cols() ; ++j )
      EQ( A(i,j), B(i,j) );
}

// =================================================================================================

#endif
