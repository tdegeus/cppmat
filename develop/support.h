
#ifndef SUPPORT_H
#define SUPPORT_H

#include <Eigen/Eigen>

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

#endif
