/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_MISC_MATRIX_HPP
#define CPPMAT_FIX_MISC_MATRIX_HPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// =================================================================================================
// extra arithmetic operators : cppmat::tiny::matrix
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N>&
cppmat::tiny::matrix<X,M,N>::operator*= (
  const cppmat::tiny::symmetric::matrix<X,M,N> &B
)
{
  Assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      // - extract value
      X b = B[i*N-(i-1)*i/2+j-i];
      // - store symmetrically
                    this->mData[i*N+j] *= b;
      if ( i != j ) this->mData[j*N+i] *= b;
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N>&
cppmat::tiny::matrix<X,M,N>::operator/= (
  const cppmat::tiny::symmetric::matrix<X,M,N> &B
)
{
  Assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      // - extract value
      X b = B[i*N-(i-1)*i/2+j-i];
      // - store symmetrically
                    this->mData[i*N+j] /= b;
      if ( i != j ) this->mData[j*N+i] /= b;
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N>&
cppmat::tiny::matrix<X,M,N>::operator+= (
  const cppmat::tiny::symmetric::matrix<X,M,N> &B
)
{
  Assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      // - extract value
      X b = B[i*N-(i-1)*i/2+j-i];
      // - store symmetrically
                    this->mData[i*N+j] += b;
      if ( i != j ) this->mData[j*N+i] += b;
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N>&
cppmat::tiny::matrix<X,M,N>::operator-= (
  const cppmat::tiny::symmetric::matrix<X,M,N> &B
)
{
  Assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      // - extract value
      X b = B[i*N-(i-1)*i/2+j-i];
      // - store symmetrically
                    this->mData[i*N+j] -= b;
      if ( i != j ) this->mData[j*N+i] -= b;
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N>&
cppmat::tiny::matrix<X,M,N>::operator*= (
  const cppmat::tiny::diagonal::matrix<X,M,N> &B
)
{
  Assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
      if ( i == j ) this->mData[i*N+i] *= B[i];
      else          this->mData[i*N+j]  = static_cast<X>(0);
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N>&
cppmat::tiny::matrix<X,M,N>::operator+= (
  const cppmat::tiny::diagonal::matrix<X,M,N> &B
)
{
  Assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N ; ++i )
    this->mData[i*N+i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N>&
cppmat::tiny::matrix<X,M,N>::operator-= (
  const cppmat::tiny::diagonal::matrix<X,M,N> &B
)
{
  Assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N ; ++i )
    this->mData[i*N+i] -= B[i];

  return *this;
}

// =================================================================================================
// extra arithmetic operators : cppmat::tiny::symmetric::matrix
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::symmetric::matrix<X,M,N>&
cppmat::tiny::symmetric::matrix<X,M,N>::operator*= (
  const cppmat::tiny::diagonal::matrix<X,M,N> &B
)
{
  Assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N; ++j ) {
      if ( i == j ) this->mData[i*N-(i-1)*i/2    ] *= B[i];
      else          this->mData[i*N-(i-1)*i/2+j-i]  = static_cast<X>(0);
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::symmetric::matrix<X,M,N>&
cppmat::tiny::symmetric::matrix<X,M,N>::operator+= (
  const cppmat::tiny::diagonal::matrix<X,M,N> &B
)
{
  Assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N ; ++i )
    this->mData[i*N-(i-1)*i/2] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::symmetric::matrix<X,M,N>&
cppmat::tiny::symmetric::matrix<X,M,N>::operator-= (
  const cppmat::tiny::diagonal::matrix<X,M,N> &B
)
{
  Assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N ; ++i )
    this->mData[i*N-(i-1)*i/2] -= B[i];

  return *this;
}

// =================================================================================================
// extra arithmetic operators : cppmat::tiny::diagonal::matrix
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::diagonal::matrix<X,M,N>&
cppmat::tiny::diagonal::matrix<X,M,N>::operator*= (
  const cppmat::tiny::matrix<X,M,N> &B
)
{
  Assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N; ++i )
    this->mData[i] *= B[i*N+i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::diagonal::matrix<X,M,N>&
cppmat::tiny::diagonal::matrix<X,M,N>::operator/= (
  const cppmat::tiny::matrix<X,M,N> &B
)
{
  Assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N; ++i )
    this->mData[i] /= B[i*N+i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::diagonal::matrix<X,M,N>&
cppmat::tiny::diagonal::matrix<X,M,N>::operator*= (
  const cppmat::tiny::symmetric::matrix<X,M,N> &B
)
{
  Assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N ; ++i )
    this->mData[i] *= B[i*N-(i-1)*i/2];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::diagonal::matrix<X,M,N>&
cppmat::tiny::diagonal::matrix<X,M,N>::operator/= (
  const cppmat::tiny::symmetric::matrix<X,M,N> &B
)
{
  Assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N ; ++i )
    this->mData[i] /= B[i*N-(i-1)*i/2];

  return *this;
}

// =================================================================================================
// extra external arithmetic operators -> cppmat::tiny::matrix
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N> operator* (
  const cppmat::tiny::matrix<X,M,N> &A,
  const cppmat::tiny::symmetric::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      // - extract value
      X b = B[ i*N - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*N + j ] = A[ i*N + j ] * b;
      if ( i != j ) C[ j*N + i ] = A[ j*N + i ] * b;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N> operator/ (
  const cppmat::tiny::matrix<X,M,N> &A,
  const cppmat::tiny::symmetric::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      // - extract value
      X b = B[ i*N - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*N + j ] = A[ i*N + j ] / b;
      if ( i != j ) C[ j*N + i ] = A[ j*N + i ] / b;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N> operator+ (
  const cppmat::tiny::matrix<X,M,N> &A,
  const cppmat::tiny::symmetric::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      // - extract value
      X b = B[ i*N - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*N + j ] = A[ i*N + j ] + b;
      if ( i != j ) C[ j*N + i ] = A[ j*N + i ] + b;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N> operator- (
  const cppmat::tiny::matrix<X,M,N> &A,
  const cppmat::tiny::symmetric::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      // - extract value
      X b = B[ i*N - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*N + j ] = A[ i*N + j ] - b;
      if ( i != j ) C[ j*N + i ] = A[ j*N + i ] - b;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N> operator+ (
  const cppmat::tiny::matrix<X,M,N> &A,
  const cppmat::tiny::diagonal::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
      if ( i == j ) C[ i*N + j ] = A[ i*N + j ] + B[ i ];
      else          C[ i*N + j ] = A[ i*N + j ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N> operator- (
  const cppmat::tiny::matrix<X,M,N> &A,
  const cppmat::tiny::diagonal::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
      if ( i == j ) C[ i*N + j ] = A[ i*N + j ] - B[ i ];
      else          C[ i*N + j ] = A[ i*N + j ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N> operator* (
  const cppmat::tiny::symmetric::matrix<X,M,N> &A,
  const cppmat::tiny::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      // - extract value
      X a = A[ i*N - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*N + j ] = a * B[ i*N + j ];
      if ( i != j ) C[ j*N + i ] = a * B[ j*N + i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N> operator/ (
  const cppmat::tiny::symmetric::matrix<X,M,N> &A,
  const cppmat::tiny::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      // - extract value
      X a = A[ i*N - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*N + j ] = a / B[ i*N + j ];
      if ( i != j ) C[ j*N + i ] = a / B[ j*N + i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N> operator+ (
  const cppmat::tiny::symmetric::matrix<X,M,N> &A,
  const cppmat::tiny::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      // - extract value
      X a = A[ i*N - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*N + j ] = a + B[ i*N + j ];
      if ( i != j ) C[ j*N + i ] = a + B[ j*N + i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N> operator- (
  const cppmat::tiny::symmetric::matrix<X,M,N> &A,
  const cppmat::tiny::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      // - extract value
      X a = A[ i*N - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*N + j ] = a - B[ i*N + j ];
      if ( i != j ) C[ j*N + i ] = a - B[ j*N + i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N> operator+ (
  const cppmat::tiny::diagonal::matrix<X,M,N> &A,
  const cppmat::tiny::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
      if ( i == j ) C[ i*N + j ] = A[ i ] + B[ i*N + j ];
      else          C[ i*N + j ] =          B[ i*N + j ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::matrix<X,M,N> operator- (
  const cppmat::tiny::diagonal::matrix<X,M,N> &A,
  const cppmat::tiny::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
      if ( i == j ) C[ i*N + j ] = A[ i ] - B[ i*N + j ];
      else          C[ i*N + j ] =        - B[ i*N + j ];
    }
  }

  return C;
}

// =================================================================================================
// extra external arithmetic operators -> cppmat::tiny::symmetric::matrix
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::symmetric::matrix<X,M,N> operator+ (
  const cppmat::tiny::symmetric::matrix<X,M,N> &A,
  const cppmat::tiny::diagonal::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::symmetric::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      if ( i == j ) C[ i*N - (i-1)*i/2         ] = A[ i*N - (i-1)*i/2         ] + B[ i ];
      else          C[ i*N - (i-1)*i/2 + j - i ] = A[ i*N - (i-1)*i/2 + j - i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::symmetric::matrix<X,M,N> operator- (
  const cppmat::tiny::symmetric::matrix<X,M,N> &A,
  const cppmat::tiny::diagonal::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::symmetric::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      if ( i == j ) C[ i*N - (i-1)*i/2         ] = A[ i*N - (i-1)*i/2         ] - B[ i ];
      else          C[ i*N - (i-1)*i/2 + j - i ] = A[ i*N - (i-1)*i/2 + j - i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::symmetric::matrix<X,M,N> operator+ (
  const cppmat::tiny::diagonal::matrix<X,M,N> &A,
  const X &B
)
{
  cppmat::tiny::symmetric::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      if ( i == j ) C[ i*N - (i-1)*i/2         ] = A[ i ] + B;
      else          C[ i*N - (i-1)*i/2 + j - i ] =          B;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::symmetric::matrix<X,M,N> operator- (
  const cppmat::tiny::diagonal::matrix<X,M,N> &A,
  const X &B
)
{
  cppmat::tiny::symmetric::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      if ( i == j ) C[ i*N - (i-1)*i/2         ] = A[ i ] - B;
      else          C[ i*N - (i-1)*i/2 + j - i ] =        - B;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::symmetric::matrix<X,M,N> operator+ (
  const cppmat::tiny::diagonal::matrix<X,M,N> &A,
  const cppmat::tiny::symmetric::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::symmetric::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      if ( i == j ) C[ i*N - (i-1)*i/2         ] = A[ i ] + B[ i*N - (i-1)*i/2         ];
      else          C[ i*N - (i-1)*i/2 + j - i ] =          B[ i*N - (i-1)*i/2 + j - i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::symmetric::matrix<X,M,N> operator- (
  const cppmat::tiny::diagonal::matrix<X,M,N> &A,
  const cppmat::tiny::symmetric::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::symmetric::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      if ( i == j ) C[ i*N - (i-1)*i/2         ] = A[ i ] - B[ i*N - (i-1)*i/2         ];
      else          C[ i*N - (i-1)*i/2 + j - i ] =        - B[ i*N - (i-1)*i/2 + j - i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::symmetric::matrix<X,M,N> operator+ (
  const X &A,
  const cppmat::tiny::diagonal::matrix<X,M,N> &B
)
{
  cppmat::tiny::symmetric::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      if ( i == j ) C[ i*N - (i-1)*i/2         ] = A + B[ i ];
      else          C[ i*N - (i-1)*i/2 + j - i ] = A;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::symmetric::matrix<X,M,N> operator- (
  const X &A,
  const cppmat::tiny::diagonal::matrix<X,M,N> &B
)
{
  cppmat::tiny::symmetric::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      if ( i == j ) C[ i*N - (i-1)*i/2         ] = A - B[ i ];
      else          C[ i*N - (i-1)*i/2 + j - i ] = A;
    }
  }

  return C;
}

// =================================================================================================
// extra external arithmetic operators -> cppmat::tiny::diagonal::matrix
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::diagonal::matrix<X,M,N> operator* (
  const cppmat::tiny::diagonal::matrix<X,M,N> &A,
  const cppmat::tiny::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::diagonal::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i )
    C[i] = A[i] * B[ i*N + i ];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::diagonal::matrix<X,M,N> operator/ (
  const cppmat::tiny::diagonal::matrix<X,M,N> &A,
  const cppmat::tiny::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::diagonal::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i )
    C[i] = A[i] / B[ i*N + i ];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::diagonal::matrix<X,M,N> operator* (
  const cppmat::tiny::diagonal::matrix<X,M,N> &A,
  const cppmat::tiny::symmetric::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::diagonal::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i )
    C[i] = A[i] * B[ i*N - (i-1)*i/2 ];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::diagonal::matrix<X,M,N> operator/ (
  const cppmat::tiny::diagonal::matrix<X,M,N> &A,
  const cppmat::tiny::symmetric::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::diagonal::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i )
    C[i] = A[i] / B[ i*N - (i-1)*i/2 ];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::diagonal::matrix<X,M,N> operator* (
  const cppmat::tiny::matrix<X,M,N> &A,
  const cppmat::tiny::diagonal::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::diagonal::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i )
    C[i] = A[ i*N + i ] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
cppmat::tiny::diagonal::matrix<X,M,N> operator* (
  const cppmat::tiny::symmetric::matrix<X,M,N> &A,
  const cppmat::tiny::diagonal::matrix<X,M,N> &B
)
{
  Assert( A.shape() == B.shape() );

  cppmat::tiny::diagonal::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < N ; ++i )
    C[i] = A[ i*N - (i-1)*i/2 ] * B[i];

  return C;
}

// =================================================================================================

#endif

