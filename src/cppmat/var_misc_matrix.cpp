/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_MISC_MATRIX_CPP
#define CPPMAT_VAR_MISC_MATRIX_CPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// =================================================================================================
// extra arithmetic operators : cppmat::matrix
// =================================================================================================

template<class X>
inline
cppmat::matrix<X>&
cppmat::matrix<X>::operator*= (
  const cppmat::symmetric::matrix<X> &B
)
{
  assert( this->shape() == B.shape() );

  size_t N = this->mShape[0];

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

template<class X>
inline
cppmat::matrix<X>&
cppmat::matrix<X>::operator/= (
  const cppmat::symmetric::matrix<X> &B
)
{
  assert( this->shape() == B.shape() );

  size_t N = this->mShape[0];

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

template<class X>
inline
cppmat::matrix<X>&
cppmat::matrix<X>::operator+= (
  const cppmat::symmetric::matrix<X> &B
)
{
  assert( this->shape() == B.shape() );

  size_t N = this->mShape[0];

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

template<class X>
inline
cppmat::matrix<X>&
cppmat::matrix<X>::operator-= (
  const cppmat::symmetric::matrix<X> &B
)
{
  assert( this->shape() == B.shape() );

  size_t N = this->mShape[0];

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

template<class X>
inline
cppmat::matrix<X>&
cppmat::matrix<X>::operator*= (
  const cppmat::diagonal::matrix<X> &B
)
{
  assert( this->shape() == B.shape() );

  size_t N = this->mShape[0];

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
      if ( i == j ) this->mData[i*N+i] *= B[i];
      else          this->mData[i*N+j]  = static_cast<X>(0);
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::matrix<X>&
cppmat::matrix<X>::operator+= (
  const cppmat::diagonal::matrix<X> &B
)
{
  assert( this->shape() == B.shape() );

  size_t N = this->mShape[0];

  for ( size_t i = 0 ; i < N ; ++i )
    this->mData[i*N+i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::matrix<X>&
cppmat::matrix<X>::operator-= (
  const cppmat::diagonal::matrix<X> &B
)
{
  assert( this->shape() == B.shape() );

  size_t N = this->mShape[0];

  for ( size_t i = 0 ; i < N ; ++i )
    this->mData[i*N+i] -= B[i];

  return *this;
}

// =================================================================================================
// extra arithmetic operators : cppmat::symmetric::matrix
// =================================================================================================

template<class X>
inline
cppmat::symmetric::matrix<X>&
cppmat::symmetric::matrix<X>::operator*= (
  const cppmat::diagonal::matrix<X> &B
)
{
  assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N; ++j ) {
      if ( i == j ) this->mData[i*N-(i-1)*i/2    ] *= B[i];
      else          this->mData[i*N-(i-1)*i/2+j-i]  = static_cast<X>(0);
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::symmetric::matrix<X>&
cppmat::symmetric::matrix<X>::operator+= (
  const cppmat::diagonal::matrix<X> &B
)
{
  assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N ; ++i )
    this->mData[i*N-(i-1)*i/2] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::symmetric::matrix<X>&
cppmat::symmetric::matrix<X>::operator-= (
  const cppmat::diagonal::matrix<X> &B
)
{
  assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N ; ++i )
    this->mData[i*N-(i-1)*i/2] -= B[i];

  return *this;
}

// =================================================================================================
// extra arithmetic operators : cppmat::diagonal::matrix
// =================================================================================================

template<class X>
inline
cppmat::diagonal::matrix<X>&
cppmat::diagonal::matrix<X>::operator*= (
  const cppmat::matrix<X> &B
)
{
  assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N; ++i )
    this->mData[i] *= B[i*N+i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::diagonal::matrix<X>&
cppmat::diagonal::matrix<X>::operator/= (
  const cppmat::matrix<X> &B
)
{
  assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N; ++i )
    this->mData[i] /= B[i*N+i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::diagonal::matrix<X>&
cppmat::diagonal::matrix<X>::operator*= (
  const cppmat::symmetric::matrix<X> &B
)
{
  assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N ; ++i )
    this->mData[i] *= B[i*N-(i-1)*i/2];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::diagonal::matrix<X>&
cppmat::diagonal::matrix<X>::operator/= (
  const cppmat::symmetric::matrix<X> &B
)
{
  assert( this->shape() == B.shape() );

  for ( size_t i = 0 ; i < N ; ++i )
    this->mData[i] /= B[i*N-(i-1)*i/2];

  return *this;
}

// =================================================================================================
// extra external arithmetic operators -> cppmat::matrix
// =================================================================================================

template<class X>
inline
cppmat::matrix<X> operator* (
  const cppmat::matrix<X> &A,
  const cppmat::symmetric::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::matrix<X> C(N,N);

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

template<class X>
inline
cppmat::matrix<X> operator/ (
  const cppmat::matrix<X> &A,
  const cppmat::symmetric::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::matrix<X> C(N,N);

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

template<class X>
inline
cppmat::matrix<X> operator+ (
  const cppmat::matrix<X> &A,
  const cppmat::symmetric::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::matrix<X> C(N,N);

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

template<class X>
inline
cppmat::matrix<X> operator- (
  const cppmat::matrix<X> &A,
  const cppmat::symmetric::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::matrix<X> C(N,N);

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

template<class X>
inline
cppmat::matrix<X> operator+ (
  const cppmat::matrix<X> &A,
  const cppmat::diagonal::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
      if ( i == j ) C[ i*N + j ] = A[ i*N + j ] + B[ i ];
      else          C[ i*N + j ] = A[ i*N + j ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::matrix<X> operator- (
  const cppmat::matrix<X> &A,
  const cppmat::diagonal::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
      if ( i == j ) C[ i*N + j ] = A[ i*N + j ] - B[ i ];
      else          C[ i*N + j ] = A[ i*N + j ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::matrix<X> operator* (
  const cppmat::symmetric::matrix<X> &A,
  const cppmat::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::matrix<X> C(N,N);

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

template<class X>
inline
cppmat::matrix<X> operator/ (
  const cppmat::symmetric::matrix<X> &A,
  const cppmat::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::matrix<X> C(N,N);

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

template<class X>
inline
cppmat::matrix<X> operator+ (
  const cppmat::symmetric::matrix<X> &A,
  const cppmat::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::matrix<X> C(N,N);

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

template<class X>
inline
cppmat::matrix<X> operator- (
  const cppmat::symmetric::matrix<X> &A,
  const cppmat::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::matrix<X> C(N,N);

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

template<class X>
inline
cppmat::matrix<X> operator+ (
  const cppmat::diagonal::matrix<X> &A,
  const cppmat::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
      if ( i == j ) C[ i*N + j ] = A[ i ] + B[ i*N + j ];
      else          C[ i*N + j ] =          B[ i*N + j ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::matrix<X> operator- (
  const cppmat::diagonal::matrix<X> &A,
  const cppmat::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
      if ( i == j ) C[ i*N + j ] = A[ i ] - B[ i*N + j ];
      else          C[ i*N + j ] =        - B[ i*N + j ];
    }
  }

  return C;
}

// =================================================================================================
// extra external arithmetic operators -> cppmat::symmetric::matrix
// =================================================================================================

template<class X>
inline
cppmat::symmetric::matrix<X> operator+ (
  const cppmat::symmetric::matrix<X> &A,
  const cppmat::diagonal::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::symmetric::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      if ( i == j ) C[ i*N - (i-1)*i/2         ] = A[ i*N - (i-1)*i/2         ] + B[ i ];
      else          C[ i*N - (i-1)*i/2 + j - i ] = A[ i*N - (i-1)*i/2 + j - i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::symmetric::matrix<X> operator- (
  const cppmat::symmetric::matrix<X> &A,
  const cppmat::diagonal::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::symmetric::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      if ( i == j ) C[ i*N - (i-1)*i/2         ] = A[ i*N - (i-1)*i/2         ] - B[ i ];
      else          C[ i*N - (i-1)*i/2 + j - i ] = A[ i*N - (i-1)*i/2 + j - i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::symmetric::matrix<X> operator+ (
  const cppmat::diagonal::matrix<X> &A,
  const X &B
)
{
  size_t N = A.shape(0);

  cppmat::symmetric::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      if ( i == j ) C[ i*N - (i-1)*i/2         ] = A[ i ] + B;
      else          C[ i*N - (i-1)*i/2 + j - i ] =          B;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::symmetric::matrix<X> operator- (
  const cppmat::diagonal::matrix<X> &A,
  const X &B
)
{
  size_t N = A.shape(0);

  cppmat::symmetric::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      if ( i == j ) C[ i*N - (i-1)*i/2         ] = A[ i ] - B;
      else          C[ i*N - (i-1)*i/2 + j - i ] =        - B;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::symmetric::matrix<X> operator+ (
  const cppmat::diagonal::matrix<X> &A,
  const cppmat::symmetric::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::symmetric::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      if ( i == j ) C[ i*N - (i-1)*i/2         ] = A[ i ] + B[ i*N - (i-1)*i/2         ];
      else          C[ i*N - (i-1)*i/2 + j - i ] =          B[ i*N - (i-1)*i/2 + j - i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::symmetric::matrix<X> operator- (
  const cppmat::diagonal::matrix<X> &A,
  const cppmat::symmetric::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::symmetric::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      if ( i == j ) C[ i*N - (i-1)*i/2         ] = A[ i ] - B[ i*N - (i-1)*i/2         ];
      else          C[ i*N - (i-1)*i/2 + j - i ] =        - B[ i*N - (i-1)*i/2 + j - i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::symmetric::matrix<X> operator+ (
  const X &A,
  const cppmat::diagonal::matrix<X> &B
)
{
  size_t N = B.shape(0);

  cppmat::symmetric::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      if ( i == j ) C[ i*N - (i-1)*i/2         ] = A + B[ i ];
      else          C[ i*N - (i-1)*i/2 + j - i ] = A;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::symmetric::matrix<X> operator- (
  const X &A,
  const cppmat::diagonal::matrix<X> &B
)
{
  size_t N = B.shape(0);

  cppmat::symmetric::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = i ; j < N ; ++j ) {
      if ( i == j ) C[ i*N - (i-1)*i/2         ] = A - B[ i ];
      else          C[ i*N - (i-1)*i/2 + j - i ] = A;
    }
  }

  return C;
}

// =================================================================================================
// extra external arithmetic operators -> cppmat::diagonal::matrix
// =================================================================================================

template<class X>
inline
cppmat::diagonal::matrix<X> operator* (
  const cppmat::diagonal::matrix<X> &A,
  const cppmat::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::diagonal::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i )
    C[i] = A[i] * B[ i*N + i ];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::diagonal::matrix<X> operator/ (
  const cppmat::diagonal::matrix<X> &A,
  const cppmat::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::diagonal::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i )
    C[i] = A[i] / B[ i*N + i ];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::diagonal::matrix<X> operator* (
  const cppmat::diagonal::matrix<X> &A,
  const cppmat::symmetric::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::diagonal::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i )
    C[i] = A[i] * B[ i*N - (i-1)*i/2 ];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::diagonal::matrix<X> operator/ (
  const cppmat::diagonal::matrix<X> &A,
  const cppmat::symmetric::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::diagonal::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i )
    C[i] = A[i] / B[ i*N - (i-1)*i/2 ];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::diagonal::matrix<X> operator* (
  const cppmat::matrix<X> &A,
  const cppmat::diagonal::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::diagonal::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i )
    C[i] = A[ i*N + i ] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::diagonal::matrix<X> operator* (
  const cppmat::symmetric::matrix<X> &A,
  const cppmat::diagonal::matrix<X> &B
)
{
  assert( A.shape() == B.shape() );

  size_t N = A.shape(0);

  cppmat::diagonal::matrix<X> C(N,N);

  for ( size_t i = 0 ; i < N ; ++i )
    C[i] = A[ i*N - (i-1)*i/2 ] * B[i];

  return C;
}

// =================================================================================================

#endif

