/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_SYMMETRIC_MATRIX_CPP
#define CPPMAT_VAR_SYMMETRIC_MATRIX_CPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace symmetric {

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline
matrix<X>::matrix(size_t m, size_t n)
{
  assert( m == n );

  resize(m,n);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X>::matrix(const matrix<X> &A)
{
  resize(A.shape(0), A.shape(1));

  setCopy(A.begin(), A.end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> matrix<X>::Arange(size_t m, size_t n)
{
  assert( m == n );

  matrix<X> out(m,n);

  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> matrix<X>::Zero(size_t m, size_t n)
{
  assert( m == n );

  matrix<X> out(m,n);

  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> matrix<X>::Ones(size_t m, size_t n)
{
  assert( m == n );

  matrix<X> out(m,n);

  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> matrix<X>::Constant(size_t m, size_t n, X D)
{
  assert( m == n );

  matrix<X> out(m,n);

  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
matrix<X> matrix<X>::Copy(size_t m, size_t n, Iterator first)
{
  assert( m == n );

  matrix<X> out(m,n);

  out.setCopy(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
matrix<X> matrix<X>::Copy(size_t m, size_t n, Iterator first, Iterator last)
{
  assert( m == n );

  matrix<X> out(m,n);

  out.setCopy(first,last);

  return out;
}

// =================================================================================================
// copy from dense storage
// =================================================================================================

template<class X>
template<typename Iterator>
inline
matrix<X> matrix<X>::CopyDense(size_t m, size_t n, Iterator first)
{
  assert( m == n );

  matrix<X> out(m,n);

  out.setCopyDense(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
matrix<X> matrix<X>::CopyDense(size_t m, size_t n, Iterator first, Iterator last)
{
  assert( m == n );

  matrix<X> out(m,n);

  out.setCopyDense(first,last);

  return out;
}

// =================================================================================================
// copy constructor
// =================================================================================================

#ifndef CPPMAT_NOCONVERT
template<class X>
inline
matrix<X>::operator cppmat::matrix<X> () const
{
  cppmat::matrix<X> out(N,N);

  for ( size_t i = 0 ; i < N ; ++i )
    for ( size_t j = i ; j < N ; ++j )
      out[i*N+j] = out[j*N+i] = mData[i*N-(i-1)*i/2+j-i];

  return out;
}
#endif

// =================================================================================================
// resize
// =================================================================================================

template<class X>
inline
void matrix<X>::resize(size_t m, size_t n)
{
  assert( m == n );

  // store old size
  size_t size = mSize;

  // copy to class member
  N = m;

  // set number of dimensions and total size
  mSize = (N+1)*N/2;

  // resize data container
  if ( mSize != size ) mData.resize(mSize);
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X>
inline
size_t matrix<X>::size() const
{
  return mSize;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t matrix<X>::rank() const
{
  return mRank;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t matrix<X>::shape(int i) const
{
  // check axis: (0,1,...,rank-1) or (-1,-2,...,-rank)
  assert( i  <      static_cast<int>(mRank) );
  assert( i >= -1 * static_cast<int>(mRank) );

  // return shape
  return N;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t matrix<X>::shape(size_t i) const
{
  // check axis: (0,1,...,rank-1)
  assert( i < mRank );

  // return shape
  return N;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
std::vector<size_t> matrix<X>::shape() const
{
  std::vector<size_t> shape(mRank);

  std::fill(shape.begin(), shape.end(), N);

  return shape;
}

// =================================================================================================
// index operators : operator[...]
// =================================================================================================

template<class X>
inline
X& matrix<X>::operator[](size_t i)
{
  assert( i < mSize );

  return mData[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
const X& matrix<X>::operator[](size_t i) const
{
  assert( i < mSize );

  return mData[i];
}

// =================================================================================================
// index operators : operator(...)
// =================================================================================================

template<class X>
inline
X& matrix<X>::operator()(size_t a, size_t b)
{
  assert( a < N );
  assert( b < N );

  if (a <= b) return mData[ a*N - (a-1)*a/2 + b - a ];
  else        return mData[ b*N - (b-1)*b/2 + a - b ];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
const X& matrix<X>::operator()(size_t a, size_t b) const
{
  assert( a < N );
  assert( b < N );

  if (a <= b) return mData[ a*N - (a-1)*a/2 + b - a ];
  else        return mData[ b*N - (b-1)*b/2 + a - b ];
}

// =================================================================================================
// index operators : compress(...)
// =================================================================================================

template<class X>
inline
size_t matrix<X>::compress(size_t a, size_t b) const
{
  assert( a < N );
  assert( b < N );

  if (a <= b) return a*N - (a-1)*a/2 + b - a;
  else        return b*N - (b-1)*b/2 + a - b;
}

// =================================================================================================
// index operators : decompress(...)
// =================================================================================================

template<class X>
inline
std::vector<size_t> matrix<X>::decompress(size_t i) const
{
  // check input
  assert( i < mSize );

  // allocate matrix-index
  std::vector<size_t> idx(mRank);

  // reconstruct
  idx[0] = 0;
  size_t keyafter;
  do
  {
    idx[0]++;
    keyafter = idx[0] * N - (idx[0] - 1) * idx[0] / 2;
  } while ( i >= keyafter );
  idx[0]--;
  idx[1] = N - keyafter + i;

  return idx;
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<class X>
inline
X* matrix<X>::data()
{
  return mData.data();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
const X* matrix<X>::data() const
{
  return mData.data();
}

// =================================================================================================
// iterators : begin() and end()
// =================================================================================================

template<class X>
inline
auto matrix<X>::begin()
{
  return mData.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto matrix<X>::begin() const
{
  return mData.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto matrix<X>::end()
{
  return mData.end();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto matrix<X>::end() const
{
  return mData.end();
}

// =================================================================================================
// iterators : index()
// =================================================================================================

template<class X>
inline
auto matrix<X>::index(size_t i)
{
  assert( i < mSize );

  return begin() + i;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto matrix<X>::index(size_t i) const
{
  assert( i < mSize );

  return begin() + i;
}

// =================================================================================================
// iterators : item()
// =================================================================================================

template<class X>
inline
auto matrix<X>::item(size_t a, size_t b)
{
  assert( a < N );
  assert( b < N );

  if (a <= b) return begin() + ( a*N - (a-1)*a/2 + b - a );
  else        return begin() + ( b*N - (b-1)*b/2 + a - b );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto matrix<X>::item(size_t a, size_t b) const
{
  assert( a < N );
  assert( b < N );

  if (a <= b) return begin() + ( a*N - (a-1)*a/2 + b - a );
  else        return begin() + ( b*N - (b-1)*b/2 + a - b );
}

// =================================================================================================
// initialize
// =================================================================================================

template<class X>
inline
void matrix<X>::setArange()
{
  for ( size_t i = 0 ; i < mSize ; ++i ) mData[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void matrix<X>::setZero()
{
  std::fill(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void matrix<X>::setOnes()
{
  std::fill(begin(), end(), static_cast<X>(1));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void matrix<X>::setConstant(X D)
{
  std::fill(begin(), end(), D);
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline
void matrix<X>::setCopy(Iterator first)
{
  std::copy(first, first+mSize, begin());
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline
void matrix<X>::setCopy(Iterator first, Iterator last)
{
  assert( mSize == last-first );

  std::copy(first, last, begin());
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
void matrix<X>::setCopyDense(Iterator first)
{
  // check for symmetry
  #ifndef NDEBUG
  for ( size_t i = 0 ; i < N ; ++i )
    for ( size_t j = i+1 ; j < N ; ++j )
      assert( first[i*N+j] == first[j*N+i] );
  #endif

  // copy from input (ignores lower diagonal terms)
  for ( size_t i = 0 ; i < N ; ++i )
    for ( size_t j = i ; j < N ; ++j )
      mData[i*N-(i-1)*i/2+j-i] = first[i*N+j];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
void matrix<X>::setCopyDense(Iterator first, Iterator last)
{
  // avoid compiler warning
  UNUSED(last);

  // check size
  assert( N*N == last-first );

  // check for symmetry
  #ifndef NDEBUG
  for ( size_t i = 0 ; i < N ; ++i )
    for ( size_t j = i+1 ; j < N ; ++j )
      assert( first[i*N+j] == first[j*N+i] );
  #endif

  // copy from input (ignores lower diagonal terms)
  for ( size_t i = 0 ; i < N ; ++i )
    for ( size_t j = i ; j < N ; ++j )
      mData[i*N-(i-1)*i/2+j-i] = first[i*N+j];
}

// =================================================================================================
// copy to target
// =================================================================================================

template<class X>
template<class Iterator>
inline
void matrix<X>::copyTo(Iterator first) const
{
  std::copy(begin(), end(), first);
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline
void matrix<X>::copyTo(Iterator first, Iterator last) const
{
  assert( mSize == last-first );

  UNUSED(last);

  std::copy(begin(), end(), first);
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline
void matrix<X>::copyToDense(Iterator first) const
{
  for ( size_t i = 0 ; i < N ; ++i )
    for ( size_t j = i ; j < N ; ++j )
      first[i*N+j] = first[j*N+i] = mData[i*N-(i-1)*i/2+j-i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline
void matrix<X>::copyToDense(Iterator first, Iterator last) const
{
  assert( N*N == last-first );

  UNUSED(last);

  for ( size_t i = 0 ; i < N ; ++i )
    for ( size_t j = i ; j < N ; ++j )
      first[i*N+j] = first[j*N+i] = mData[i*N-(i-1)*i/2+j-i];
}

// =================================================================================================
// sign change
// =================================================================================================

template<class X>
inline
matrix<X> matrix<X>::operator- () const
{
  matrix<X> out(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = -mData[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> matrix<X>::operator+ () const
{
  matrix<X> out(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = mData[i];

  return out;
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X>
inline
matrix<X>& matrix<X>::operator*= (const matrix<X> &B)
{
  assert( shape() == B.shape() );
  assert( rank () == B.rank () );
  assert( size () == B.size () );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X>& matrix<X>::operator/= (const matrix<X> &B)
{
  assert( shape() == B.shape() );
  assert( rank () == B.rank () );
  assert( size () == B.size () );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X>& matrix<X>::operator+= (const matrix<X> &B)
{
  assert( shape() == B.shape() );
  assert( rank () == B.rank () );
  assert( size () == B.size () );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X>& matrix<X>::operator-= (const matrix<X> &B)
{
  assert( shape() == B.shape() );
  assert( rank () == B.rank () );
  assert( size () == B.size () );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X>& matrix<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X>& matrix<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X>& matrix<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X>& matrix<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B;

  return *this;
}

// =================================================================================================
// absolute value
// =================================================================================================

template<class X>
inline
matrix<X> matrix<X>::abs() const
{
  matrix<X> out(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = std::abs(mData[i]);

  return out;
}

// =================================================================================================
// norm
// =================================================================================================

template<class X>
inline
X matrix<X>::norm() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += std::abs(mData[i]);

  return out;
}

// =================================================================================================
// location of the minimum/maximum
// =================================================================================================

template<class X>
inline
std::vector<size_t> matrix<X>::argmin() const
{
  return decompress( std::min_element(begin(), end()) - begin() );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
std::vector<size_t> matrix<X>::argmax() const
{
  return decompress( std::max_element(begin(), end()) - begin() );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t matrix<X>::argminIndex() const
{
  return std::min_element(begin(), end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t matrix<X>::argmaxIndex() const
{
  return std::max_element(begin(), end()) - begin();
}

// =================================================================================================
// minimum
// =================================================================================================

template<class X>
inline
X matrix<X>::minCoeff() const
{
  return *std::min_element(begin(),end());
}

// =================================================================================================
// maximum
// =================================================================================================

template<class X>
inline
X matrix<X>::maxCoeff() const
{
  return *std::max_element(begin(),end());
}

// =================================================================================================
// sum
// =================================================================================================

template<class X>
inline
X matrix<X>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t a = 0 ; a < N ; ++a ) {
    for ( size_t b = a ; b < N ; ++b ) {
      if (a == b) out += mData[ a*N - (a-1)*a/2 + b - a ];
      else        out += mData[ a*N - (a-1)*a/2 + b - a ] * static_cast<X>(2);
    }
  }

  return out;
}

// =================================================================================================
// mean
// =================================================================================================

template<class X>
inline
double matrix<X>::mean() const
{
  return static_cast<double>(sum())/static_cast<double>(N*N);
}

// =================================================================================================
// weighted average
// =================================================================================================

template<class X>
inline
double matrix<X>::average(const matrix<X> &weights, bool norm) const
{
  if ( norm ) return static_cast<double>((weights*(*this)).sum())/static_cast<double>(weights.sum());
  else        return static_cast<double>((weights*(*this)).sum());
}

// =================================================================================================
// find all non-zero entries
// =================================================================================================

template<class X>
inline
std::vector<size_t> matrix<X>::where() const
{
  size_t nnz = 0;

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] )
      ++nnz;

  std::vector<size_t> out(nnz);

  size_t j = 0;

  for ( size_t i = 0 ; i < mSize ; ++i ) {
    if ( mData[i] ) {
      out[j] = i;
      ++j;
    }
  }

  return out;
}

// =================================================================================================
// formatted print
// =================================================================================================

template<class X>
inline
std::ostream& operator<<(std::ostream& out, const matrix<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  for ( size_t i = 0 ; i < src.shape(0) ; ++i ) {
    for ( size_t j = 0 ; j < src.shape(1) ; ++j ) {
      out << std::setw(w) << std::setprecision(p) << src(i,j);
      if      ( j != src.shape(1)-1 ) out << ", ";
      else if ( i != src.shape(0)-1 ) out << ";" << std::endl;
      else                            out << ";";
    }
  }

  return out;
}

// =================================================================================================
// arithmetic operators: external
// =================================================================================================

template<class X>
inline
matrix<X> operator* (matrix<X> A, const matrix<X> &B)
{
  A *= B;

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> operator/ (matrix<X> A, const matrix<X> &B)
{
  A /= B;

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> operator+ (matrix<X> A, const matrix<X> &B)
{
  A += B;

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> operator- (matrix<X> A, const matrix<X> &B)
{
  A -= B;

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> operator* (matrix<X> A, const X &B)
{
  A *= B;

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> operator/ (matrix<X> A, const X &B)
{
  A /= B;

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> operator+ (matrix<X> A, const X &B)
{
  A += B;

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> operator- (matrix<X> A, const X &B)
{
  A -= B;

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> operator* (const X &A, matrix<X> B)
{
  B *= A;

  return B;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> operator/ (const X &A, const matrix<X> &B)
{
  matrix<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> operator+ (const X &A, matrix<X> B)
{
  B += A;

  return B;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> operator- (const X &A, const matrix<X> &B)
{
  matrix<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

