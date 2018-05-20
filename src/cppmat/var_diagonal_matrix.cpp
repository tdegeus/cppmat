/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_DIAGONAL_MATRIX_CPP
#define CPPMAT_VAR_DIAGONAL_MATRIX_CPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace diagonal {

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline
matrix<X>::matrix(size_t m, size_t n)
{
  mZero[0] = static_cast<X>(0);

  assert( m == n );

  resize(m,n);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X>::matrix(const matrix<X> &A)
{
  mZero[0] = static_cast<X>(0);

  resize(A.shape(0), A.shape(1));

  setCopy(A.begin(), A.end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> matrix<X>::Random(size_t m, size_t n, X lower, X upper)
{
  matrix<X> out(m,n);

  out.setRandom(lower, upper);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> matrix<X>::Arange(size_t m, size_t n)
{
  matrix<X> out(m,n);

  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> matrix<X>::Zero(size_t m, size_t n)
{
  matrix<X> out(m,n);

  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> matrix<X>::Ones(size_t m, size_t n)
{
  matrix<X> out(m,n);

  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
matrix<X> matrix<X>::Constant(size_t m, size_t n, X D)
{
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
  matrix<X> out(m,n);

  out.setCopy(first,last);

  return out;
}

// =================================================================================================
// return plain storage as vector
// =================================================================================================

template<class X>
inline
std::vector<X> matrix<X>::asVector() const
{
  return mData;
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
  mSize = N;

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

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t matrix<X>::rows() const
{
  return N;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t matrix<X>::cols() const
{
  return N;
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

  if (a == b) return mData[a];
  else        return mZero[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
const X& matrix<X>::operator()(size_t a, size_t b) const
{
  assert( a < N );
  assert( b < N );

  if (a == b) return mData[a];
  else        return mZero[0];
}

// =================================================================================================
// index operators : compress(...)
// =================================================================================================

template<class X>
inline
size_t matrix<X>::compress(size_t a, size_t b) const
{
  assert( a < N  );
  assert( b < N  );
  assert( a == b );

  return a;
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
  std::fill(idx.begin(), idx.end(), i);

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
  assert( a < N  );
  assert( b < N  );
  assert( a == b );

  return this->begin() + a;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto matrix<X>::item(size_t a, size_t b) const
{
  assert( a < N  );
  assert( b < N  );
  assert( a == b );

  return this->begin() + a;
}

// =================================================================================================
// initialize
// =================================================================================================

template<class X>
inline
void matrix<X>::setRandom(X lower, X upper)
{
  // type of random number distribution
  std::uniform_real_distribution<X> dist(lower, upper);

  // Mersenne Twister: Good quality random number generator
  std::mt19937 rng;

  // Initialize with non-deterministic seeds
  rng.seed(std::random_device{}());

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] = dist(rng);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void matrix<X>::setArange()
{
  std::iota(mData.begin(), mData.end(), static_cast<X>(0));
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
  assert( mSize == static_cast<size_t>(last-first) );

  std::copy(first, last, begin());
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
void matrix<X>::setCopyDense(Iterator first)
{
  // check the input to be diagonal
  #ifndef NDEBUG
  for ( size_t i = 0 ; i < N ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      if ( i !=j )
        assert( !first[i*N+j] );
  #endif

  // copy from input (ignores off-diagonal terms)
  for ( size_t i = 0 ; i < N ; ++i )
    mData[i] = first[i*N+i];
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
  assert( N*N == static_cast<size_t>(last-first) );

  // check the input to be diagonal
  #ifndef NDEBUG
  for ( size_t i = 0 ; i < N ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      if ( i !=j )
        assert( !first[i*N+j] );
  #endif

  // copy from input (ignores off-diagonal terms)
  for ( size_t i = 0 ; i < N ; ++i )
    mData[i] = first[i*N+i];
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
  assert( mSize == static_cast<size_t>(last-first) );

  UNUSED(last);

  std::copy(begin(), end(), first);
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline
void matrix<X>::copyToDense(Iterator first) const
{
  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
	  if ( i == j ) first[i*N+j] = mData[i];
	  else          first[i*N+j] = static_cast<X>(0);
	}
  }
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline
void matrix<X>::copyToDense(Iterator first, Iterator last) const
{
  assert( N*N == static_cast<size_t>(last-first) );

  UNUSED(last);

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
	  if ( i == j ) first[i*N+j] = mData[i];
	  else          first[i*N+j] = static_cast<X>(0);
	}
  }
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

// =================================================================================================
// absolute value
// =================================================================================================

template<class X>
inline
matrix<X> matrix<X>::abs() const
{
  matrix<X> out(N, N);

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
// return the indices that would sort the matrix
// =================================================================================================

template<class X>
inline
matrix<size_t> matrix<X>::argsort(bool ascending) const
{
  return matrix<size_t>(N, N, cppmat::argsort(mData, ascending));
}

// =================================================================================================
// location of the minimum/maximum
// =================================================================================================

template<class X>
inline
size_t matrix<X>::argmin() const
{
  return std::min_element(begin(), end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t matrix<X>::argmax() const
{
  return std::max_element(begin(), end()) - begin();
}

// =================================================================================================
// minimum
// =================================================================================================

template<class X>
inline
X matrix<X>::min() const
{
  return *std::min_element(begin(),end());
}

// =================================================================================================
// maximum
// =================================================================================================

template<class X>
inline
X matrix<X>::max() const
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
  return std::accumulate(begin(), end(), static_cast<X>(0));
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
// find the plain storage indices of all non-zero entries
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
// find the plain storage indices of all entries equal to some constant
// =================================================================================================

template<class X>
inline
std::vector<size_t> matrix<X>::where(X D) const
{
  size_t nnz = 0;

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] == D )
      ++nnz;

  std::vector<size_t> out(nnz);

  size_t j = 0;

  for ( size_t i = 0 ; i < mSize ; ++i ) {
    if ( mData[i] == D ) {
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
matrix<X> operator* (const X &A, matrix<X> B)
{
  B *= A;

  return B;
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

