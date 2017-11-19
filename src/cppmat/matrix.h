/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MATRIX_H
#define CPPMAT_MATRIX_H

#include "macros.h"

namespace cppmat {

#define MAX_DIM 6

// =================================================================================================
// cppmat::matrix
// =================================================================================================

template <class X>
class matrix
{
private:

  std::vector<X> m_data;             // data container
  size_t         m_ndim=0;           // actual number of dimensions
  size_t         m_size=0;           // total number of entries == data.size() == prod(shape)
  size_t         m_shape[MAX_DIM];   // number of entries in each dimensions
  size_t         m_strides[MAX_DIM]; // stride length for each index

public:

  // constructors
  matrix(){};
  matrix(const std::vector<size_t> &shape);
  matrix(const std::vector<size_t> &shape, X D);
  matrix(const std::vector<size_t> &shape, const X *D);

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  size_t shape(size_t i) const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // resize matrix
  void resize (const std::vector<size_t> &shape);
  void reshape(const std::vector<size_t> &shape);
  void chdim  (size_t ndim);

  // index operators: access plain storage
  X&       operator[](size_t i);
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  X&       operator()(size_t a);
  const X& operator()(size_t a) const;
  X&       operator()(size_t a, size_t b);
  const X& operator()(size_t a, size_t b) const;
  X&       operator()(size_t a, size_t b, size_t c);
  const X& operator()(size_t a, size_t b, size_t c) const;
  X&       operator()(size_t a, size_t b, size_t c, size_t d);
  const X& operator()(size_t a, size_t b, size_t c, size_t d) const;
  X&       operator()(size_t a, size_t b, size_t c, size_t d, size_t e);
  const X& operator()(size_t a, size_t b, size_t c, size_t d, size_t e) const;
  X&       operator()(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f);
  const X& operator()(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f) const;

  // arithmetic operators
  matrix<X>& operator*= (const matrix<X> &B);
  matrix<X>& operator/= (const matrix<X> &B);
  matrix<X>& operator+= (const matrix<X> &B);
  matrix<X>& operator-= (const matrix<X> &B);
  matrix<X>& operator*= (X B);
  matrix<X>& operator/= (X B);
  matrix<X>& operator+= (X B);
  matrix<X>& operator-= (X B);

  // pointer / iterators
  X*       data();
  const X* data() const;
  auto     begin();
  auto     begin() const;
  auto     end();
  auto     end() const;

  // conversion operators
  template<typename U,typename V=X,typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator matrix<U> () const;

  template<typename U,typename V=X,typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator std::vector<U> () const;

  #ifdef CPPMAT_EIGEN
  template<typename U,typename V=X,typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator Eigen::Matrix<U,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> () const;
  #endif

  #ifdef CPPMAT_EIGEN
  template<typename U,typename V=X,typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator Eigen::Matrix<U,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> () const;
  #endif

  // basic algebra
  X              sum() const;
  X              min() const;
  X              max() const;
  double         mean() const;
  double         average(const matrix<X> &weights) const;
  matrix<double> average(const matrix<X> &weights, size_t axis) const;

  // basic initialization
  void setConstant(X D);
  void setZero();
  void setOnes();
  void zeros();
  void ones();

  // formatted print; NB also "operator<<" is defined below
  void printf(std::string fmt) const;

}; // class matrix

// arithmetic operators
template <class X> matrix<X> operator* (const matrix<X> &A, const matrix<X> &B);
template <class X> matrix<X> operator/ (const matrix<X> &A, const matrix<X> &B);
template <class X> matrix<X> operator+ (const matrix<X> &A, const matrix<X> &B);
template <class X> matrix<X> operator- (const matrix<X> &A, const matrix<X> &B);
template <class X> matrix<X> operator* (const matrix<X> &A, const        X  &B);
template <class X> matrix<X> operator/ (const matrix<X> &A, const        X  &B);
template <class X> matrix<X> operator+ (const matrix<X> &A, const        X  &B);
template <class X> matrix<X> operator- (const matrix<X> &A, const        X  &B);
template <class X> matrix<X> operator* (const        X  &A, const matrix<X> &B);
template <class X> matrix<X> operator/ (const        X  &A, const matrix<X> &B);
template <class X> matrix<X> operator+ (const        X  &A, const matrix<X> &B);
template <class X> matrix<X> operator- (const        X  &A, const matrix<X> &B);

// =================================================================================================
// constructors
// =================================================================================================

template <class X>
matrix<X>::matrix(const std::vector<size_t> &shape)
{
  resize(shape);
}

template <class X>
matrix<X>::matrix(const std::vector<size_t> &shape, X D)
{
  resize(shape);

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] = D;
}

template <class X>
matrix<X>::matrix(const std::vector<size_t> &shape, const X *D)
{
  resize(shape);

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] = D[i];
}

// =================================================================================================
// get dimensions
// =================================================================================================

template <class X>
size_t matrix<X>::size() const
{
  return m_size;
}

// -------------------------------------------------------------------------------------------------

template <class X>
size_t matrix<X>::ndim() const
{
  return m_ndim;
}

// -------------------------------------------------------------------------------------------------

template <class X>
size_t matrix<X>::shape(size_t i) const
{
  assert( i < MAX_DIM );

  return m_shape[i];
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::vector<size_t> matrix<X>::shape() const
{
  std::vector<size_t> ret(m_ndim);

  for ( size_t i = 0 ; i < m_ndim ; ++i ) ret[i] = m_shape[i];

  return ret;
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::vector<size_t> matrix<X>::strides(bool bytes) const
{
  std::vector<size_t> ret(m_ndim);

  for ( size_t i = 0 ; i < m_ndim ; ++i )
    ret[i] = m_strides[i];

  if ( bytes )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      ret[i] *= sizeof(X);

  return ret;
}

// =================================================================================================
// resize matrix
// =================================================================================================

template <class X>
void matrix<X>::resize(const std::vector<size_t> &shape)
{
  assert( shape.size()  > 0       );
  assert( shape.size() <= MAX_DIM );

  for ( size_t i = 0 ; i < MAX_DIM ; ++i ) {
    m_shape  [i] = 1;
    m_strides[i] = 1;
  }

  m_ndim = shape.size();
  m_size = 1;

  for ( size_t i = 0 ; i < m_ndim ; ++i ) {
    m_shape[i] = shape[i];
    m_size    *= shape[i];
  }

  for ( size_t i = 0 ; i < m_ndim ; ++i )
    for ( size_t j = i+1 ; j < m_ndim ; ++j )
      m_strides[i] *= m_shape[j];

  m_data.resize(m_size);
}

// -------------------------------------------------------------------------------------------------

template <class X>
void matrix<X>::chdim(size_t ndim)
{
  #ifndef NDEBUG
    for ( size_t i = ndim ; i < MAX_DIM ; ++i ) assert( m_shape[i] == 1 );
  #endif

  m_ndim = ndim;
}

// -------------------------------------------------------------------------------------------------

template <class X>
void matrix<X>::reshape(const std::vector<size_t> &shape)
{
  #ifndef NDEBUG
    size_t n = 1;

    for ( auto &i : shape ) n *= i;

    assert( n == m_size );
  #endif

  resize(shape);
}

// =================================================================================================
// index operators
// =================================================================================================

template <class X>
X& matrix<X>::operator[](size_t i)
{
  return m_data[i];
}

// -------------------------------------------------------------------------------------------------

template <class X>
const X& matrix<X>::operator[](size_t i) const
{
  return m_data[i];
}

// -------------------------------------------------------------------------------------------------

template <class X>
X& matrix<X>::operator()(size_t a)
{
  return m_data[a*m_strides[0]];
}

// -------------------------------------------------------------------------------------------------

template <class X>
const X& matrix<X>::operator()(size_t a) const
{
  return m_data[a*m_strides[0]];
}

// -------------------------------------------------------------------------------------------------

template <class X>
X& matrix<X>::operator()(size_t a, size_t b)
{
  return m_data[a*m_strides[0]+b*m_strides[1]];
}

// -------------------------------------------------------------------------------------------------

template <class X>
const X& matrix<X>::operator()(size_t a, size_t b) const
{
  return m_data[a*m_strides[0]+b*m_strides[1]];
}

// -------------------------------------------------------------------------------------------------

template <class X>
X& matrix<X>::operator()(size_t a, size_t b, size_t c)
{
  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]];
}

// -------------------------------------------------------------------------------------------------

template <class X>
const X& matrix<X>::operator()(size_t a, size_t b, size_t c) const
{
  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]];
}

// -------------------------------------------------------------------------------------------------

template <class X>
X& matrix<X>::operator()(size_t a, size_t b, size_t c, size_t d)
{
  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]];
}

// -------------------------------------------------------------------------------------------------

template <class X>
const X& matrix<X>::operator()(size_t a, size_t b, size_t c, size_t d) const
{
  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]];
}

// -------------------------------------------------------------------------------------------------

template <class X>
X& matrix<X>::operator()(size_t a, size_t b, size_t c, size_t d, size_t e)
{
  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]];
}

// -------------------------------------------------------------------------------------------------

template <class X>
const X& matrix<X>::operator()(size_t a, size_t b, size_t c, size_t d, size_t e) const
{
  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]];
}

// -------------------------------------------------------------------------------------------------

template <class X>
X& matrix<X>::operator()(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f)
{
  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]+f*m_strides[5]];
}

// -------------------------------------------------------------------------------------------------

template <class X>
const X& matrix<X>::operator()(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f) const
{
  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]+f*m_strides[5]];
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template <class X>
matrix<X>& matrix<X>::operator*= (const matrix<X> &B)
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X>& matrix<X>::operator/= (const matrix<X> &B)
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X>& matrix<X>::operator+= (const matrix<X> &B)
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X>& matrix<X>::operator-= (const matrix<X> &B)
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X>& matrix<X>::operator*= (X B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X>& matrix<X>::operator/= (X B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X>& matrix<X>::operator+= (X B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X>& matrix<X>::operator-= (X B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X> operator* (const matrix<X> &A, const matrix<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  matrix<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X> operator/ (const matrix<X> &A, const matrix<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  matrix<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X> operator+ (const matrix<X> &A, const matrix<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  matrix<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X> operator- (const matrix<X> &A, const matrix<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  matrix<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X> operator* (const matrix<X> &A, const X &B)
{
  matrix<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X> operator/ (const matrix<X> &A, const X &B)
{
  matrix<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X> operator+ (const matrix<X> &A, const X &B)
{
  matrix<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X> operator- (const matrix<X> &A, const X &B)
{
  matrix<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X> operator* (const X &A, const matrix<X> &B)
{
  matrix<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X> operator/ (const X &A, const matrix<X> &B)
{
  matrix<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X> operator+ (const X &A, const matrix<X> &B)
{
  matrix<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<X> operator- (const X &A, const matrix<X> &B)
{
  matrix<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// pointers / iterators
// =================================================================================================

template<class X>
X* matrix<X>::data ()
{
  return m_data.data();
}

// -------------------------------------------------------------------------------------------------

template<class X>
auto matrix<X>::begin()
{
  return m_data.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
auto matrix<X>::end()
{
  return m_data.end();
}

// -------------------------------------------------------------------------------------------------

template<class X>
const X* matrix<X>::data() const
{
  return m_data.data();
}

// -------------------------------------------------------------------------------------------------

template<class X>
auto matrix<X>::begin() const
{
  return m_data.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
auto matrix<X>::end() const
{
  return m_data.end();
}

// =================================================================================================
// conversion operators
// =================================================================================================

template<class X>
template<typename U, typename V, typename E>
matrix<X>::operator matrix<U> () const
{
  matrix<U> out(shape());

  for ( size_t i = 0 ; i < size() ; ++i )
    out[i] = static_cast<U>( m_data[i] );

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename U, typename V, typename E>
matrix<X>::operator std::vector<U> () const
{
  std::vector<U> out(m_size);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out[i] = m_data[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

#ifdef CPPMAT_EIGEN
template<class X>
template<typename U, typename V, typename E>
matrix<X>::operator Eigen::Matrix<U,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> () const
{
  assert( m_ndim == 2 );

  Eigen::Matrix<U,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> out(m_shape[0],m_shape[1]);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out(i) = m_data[i];

  return out;
}
#endif

// -------------------------------------------------------------------------------------------------

#ifdef CPPMAT_EIGEN
template<class X>
template<typename U, typename V, typename E>
matrix<X>::operator Eigen::Matrix<U,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> () const
{
  assert( m_ndim == 2 );

  Eigen::Matrix<U,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> out(m_shape[0],m_shape[1]);

  for ( size_t i = 0 ; i < m_shape[0] ; ++i )
    for ( size_t j = 0 ; j < m_shape[1] ; ++j )
      out(i,j) = m_data[i*m_strides[0]+j];

  return out;
}
#endif

// =================================================================================================
// basic algebra
// =================================================================================================

template <class X>
X matrix<X>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out += m_data[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template <class X>
double matrix<X>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(m_size);
}

// -------------------------------------------------------------------------------------------------

template <class X>
X matrix<X>::min() const
{
  return *std::min_element(m_data.begin(),m_data.end());
}

// -------------------------------------------------------------------------------------------------

template <class X>
X matrix<X>::max() const
{
  return *std::max_element(m_data.begin(),m_data.end());
}

// -------------------------------------------------------------------------------------------------

template <class X>
double matrix<X>::average(const matrix<X> &weights) const
{
  assert( size() == weights.size() );
  assert( ndim() == weights.ndim() );

  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out += m_data[i] * weights[i];

  return static_cast<double>(out)/static_cast<double>(weights.sum());
}

// -------------------------------------------------------------------------------------------------

template <class X>
matrix<double> matrix<X>::average(const matrix<X> &weights, size_t axis) const
{
  assert( size() == weights.size() );
  assert( m_ndim == weights.ndim() );
  assert( axis   <  m_ndim         );
  assert( m_ndim == 2              ); // TODO: generalize
  assert( axis   == 1              ); // TODO: generalize


  matrix<X> out ({m_shape[0]}, static_cast<X>(0));
  matrix<X> norm({m_shape[0]}, static_cast<X>(0));

  for ( size_t i = 0 ; i < m_shape[0] ; ++i ) {
    for ( size_t j = 0 ; j < m_shape[1] ; ++j ) {
      out (i) += (*this)(i,j) * weights(i,j);
      norm(i) +=                weights(i,j);
    }
  }

  matrix<double> normlized({m_shape[0]});

  for ( size_t i = 0 ; i < m_shape[0] ; ++i )
    normlized[i] = static_cast<double>(out[i]) / static_cast<double>(norm[i]);

  return normlized;
}

// =================================================================================================
// basic initialization
// =================================================================================================

template <class X>
void matrix<X>::setConstant(X D)
{
  for ( size_t i=0; i<m_size; ++i )
    m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template <class X>
void matrix<X>::setZero()
{
  for ( size_t i=0; i<m_size; ++i )
    m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template <class X>
void matrix<X>::setOnes()
{
  for ( size_t i=0; i<m_size; ++i )
    m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template <class X>
void matrix<X>::zeros()
{
  for ( size_t i=0; i<m_size; ++i )
    m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template <class X>
void matrix<X>::ones()
{
  for ( size_t i=0; i<m_size; ++i )
    m_data[i] = static_cast<X>(1);
}

// =================================================================================================
// formatted print
// =================================================================================================

template <class X>
void matrix<X>::printf(std::string fmt) const
{
  std::vector<size_t> s = strides();

  if ( m_ndim == 1 )
  {
    for ( size_t h = 0 ; h < shape(0)-1 ; ++h ) {
      std::printf((fmt+",").c_str(),m_data[h]);
    }
    std::printf((fmt+"\n").c_str(),m_data[shape(0)-1]);
  }
  else if ( m_ndim == 2 )
  {
    for ( size_t h = 0 ; h < shape(0) ; ++h ) {
      for ( size_t i = 0 ; i < shape(1)-1 ; ++i ) {
        std::printf((fmt+",").c_str(),m_data[h*s[0]+i*s[1]]);
      }
      std::printf((fmt+";\n").c_str(),m_data[h*s[0]+(shape(1)-1)*s[1]]);
    }
  }
  else if ( m_ndim == 3 )
  {
    for ( size_t h = 0 ; h < shape(0) ; ++h ) {
      for ( size_t i = 0 ; i < shape(1) ; ++i ) {
        for ( size_t j = 0 ; j < shape(2)-1 ; ++j ) {
          std::printf((fmt+",").c_str(),m_data[h*s[0]+i*s[1]+j*s[2]]);
        }
        std::printf((fmt+";\n").c_str(),m_data[h*s[0]+i*s[1]+(shape(2)-1)*s[2]]);
      }
      if ( h<shape(0)-1 )
        std::printf("\n");
    }
  }
  else
  {
    std::cout << "cppmat::matrix[";

    for ( size_t i = 0 ; i < m_ndim-1 ; ++i )
      std::cout << m_shape[i] << ",";

    std::cout << m_shape[m_ndim-1] << "]" << std::endl;
  }
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, matrix<X>& src)
{
  if ( src.ndim() == 1 )
  {
    for ( size_t i = 0 ; i < src.shape(0)-1 ; ++i ) {
      out << src(i) << " , ";
    }
    out << src(src.shape(0)-1) << std::endl;
  }
  else if ( src.ndim() == 2 )
  {
    for ( size_t i = 0 ; i < src.shape(0) ; ++i ) {
      for ( size_t j = 0 ; j < src.shape(1)-1 ; ++j ) {
        out << src(i,j) << ", ";
      }
      out << src(i,src.shape(1)-1) << "; " << std::endl;
    }
  }
  else if ( src.ndim() == 3 )
  {
    for ( size_t h = 0 ; h < src.shape(0) ; ++h ) {
      for ( size_t i = 0 ; i < src.shape(1) ; ++i ) {
        for ( size_t j = 0 ; j < src.shape(2)-1 ; ++j ) {
          out << src(h,i,j) << ", ";
        }
        out << src(h,i,src.shape(2)-1) << "; " << std::endl;
      }
      if ( h < src.shape(0)-1 )
        out << std::endl;
    }
  }
  else
    {
      std::cout << "cppmat::matrix[";

      for ( size_t i = 0 ; i < src.ndim()-1 ; ++i )
        std::cout << src.shape(i) << ",";

      std::cout << src.shape(src.ndim()-1) << "]" << std::endl;
    }

  return out;
}

// =================================================================================================

} // namespace cppmat

#endif

