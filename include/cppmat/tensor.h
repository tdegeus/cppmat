
#ifndef TENSOR_H
#define TENSOR_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

namespace cppmat {

// =================================================================================================
// forward declaration
// =================================================================================================

template<class X> class tensor4;
template<class X> class tensor2;
template<class X> class vector ;

// =================================================================================================
// return strides (generic routine used by all tensor-classes)
// =================================================================================================

template<class X> std::vector<size_t> inline _strides(size_t rank, size_t ndim, bool bytes=false)
{
  std::vector<size_t> out(rank,1);

  for ( size_t i=0; i<rank; ++i )
    for ( size_t j=i+1; j<rank; ++j )
      out[i] *= ndim;

  if ( bytes )
    for ( auto &i: out )
      i *= sizeof(X);

  return out;
};

// =================================================================================================
// tensor::tensor4
// =================================================================================================

template<class X> class tensor4
{
private:

  std::vector<X> m_data; // data array
  size_t         m_nd  ; // number of dimensions

public:

  // constructors
  // ------------

  // implicit constructor
  tensor4(){};
  // explicit constructor
  tensor4(size_t nd            ) { resize(nd);                                                   };
  tensor4(size_t nd,       X  D) { resize(nd); for ( size_t i=0; i<size(); ++i ) m_data[i]=D;    };
  tensor4(size_t nd, const X *D) { resize(nd); for ( size_t i=0; i<size(); ++i ) m_data[i]=D[i]; };

  // resize: change number of dimensions
  void resize(size_t nd) { m_nd = nd; m_data.resize(nd*nd*nd*nd); };

  // return strides array
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(4,m_nd,bytes); };

  // print to screen
  void printf() const;

  // copy constructor + change data-type
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor4<U> ()
  {
    // - allocate copy
    tensor4<U> out(m_nd);
    // - copy all items (+ change data-type)
    for ( size_t i=0; i<size(); ++i ) out[i] = static_cast<X>(m_data[i]);
    // - return copy
    return out;
  };

  // iterators / pointer
  // -------------------

  const X* data () const { return m_data.data (); };
  auto     begin()       { return m_data.begin(); };
  auto     end  ()       { return m_data.end  (); };

  // dimensions
  // ----------

  size_t size() const { return m_data.size(); };
  size_t ndim() const { return m_nd;          };

  // initialize to zero/one
  // ----------------------

  void zeros() { for ( auto &i: m_data ) i = static_cast<X>(0); };
  void ones () { for ( auto &i: m_data ) i = static_cast<X>(1); };

  // tensor products
  // ---------------

  tensor4<X> inline ddot(const tensor4<X> &B); // double contraction: C_ijmn = A_ijkl * B_lkmn
  tensor2<X> inline ddot(const tensor2<X> &B); // double contraction: C_ij   = A_ijkl * B_lk
  tensor4<X> inline T   () const;              // transposition     : B_lkji = A_ijkl
  tensor4<X> inline RT  () const;              // transposition     : B_ijlk = A_ijkl
  tensor4<X> inline LT  () const;              // transposition     : B_jikl = A_ijkl

  // operators
  // ---------

  X& operator[](size_t i)
  { return m_data[i]; };

  const X& operator[](size_t i) const
  { return m_data[i]; };

  X& operator()(size_t i, size_t j, size_t k, size_t l)
  { return m_data[i*m_nd*m_nd*m_nd+j*m_nd*m_nd+k*m_nd+l]; };

  const X& operator()(size_t i, size_t j, size_t k, size_t l) const
  { return m_data[i*m_nd*m_nd*m_nd+j*m_nd*m_nd+k*m_nd+l]; };

  tensor4<X>& operator*= (const tensor4<X> &rhs)
  { for ( size_t i=0; i<size(); ++i ) m_data[i] *= rhs.m_data[i]; return *this; };

  tensor4<X>& operator/= (const tensor4<X> &rhs)
  { for ( size_t i=0; i<size(); ++i ) m_data[i] /= rhs.m_data[i]; return *this; };

  tensor4<X>& operator+= (const tensor4<X> &rhs)
  { for ( size_t i=0; i<size(); ++i ) m_data[i] += rhs.m_data[i]; return *this; };

  tensor4<X>& operator-= (const tensor4<X> &rhs)
  { for ( size_t i=0; i<size(); ++i ) m_data[i] -= rhs.m_data[i]; return *this; };

  tensor4<X>& operator*= (const X &rhs)
  { for ( auto &i: m_data ) i *= rhs; return *this; };

  tensor4<X>& operator/= (const X &rhs)
  { for ( auto &i: m_data ) i /= rhs; return *this; };

  tensor4<X>& operator+= (const X &rhs)
  { for ( auto &i: m_data ) i += rhs; return *this; };

  tensor4<X>& operator-= (const X &rhs)
  { for ( auto &i: m_data ) i -= rhs; return *this; };

}; // class tensor4

template <class X>
tensor4<X> operator* (const tensor4<X> &A, const tensor4<X> &B) { tensor4<X> C=A; return C*=B; };

template <class X>
tensor4<X> operator* (const tensor4<X> &A, const         X  &B) { tensor4<X> C=A; return C*=B; };

template <class X>
tensor4<X> operator* (const         X  &A, const tensor4<X> &B) { tensor4<X> C=B; return C*=A; };

template <class X>
tensor4<X> operator/ (const tensor4<X> &A, const tensor4<X> &B) { tensor4<X> C=A; return C/=B; };

template <class X>
tensor4<X> operator/ (const tensor4<X> &A, const         X  &B) { tensor4<X> C=A; return C/=B; };

template <class X>
tensor4<X> operator/ (const         X  &A, const tensor4<X> &B) { tensor4<X> C=B; return C/=A; };

template <class X>
tensor4<X> operator+ (const tensor4<X> &A, const tensor4<X> &B) { tensor4<X> C=A; return C+=B; };

template <class X>
tensor4<X> operator+ (const tensor4<X> &A, const         X  &B) { tensor4<X> C=A; return C+=B; };

template <class X>
tensor4<X> operator+ (const         X  &A, const tensor4<X> &B) { tensor4<X> C=B; return C+=A; };

template <class X>
tensor4<X> operator- (const tensor4<X> &A, const tensor4<X> &B) { tensor4<X> C=A; return C-=B; };

template <class X>
tensor4<X> operator- (const tensor4<X> &A, const         X  &B) { tensor4<X> C=A; return C-=B; };

template <class X>
tensor4<X> operator- (const         X  &A, const tensor4<X> &B) { tensor4<X> C=B; return C-=A; };

// =================================================================================================
// tensor::tensor2
// =================================================================================================

template<class X> class tensor2
{
private:

  std::vector<X> m_data; // data array
  size_t         m_nd  ; // number of dimensions

public:

  // constructors
  // ------------

  // implicit constructor
  tensor2(){};
  // explicit constructor
  tensor2(size_t nd            ) { resize(nd);                                                   };
  tensor2(size_t nd,       X  D) { resize(nd); for ( size_t i=0; i<size(); ++i ) m_data[i]=D;    };
  tensor2(size_t nd, const X *D) { resize(nd); for ( size_t i=0; i<size(); ++i ) m_data[i]=D[i]; };

  // resize: change number of dimensions
  void resize(size_t nd) { m_nd = nd; m_data.resize(nd*nd); };

  // return strides array
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(2,m_nd,bytes); };

  // print to screen
  void printf(std::string fmt) const;

  // copy constructor + change data-type
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor2<U> ()
  {
    // - allocate copy
    tensor2<U> out(m_nd);
    // - copy all items (+ change data-type)
    for ( size_t i=0; i<size(); ++i ) out[i] = static_cast<X>(m_data[i]);
    // - return copy
    return out;
  };

  // iterators / pointer
  // -------------------

  const X* data () const { return m_data.data (); };
  auto     begin()       { return m_data.begin(); };
  auto     end  ()       { return m_data.end  (); };

  // dimensions
  // ----------

  size_t size() const { return m_data.size(); };
  size_t ndim() const { return m_nd;          };

  // initialize to zero/one
  // ----------------------

  void zeros() { for ( auto &i: m_data ) i = static_cast<X>(0); };
  void ones () { for ( auto &i: m_data ) i = static_cast<X>(1); };

  // tensor products
  // ---------------

  tensor2<X> inline dot   (const tensor2<X> &B); // single contraction: C_ik   = A_ij * B_jk
  vector <X> inline dot   (const vector <X> &B); // single contraction: C_i    = A_ij * B_j
  tensor2<X> inline ddot  (const tensor4<X> &B); // double contraction: C_kl   = A_ij * B_jikl
  X          inline ddot  (const tensor2<X> &B); // double contraction: C      = A_ij * B_ji
  tensor4<X> inline dyadic(const tensor2<X> &B); // dyadic product    : C_ijkl = A_ij * B_kl
  tensor2<X> inline T     () const;              // transpose         : B_ij   = A_ji
  X          inline trace () const;              // trace             : A_ii
  X          inline det   () const;              // determinant (only in 2D/3D)
  tensor2<X> inline inv   () const;              // inverse     (only in 2D/3D)

  // operators
  // ---------

  X&       operator[](size_t i          )       { return m_data[i];        };
  const X& operator[](size_t i          ) const { return m_data[i];        };
  X&       operator()(size_t i, size_t j)       { return m_data[i*m_nd+j]; };
  const X& operator()(size_t i, size_t j) const { return m_data[i*m_nd+j]; };

  tensor2<X>& operator*= (const tensor2<X> &rhs)
  { for ( size_t i=0; i<size(); ++i ) m_data[i] *= rhs.m_data[i]; return *this; };

  tensor2<X>& operator/= (const tensor2<X> &rhs)
  { for ( size_t i=0; i<size(); ++i ) m_data[i] /= rhs.m_data[i]; return *this; };

  tensor2<X>& operator+= (const tensor2<X> &rhs)
  { for ( size_t i=0; i<size(); ++i ) m_data[i] += rhs.m_data[i]; return *this; };

  tensor2<X>& operator-= (const tensor2<X> &rhs)
  { for ( size_t i=0; i<size(); ++i ) m_data[i] -= rhs.m_data[i]; return *this; };

  tensor2<X>& operator*= (const X &rhs)
  { for ( auto &i: m_data ) i *= rhs; return *this; };

  tensor2<X>& operator/= (const X &rhs)
  { for ( auto &i: m_data ) i /= rhs; return *this; };

  tensor2<X>& operator+= (const X &rhs)
  { for ( auto &i: m_data ) i += rhs; return *this; };

  tensor2<X>& operator-= (const X &rhs)
  { for ( auto &i: m_data ) i -= rhs; return *this; };

}; // class tensor2

// arithmetic operators
// --------------------

template <class X>
tensor2<X> operator* (const tensor2<X> &A, const tensor2<X> &B) { tensor2<X> C=A; return C*=B; };

template <class X>
tensor2<X> operator* (const tensor2<X> &A, const         X  &B) { tensor2<X> C=A; return C*=B; };

template <class X>
tensor2<X> operator* (const         X  &A, const tensor2<X> &B) { tensor2<X> C=B; return C*=A; };

template <class X>
tensor2<X> operator/ (const tensor2<X> &A, const tensor2<X> &B) { tensor2<X> C=A; return C/=B; };

template <class X>
tensor2<X> operator/ (const tensor2<X> &A, const         X  &B) { tensor2<X> C=A; return C/=B; };

template <class X>
tensor2<X> operator/ (const         X  &A, const tensor2<X> &B) { tensor2<X> C=B; return C/=A; };

template <class X>
tensor2<X> operator+ (const tensor2<X> &A, const tensor2<X> &B) { tensor2<X> C=A; return C+=B; };

template <class X>
tensor2<X> operator+ (const tensor2<X> &A, const         X  &B) { tensor2<X> C=A; return C+=B; };

template <class X>
tensor2<X> operator+ (const         X  &A, const tensor2<X> &B) { tensor2<X> C=B; return C+=A; };

template <class X>
tensor2<X> operator- (const tensor2<X> &A, const tensor2<X> &B) { tensor2<X> C=A; return C-=B; };

template <class X>
tensor2<X> operator- (const tensor2<X> &A, const         X  &B) { tensor2<X> C=A; return C-=B; };

template <class X>
tensor2<X> operator- (const         X  &A, const tensor2<X> &B) { tensor2<X> C=B; return C-=A; };

// =================================================================================================
// tensor::vector
// =================================================================================================

template<class X> class vector
{
private:

  std::vector<X> m_data; // data array
  size_t         m_nd  ; // number of dimensions

public:

  // constructors
  // ------------

  // implicit constructor
  vector(){};
  // explicit constructor
  vector(size_t nd            ) { resize(nd);                                                     };
  vector(size_t nd,       X  D) { resize(nd); for ( size_t i=0; i<size(); ++i ) m_data[i] = D;    };
  vector(size_t nd, const X *D) { resize(nd); for ( size_t i=0; i<size(); ++i ) m_data[i] = D[i]; };

  // resize: change number of dimensions
  void resize(size_t nd) { m_nd = nd; m_data.resize(nd); };

  // return strides array
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(1,m_nd,bytes); };

  // print to screen
  void printf(std::string fmt) const;

  // copy constructor + change data-type
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator vector<U> ()
  {
    // - allocate copy
    vector<U> out(m_nd);
    // - copy all items (+ change data-type)
    for ( size_t i=0; i<size(); ++i ) out[i] = static_cast<X>(m_data[i]);
    // - return copy
    return out;
  };

  // iterators / pointer
  // -------------------

  const X* data () const { return m_data.data (); };
  auto     begin()       { return m_data.begin(); };
  auto     end  ()       { return m_data.end  (); };

  // dimensions
  // ----------

  size_t size() const { return m_data.size(); };
  size_t ndim() const { return m_nd;          };

  // initialize to zero/one
  // ----------------------

  void zeros() { for ( auto &i: m_data ) i = static_cast<X>(0); };
  void ones () { for ( auto &i: m_data ) i = static_cast<X>(1); };

  // tensor products
  // ---------------

  tensor2<X> inline dyadic(const vector<X> &B); // dyadic product: C_ij = A_i*B_j
  vector <X> inline cross (const vector<X> &B); // cross product (only in 3D)

  // operators
  // ---------

  X&       operator[](size_t i)       { return m_data[i]; };
  const X& operator[](size_t i) const { return m_data[i]; };
  X&       operator()(size_t i)       { return m_data[i]; };
  const X& operator()(size_t i) const { return m_data[i]; };

  vector<X>& operator*= (const vector<X> &rhs)
  { for ( size_t i=0; i<size(); ++i ) m_data[i] *= rhs.m_data[i]; return *this; };

  vector<X>& operator/= (const vector<X> &rhs)
  { for ( size_t i=0; i<size(); ++i ) m_data[i] /= rhs.m_data[i]; return *this; };

  vector<X>& operator+= (const vector<X> &rhs)
  { for ( size_t i=0; i<size(); ++i ) m_data[i] += rhs.m_data[i]; return *this; };

  vector<X>& operator-= (const vector<X> &rhs)
  { for ( size_t i=0; i<size(); ++i ) m_data[i] -= rhs.m_data[i]; return *this; };

  vector<X>& operator*= (const X &rhs)
  { for ( auto &i: m_data ) i *= rhs; return *this; };

  vector<X>& operator/= (const X &rhs)
  { for ( auto &i: m_data ) i /= rhs; return *this; };

  vector<X>& operator+= (const X &rhs)
  { for ( auto &i: m_data ) i += rhs; return *this; };

  vector<X>& operator-= (const X &rhs)
  { for ( auto &i: m_data ) i -= rhs; return *this; };

}; // class vector

template <class X>
vector<X> operator* (const vector<X> &A, const vector<X> &B) { vector<X> C=A; return C*=B; };

template <class X>
vector<X> operator* (const vector<X> &A, const        X  &B) { vector<X> C=A; return C*=B; };

template <class X>
vector<X> operator* (const        X  &A, const vector<X> &B) { vector<X> C=B; return C*=A; };

template <class X>
vector<X> operator/ (const vector<X> &A, const vector<X> &B) { vector<X> C=A; return C/=B; };

template <class X>
vector<X> operator/ (const vector<X> &A, const        X  &B) { vector<X> C=A; return C/=B; };

template <class X>
vector<X> operator/ (const        X  &A, const vector<X> &B) { vector<X> C=B; return C/=A; };

template <class X>
vector<X> operator+ (const vector<X> &A, const vector<X> &B) { vector<X> C=A; return C+=B; };

template <class X>
vector<X> operator+ (const vector<X> &A, const        X  &B) { vector<X> C=A; return C+=B; };

template <class X>
vector<X> operator+ (const        X  &A, const vector<X> &B) { vector<X> C=B; return C+=A; };

template <class X>
vector<X> operator- (const vector<X> &A, const vector<X> &B) { vector<X> C=A; return C-=B; };

template <class X>
vector<X> operator- (const vector<X> &A, const        X  &B) { vector<X> C=A; return C-=B; };

template <class X>
vector<X> operator- (const        X  &A, const vector<X> &B) { vector<X> C=B; return C-=A; };

// =================================================================================================
// print to screen
// =================================================================================================

template<class X> void inline tensor4<X>::printf() const
{
  std::printf("tensor::tensor4 (%dD)\n",this->m_nd);
};

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, tensor4<X>& src)
{
  out << "tensor::tensor4 (" << src.ndim() << ")" << std::endl;
  return out;
};

// -------------------------------------------------------------------------------------------------

template<class X> void inline tensor2<X>::printf(std::string fmt) const
{
  size_t nd = this->m_nd;

  for ( size_t h=0; h<nd; ++h ) {
    for ( size_t i=0; i<nd-1; ++i )
      std::printf((fmt+",").c_str(),this->m_data[h*nd+i]);
    std::printf((fmt+";\n").c_str(),this->m_data[h*nd+(nd-1)]);
  }
};

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, tensor2<X>& src)
{
  for ( size_t i=0; i<src.ndim(); ++i ) {
    for ( size_t j=0; j<src.ndim()-1; ++j )
      out << src(i,j) << ", ";
    out << src(i,src.ndim()-1) << "; " << std::endl;
  }
  return out;
};

// -------------------------------------------------------------------------------------------------

template<class X> void inline vector<X>::printf(std::string fmt) const
{
  size_t nd = this->m_nd;

  for ( size_t i=0; i<nd-1; ++i )
    std::printf((fmt+",").c_str(),this->m_data[i]);
  std::printf((fmt+"\n").c_str(),this->m_data[nd-1]);
};

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, vector<X>& src)
{
  for ( size_t i=0; i<src.ndim()-1; ++i )
    out << src(i) << ", ";
  out << src(src.ndim()-1) << std::endl;
  return out;
};

// =================================================================================================
// identity tensors
// =================================================================================================

tensor4<double> inline identity4(size_t nd)
{
  tensor4<double> I(nd,0.0);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          if ( i==l and j==k )
            I(i,j,k,l) = 1.;

  return I;
};

// -------------------------------------------------------------------------------------------------

tensor4<double> inline identity4rt(size_t nd)
{
  tensor4<double> I(nd,0.0);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          if ( i==k and j==l )
            I(i,j,k,l) = 1.;

  return I;
};

// -------------------------------------------------------------------------------------------------

tensor4<double> inline identity4II(size_t nd)
{
  tensor4<double> I(nd,0.0);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          if ( i==j and k==l )
            I(i,j,k,l) = 1.;

  return I;
};

// -------------------------------------------------------------------------------------------------

tensor4<double> inline identity4s(size_t nd)
{ return (identity4(nd)+identity4rt(nd))/2.; };

// -------------------------------------------------------------------------------------------------

tensor4<double> inline identity4d(size_t nd)
{ return identity4s(nd)-identity4II(nd)/static_cast<double>(nd); };

// -------------------------------------------------------------------------------------------------

tensor2<double> inline identity2(size_t nd)
{
  tensor2<double> I(nd,0.0);

  for ( size_t i=0; i<nd; ++i )
    I(i,i) = 1.;

  return I;
};

// =================================================================================================
// tensor products
// =================================================================================================

template<class X> tensor4<X> inline tensor4<X>::ddot(const tensor4<X> &B)
{
  tensor4<X> C(this->ndim(),static_cast<X>(0));

  for ( size_t i=0; i<this->ndim(); ++i )
    for ( size_t j=0; j<this->ndim(); ++j )
      for ( size_t k=0; k<this->ndim(); ++k )
        for ( size_t l=0; l<this->ndim(); ++l )
          for ( size_t m=0; m<this->ndim(); ++m )
            for ( size_t n=0; n<this->ndim(); ++n )
              C(i,j,m,n) += (*this)(i,j,k,l)*B(l,k,m,n);

  return C;
};

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor4<X>::ddot(const tensor2<X> &B)
{
  tensor2<X> C(this->ndim(),static_cast<X>(0));

  for ( size_t i=0; i<this->ndim(); ++i )
    for ( size_t j=0; j<this->ndim(); ++j )
      for ( size_t k=0; k<this->ndim(); ++k )
        for ( size_t l=0; l<this->ndim(); ++l )
          C(i,j) += (*this)(i,j,k,l)*B(l,k);

  return C;
};

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2<X>::ddot(const tensor4<X> &B)
{
  tensor2<X> C(this->ndim(),static_cast<X>(0));

  for ( size_t i=0; i<this->ndim(); ++i )
    for ( size_t j=0; j<this->ndim(); ++j )
      for ( size_t k=0; k<this->ndim(); ++k )
        for ( size_t l=0; l<this->ndim(); ++l )
          C(k,l) += (*this)(i,j)*B(j,i,k,l);

  return C;
};

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2<X>::ddot(const tensor2<X> &B)
{
  X C = static_cast<X>(0);

  for ( size_t i=0; i<this->ndim(); ++i )
    for ( size_t j=0; j<this->ndim(); ++j )
      C += (*this)(i,j)*B(j,i);

  return C;
};

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2<X>::dot(const tensor2<X> &B)
{
  tensor2<X> C(this->ndim(),static_cast<X>(0));

  for ( size_t i=0; i<this->ndim(); ++i )
    for ( size_t j=0; j<this->ndim(); ++j )
      for ( size_t k=0; k<this->ndim(); ++k )
        C(i,k) += (*this)(i,j)*B(j,k);

  return C;
};

// -------------------------------------------------------------------------------------------------

template<class X> vector<X> inline tensor2<X>::dot(const vector<X> &B)
{
  vector<X> C(this->ndim(),static_cast<X>(0));

  for ( size_t i=0; i<this->ndim(); ++i )
    for ( size_t j=0; j<this->ndim(); ++j )
      C(i) += (*this)(i,j)*B(j);

  return C;
};

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2<X>::dyadic(const tensor2<X> &B)
{
  tensor4<X> C(this->ndim(),static_cast<X>(0));

  for ( size_t i=0; i<this->ndim(); ++i )
    for ( size_t j=0; j<this->ndim(); ++j )
      for ( size_t k=0; k<this->ndim(); ++k )
        for ( size_t l=0; l<this->ndim(); ++l )
          C(i,j,k,l) += (*this)(i,j)*B(k,l);

  return C;
};

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline vector<X>::dyadic(const vector<X> &B)
{
  tensor2<X> C(this->ndim(),static_cast<X>(0));

  for ( size_t i=0; i<this->ndim(); ++i )
    for ( size_t j=0; j<this->ndim(); ++j )
      C(i,j) += (*this)(i)*B(j);

  return C;
};

// -------------------------------------------------------------------------------------------------

template<class X> vector<X> inline vector<X>::cross(const vector<X> &B)
{
  if ( this->ndim()!=3 || B->ndim()!=3 )
    throw std::runtime_error("'cross' only implemented in 3D, use e.g. 'Eigen'");

  vector<X> C(3);

  C(0) =      (*this)(1)*B(2)-B(1)*(*this)(2) ;
  C(1) = -1.*((*this)(0)*B(2)-B(0)*(*this)(2));
  C(2) =      (*this)(0)*B(1)-B(0)*(*this)(1) ;

  return C;
};

// =================================================================================================
// transpositions
// =================================================================================================

template<class X> tensor4<X> inline tensor4<X>::T() const
{
  tensor4<X> C(this->ndim());

  for ( size_t i=0; i<this->ndim(); ++i )
    for ( size_t j=0; j<this->ndim(); ++j )
      for ( size_t k=0; k<this->ndim(); ++k )
        for ( size_t l=0; l<this->ndim(); ++l )
          C(l,k,j,i) = (*this)(i,j,k,l);

  return C;
};

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor4<X>::RT() const
{
  tensor4<X> C(this->ndim());

  for ( size_t i=0; i<this->ndim(); ++i )
    for ( size_t j=0; j<this->ndim(); ++j )
      for ( size_t k=0; k<this->ndim(); ++k )
        for ( size_t l=0; l<this->ndim(); ++l )
          C(i,j,l,k) = (*this)(i,j,k,l);

  return C;
};

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor4<X>::LT() const
{
  tensor4<X> C(this->ndim());

  for ( size_t i=0; i<this->ndim(); ++i )
    for ( size_t j=0; j<this->ndim(); ++j )
      for ( size_t k=0; k<this->ndim(); ++k )
        for ( size_t l=0; l<this->ndim(); ++l )
          C(j,i,k,l) = (*this)(i,j,k,l);

  return C;
};

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2<X>::T() const
{
  tensor2<X> C(this->ndim());

  for ( size_t i=0; i<this->ndim(); ++i )
    for ( size_t j=0; j<this->ndim(); ++j )
      C(j,i) = (*this)(i,j);

  return C;
};

// =================================================================================================
// miscellaneous
// =================================================================================================

template<class X> X inline tensor2<X>::trace() const
{
  X C = static_cast<X>(0);

  for ( size_t i=0; i<this->ndim(); ++i )
    C += (*this)(i,i);

  return C;
};

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2<X>::det() const
{
  if ( this->ndim()==2 )
   return (*this)(0,0)*(*this)(1,1)-(*this)(0,1)*(*this)(1,0);

  if ( this->ndim()==3 )
   return ((*this)(0,0)*(*this)(1,1)*(*this)(2,2)+
           (*this)(0,1)*(*this)(1,2)*(*this)(2,0)+
           (*this)(0,2)*(*this)(1,0)*(*this)(2,1))-
          ((*this)(0,2)*(*this)(1,1)*(*this)(2,0)+
           (*this)(0,1)*(*this)(1,0)*(*this)(2,2)+
           (*this)(0,0)*(*this)(1,2)*(*this)(2,1));

  throw std::runtime_error("'det' only implemented in 2D/3D, use e.g. 'Eigen'");
};

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2<X>::inv() const
{
  // compute determinant
  X D = this->det();

  // allocate result
  tensor2<X> C(this->ndim());

  if ( this->ndim()==2 ) {
    C(0,0) =     (*this)(1,1)/D;
    C(0,1) = -1.*(*this)(0,1)/D;
    C(1,0) = -1.*(*this)(1,0)/D;
    C(1,1) =     (*this)(0,0)/D;
    return C;
  }

  if ( this->ndim()==3 ) {
    C(0,0) = ((*this)(1,1)*(*this)(2,2)-(*this)(1,2)*(*this)(2,1))/D;
    C(0,1) = ((*this)(0,2)*(*this)(2,1)-(*this)(0,1)*(*this)(2,2))/D;
    C(0,2) = ((*this)(0,1)*(*this)(1,2)-(*this)(0,2)*(*this)(1,1))/D;
    C(1,0) = ((*this)(1,2)*(*this)(2,0)-(*this)(1,0)*(*this)(2,2))/D;
    C(1,1) = ((*this)(0,0)*(*this)(2,2)-(*this)(0,2)*(*this)(2,0))/D;
    C(1,2) = ((*this)(0,2)*(*this)(1,0)-(*this)(0,0)*(*this)(1,2))/D;
    C(2,0) = ((*this)(1,0)*(*this)(2,1)-(*this)(1,1)*(*this)(2,0))/D;
    C(2,1) = ((*this)(0,1)*(*this)(2,0)-(*this)(0,0)*(*this)(2,1))/D;
    C(2,2) = ((*this)(0,0)*(*this)(1,1)-(*this)(0,1)*(*this)(1,0))/D;
    return C;
  }

  throw std::runtime_error("'inv' only implemented in 2D/3D, use e.g. 'Eigen'");
};

// =================================================================================================
// create aliases to call class functions as functions, not members
// =================================================================================================

// products
template<class X> tensor4<X> inline ddot  (const tensor4<X> &A, const tensor4<X> &B)
{ return A.ddot  (B); };

template<class X> tensor2<X> inline ddot  (const tensor4<X> &A, const tensor2<X> &B)
{ return A.ddot  (B); };

template<class X> tensor2<X> inline ddot  (const tensor2<X> &A, const tensor4<X> &B)
{ return A.ddot  (B); };

template<class X>         X  inline ddot  (const tensor2<X> &A, const tensor2<X> &B)
{ return A.ddot  (B); };

template<class X> tensor2<X> inline dot   (const tensor2<X> &A, const tensor2<X> &B)
{ return A.dot   (B); };

template<class X> vector <X> inline dot   (const tensor2<X> &A, const vector <X> &B)
{ return A.dot   (B); };

template<class X> tensor4<X> inline dyadic(const tensor2<X> &A, const tensor2<X> &B)
{ return A.dyadic(B); };

template<class X> tensor2<X> inline dyadic(const vector <X> &A, const vector <X> &B)
{ return A.dyadic(B); };

template<class X> vector <X> inline cross (const vector <X> &A, const vector <X> &B)
{ return A.cross (B); };

// operations
template<class X> tensor2<X> inline transpose (const tensor2<X> &A) { return A.T    (); };
template<class X> tensor4<X> inline transpose (const tensor4<X> &A) { return A.T    (); };
template<class X> tensor4<X> inline transposeR(const tensor4<X> &A) { return A.RT   (); };
template<class X> tensor4<X> inline transposeL(const tensor4<X> &A) { return A.LT   (); };
template<class X> tensor2<X> inline inv       (const tensor2<X> &A) { return A.inv  (); };
template<class X>         X  inline det       (const tensor2<X> &A) { return A.det  (); };
template<class X>         X  inline trace     (const tensor2<X> &A) { return A.trace(); };

// =================================================================================================

}; // namespace tensor

#endif

