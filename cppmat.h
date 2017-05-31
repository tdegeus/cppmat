
#ifndef CPPMAT_H
#define CPPMAT_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

namespace mat {

// =================================================================================================
// mat::matrix
// =================================================================================================

template <class T> class matrix
{
private:

  std::vector<T>      _data;    // data array
  std::vector<size_t> _shape;   // number of entries in each dimensions
  std::vector<size_t> _strides; // stride length for each index
  // Note:
  // - Number of dimensions determined from "_shape.size()".
  // - "_strides.size()" can be bigger than "_shape.size()" to allow reading with more indices.
  //   E.g. after using "atleast_3d" also a 2D matrix can be read using 3 indices.

public:

  // (copy) constructor
  // ------------------

  matrix               (const matrix<T> &) = default;
  matrix<T>& operator= (const matrix<T> &) = default;
  matrix<T>(){};

  // explicit constructors
  // ---------------------

  matrix(std::vector<size_t> shape)
  { resize(shape); };

  matrix(std::vector<size_t> shape, T D)
  { resize(shape); for ( auto &i: _data ) i = D; };

  matrix(std::vector<size_t> shape, const T *D)
  { resize(shape); for ( size_t i=0; i<size(); ++i )_data[i] = D[i]; };

  // constructor to copy + change data type
  // --------------------------------------

  template<typename U,typename V=T,\
    typename=typename std::enable_if<std::is_convertible<T,U>::value>::type>
  operator matrix<U> ()
  {
    // - allocate copy
    matrix<U> out(shape());
    // - copy all items (+ change data-type)
    for ( size_t i=0; i<size(); ++i ) out[i] = static_cast<T>(_data[i]);
    // - return copy
    return out;
  };

  // resize matrix
  // -------------

  void resize(std::vector<size_t> shape)
  {
    if ( shape.size()<1 )
      throw std::runtime_error("Input should be >= 1-D");

    _shape  .resize(shape.size());
    _strides.resize(shape.size());
    size_t n = shape[0];

    for ( size_t i=1; i<shape.size(); ++i ) { n          *= shape[i]; }
    for ( size_t i=0; i<shape.size(); ++i ) { _shape  [i] = shape[i]; }
    for ( size_t i=0; i<shape.size(); ++i ) { _strides[i] = 1;        }

    for ( size_t i=0; i<shape.size(); ++i )
      for ( size_t j=i+1; j<shape.size(); ++j )
        _strides[i] *= _shape[j];

    _data.resize(n);
  };

  // convert matrix to view it in d >= nd
  // ------------------------------------

  void atleast_1d()         { while ( _strides.size()<1 ) { _strides.push_back(1); } }
  void atleast_2d()         { while ( _strides.size()<2 ) { _strides.push_back(1); } }
  void atleast_3d()         { while ( _strides.size()<3 ) { _strides.push_back(1); } }
  void atleast_nd(size_t n) { while ( _strides.size()<n ) { _strides.push_back(1); } }

  // index operators (now up to 6-d, extend if needed)
  // -------------------------------------------------

  T& operator[](size_t i)
  { return _data[i]; };

  T& operator()(size_t a)
  { return _data[a*_strides[0]]; };

  T& operator()(size_t a, size_t b)
  { return _data[a*_strides[0]+b*_strides[1]]; };

  T& operator()(size_t a, size_t b, size_t c)
  { return _data[a*_strides[0]+b*_strides[1]+c*_strides[2]]; };

  T& operator()(size_t a, size_t b, size_t c, size_t d)
  { return _data[a*_strides[0]+b*_strides[1]+c*_strides[2]+d*_strides[3]]; };

  T& operator()(size_t a, size_t b, size_t c, size_t d, size_t e)
  { return _data[a*_strides[0]+b*_strides[1]+c*_strides[2]+d*_strides[3]+e*_strides[4]]; };

  T& operator()(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f)
  { return _data[a*_strides[0]+b*_strides[1]+c*_strides[2]+d*_strides[3]+e*_strides[4]+f*_strides[5]]; };

  // arithmetic operators
  // --------------------

  matrix<T>& operator*= (const matrix<T> &rhs)
  { for ( size_t i=0; i<size(); ++i ) _data[i] *= rhs._data[i]; return *this; };

  matrix<T>& operator/= (const matrix<T> &rhs)
  { for ( size_t i=0; i<size(); ++i ) _data[i] /= rhs._data[i]; return *this; };

  matrix<T>& operator+= (const matrix<T> &rhs)
  { for ( size_t i=0; i<size(); ++i ) _data[i] += rhs._data[i]; return *this; };

  matrix<T>& operator-= (const matrix<T> &rhs)
  { for ( size_t i=0; i<size(); ++i ) _data[i] -= rhs._data[i]; return *this; };

  matrix<T>& operator*= (T rhs)
  { for ( auto &i: _data ) i *= rhs; return *this; };

  matrix<T>& operator/= (T rhs)
  { for ( auto &i: _data ) i /= rhs; return *this; };

  matrix<T>& operator+= (T rhs)
  { for ( auto &i: _data ) i += rhs; return *this; };

  matrix<T>& operator-= (T rhs)
  { for ( auto &i: _data ) i -= rhs; return *this; };

  // iterators / pointer
  // -------------------

  const T* data () const { return _data.data (); };
  auto     begin()       { return _data.begin(); };
  auto     end  ()       { return _data.end  (); };


  // return shape array [ndim]
  // -------------------------

  std::vector<size_t> shape(size_t nd=0) const
  {
    if ( nd==0 )
      nd = ndim();

    std::vector<size_t> ret(nd,1);

    for ( size_t i=0; i<nd; ++i )
      ret[i] = _shape[i];

    return ret;
  };

  // return strides array [ndim]
  // ---------------------------

  std::vector<size_t> strides(bool bytes=false) const
  {
    size_t nd = ndim();
    std::vector<size_t> ret(nd);

    for ( size_t i=0; i<nd; ++i )
      ret[i] = _strides[i];

    if ( bytes )
      for ( size_t i=0; i<nd; ++i )
        ret[i] *= sizeof(T);

    return ret;
  };

  // return size
  // -----------

  size_t size() const { return _data .size(); };
  size_t ndim() const { return _shape.size(); };

  // minimum / maximum / mean / sum
  // ------------------------------

  T sum() const
  {
    T out = static_cast<T>(0);

    for ( auto i: _data )
      out += i;

    return out;
  };

  double mean() const { return static_cast<double>(this->sum())/static_cast<double>(size()); };
  T      min () const { return *std::min_element(_data.begin(),_data.end()); };
  T      max () const { return *std::max_element(_data.begin(),_data.end()); };

  // initialize to zero/one
  // ----------------------

  void zeros() { for ( auto &i: _data ) i = static_cast<T>(0); };
  void ones () { for ( auto &i: _data ) i = static_cast<T>(1); };

  // print to screen
  // ---------------

  void printf(std::string fmt) const
  {
    std::vector<size_t> s = _strides;
    if ( ndim()==1 ) {
      for ( size_t h=0; h<shape()[0]-1; ++h )
        std::printf((fmt+",").c_str(),_data[h]);
      std::printf((fmt+"\n").c_str(),_data[shape()[0]-1]);
    }
    else if ( ndim()==2 ) {
      for ( size_t h=0; h<shape()[0]; ++h ) {
        for ( size_t i=0; i<shape()[1]-1; ++i )
          std::printf((fmt+",").c_str(),_data[h*s[0]+i*s[1]]);
        std::printf((fmt+";\n").c_str(),_data[h*s[0]+(shape()[1]-1)*s[1]]);
      }
    }
    else if ( ndim()==3 ) {
      for ( size_t h=0; h<shape()[0]; ++h ) {
        for ( size_t i=0; i<shape()[1]; ++i ) {
          for ( size_t j=0; j<shape()[2]-1; ++j )
            std::printf((fmt+",").c_str(),_data[h*s[0]+i*s[1]+j*s[2]]);
          std::printf((fmt+";\n").c_str(),_data[h*s[0]+i*s[1]+(shape()[2]-1)*s[2]]);
        }
        if ( h<shape()[0]-1 )
          std::printf("\n");
      }
    }
  }

}; // class matrix

// arithmetic operators
// --------------------

template<class T>
matrix<T> operator* (const matrix<T> &A, const matrix<T> &B) { matrix<T> C = A; return C *= B; };

template<class T>
matrix<T> operator* (const matrix<T> &A, const        T  &B) { matrix<T> C = A; return C *= B; };

template<class T>
matrix<T> operator* (const        T  &A, const matrix<T> &B) { matrix<T> C = B; return C *= A; };

template<class T>
matrix<T> operator/ (const matrix<T> &A, const matrix<T> &B) { matrix<T> C = A; return C /= B; };

template<class T>
matrix<T> operator/ (const matrix<T> &A, const        T  &B) { matrix<T> C = A; return C /= B; };

template<class T>
matrix<T> operator/ (const        T  &A, const matrix<T> &B) { matrix<T> C = B; return C /= A; };

template<class T>
matrix<T> operator+ (const matrix<T> &A, const matrix<T> &B) { matrix<T> C = A; return C += B; };

template<class T>
matrix<T> operator+ (const matrix<T> &A, const        T  &B) { matrix<T> C = A; return C += B; };

template<class T>
matrix<T> operator+ (const        T  &A, const matrix<T> &B) { matrix<T> C = B; return C += A; };

template<class T>
matrix<T> operator- (const matrix<T> &A, const matrix<T> &B) { matrix<T> C = A; return C -= B; };

template<class T>
matrix<T> operator- (const matrix<T> &A, const        T  &B) { matrix<T> C = A; return C -= B; };

template<class T>
matrix<T> operator- (const        T  &A, const matrix<T> &B) { matrix<T> C = B; return C -= A; };

// print to "std::cout"
// --------------------

template <class T>
std::ostream& operator<<(std::ostream& out, matrix<T>& src)
{
  if ( src.ndim()==1 ) {
    for ( size_t i=0; i<src.shape()[0]-1; ++i )
      out << src(i) << " , ";
    out << src(src.shape()[0]-1) << std::endl;
  }
  else if ( src.ndim()==2 ) {
    for ( size_t i=0; i<src.shape()[0]; ++i ) {
      for ( size_t j=0; j<src.shape()[1]-1; ++j ) {
        out << src(i,j) << ", ";
      }
      out << src(i,src.shape()[1]-1) << "; " << std::endl;
    }
  }
  else if ( src.ndim()==3 ) {
    for ( size_t h=0; h<src.shape()[0]; ++h ) {
      for ( size_t i=0; i<src.shape()[1]; ++i ) {
        for ( size_t j=0; j<src.shape()[2]-1; ++j ) {
          out << src(h,i,j) << ", ";
        }
        out << src(i,src.shape()[1]-1) << "; " << std::endl;
      }
      if ( h<src.shape()[0]-1 )
        out << std::endl;
    }
  }

  return out;
};

}; // namespace mat

#endif
