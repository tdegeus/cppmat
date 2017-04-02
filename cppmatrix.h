
#ifndef CPPMATRIX_H
#define CPPMATRIX_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

// =============================================================================
// Image::Matrix class
// =============================================================================

template <class T> class Matrix
{

  private:

    std::vector<T>      _data;     // data array
    std::vector<size_t> _shape;    // number of entries in each dimensions
    std::vector<size_t> _strides;  // stride length for each index

  public:

    // (copy) constructor
    // ------------------

    Matrix               (const Matrix<T> &) = default;
    Matrix<T>& operator= (const Matrix<T> &) = default;

    Matrix<T>(){};

    // explicit constructors
    // ---------------------

    Matrix ( std::vector<size_t> shape )
    { resize(shape); };

    Matrix ( std::vector<size_t> shape, T value )
    {
      resize(shape);

      for ( auto &i : _data )
        i = value;
    };

    Matrix ( std::vector<size_t> shape, const T *data )
    {
      resize(shape);

      for ( size_t i=0 ; i<size() ; ++i )
        _data[i] = data[i];
    };

    // constructor to copy + change data type
    // --------------------------------------

    template<\
      typename U,\
      typename V=T,\
      typename=typename std::enable_if<std::is_convertible<T,U>::value>::type\
    >
    operator Matrix<U> ()
    {
      Matrix<U> out(shape());
      for ( size_t i=0 ; i<size() ; ++i ) {
        out[i] = static_cast<T>(_data[i]);
      }
      return out;
    };

    // resize matrix
    // -------------

    void resize ( std::vector<size_t> shape )
    {
      if ( shape.size()<1 )
        throw std::runtime_error("Input should be >= 1-D");

      _shape  .resize(shape.size());
      _strides.resize(shape.size());
      size_t n = shape[0];

      for ( size_t i=1 ; i<shape.size() ; ++i ) { n          *= shape[i]; }
      for ( size_t i=0 ; i<shape.size() ; ++i ) { _shape  [i] = shape[i]; }
      for ( size_t i=0 ; i<shape.size() ; ++i ) { _strides[i] = 1;        }

      for ( size_t i=0 ; i<shape.size() ; ++i )
        for ( size_t j=i+1 ; j<shape.size() ; ++j )
          _strides[i] *= _shape[j];

      _data.resize(n);
    };

    // index operators (now up to 6-d, extend if needed)
    // -------------------------------------------------

    T& operator[] ( size_t i )
    { return _data[i]; };

    T& operator() ( size_t a )
    { return _data[a*_strides[0]]; };

    T& operator() ( size_t a, size_t b )
    { return _data[a*_strides[0]+b*_strides[1]]; };

    T& operator() ( size_t a, size_t b, size_t c )
    { return _data[a*_strides[0]+b*_strides[1]+c*_strides[2]]; };

    T& operator() ( size_t a, size_t b, size_t c, size_t d )
    { return _data[a*_strides[0]+b*_strides[1]+c*_strides[2]+d*_strides[3]]; };

    T& operator() ( size_t a, size_t b, size_t c, size_t d, size_t e )
    { return _data[a*_strides[0]+b*_strides[1]+c*_strides[2]+d*_strides[3]+e*_strides[4]]; };

    T& operator() ( size_t a, size_t b, size_t c, size_t d, size_t e, size_t f )
    { return _data[a*_strides[0]+b*_strides[1]+c*_strides[2]+d*_strides[3]+e*_strides[4]+f*_strides[5]]; };

    // arithmetic operators
    // --------------------

    Matrix<T>& operator*= (const Matrix<T> &rhs)
    {
      if ( shape()!=rhs.shape() )
        throw std::runtime_error("Matrices must have the same shape");

      for ( size_t i=0 ; i<size() ; ++i )
        _data[i] *= rhs._data[i];

      return *this;
    }

    Matrix<T>& operator/= (const Matrix<T> &rhs)
    {
      if ( shape()!=rhs.shape() )
        throw std::runtime_error("Matrices must have the same shape");

      for ( size_t i=0 ; i<size() ; ++i )
        _data[i] /= rhs._data[i];

      return *this;
    }

    Matrix<T>& operator+= (const Matrix<T> &rhs)
    {
      if ( shape()!=rhs.shape() )
        throw std::runtime_error("Matrices must have the same shape");

      for ( size_t i=0 ; i<size() ; ++i )
        _data[i] += rhs._data[i];

      return *this;
    }

    Matrix<T>& operator-= (const Matrix<T> &rhs)
    {
      if ( shape()!=rhs.shape() )
        throw std::runtime_error("Matrices must have the same shape");

      for ( size_t i=0 ; i<size() ; ++i )
        _data[i] -= rhs._data[i];

      return *this;
    }

    Matrix<T>& operator*= (T rhs)
    { for ( size_t i=0 ; i<size() ; ++i ) _data[i] *= rhs; return *this; }

    Matrix<T>& operator/= (T rhs)
    { for ( size_t i=0 ; i<size() ; ++i ) _data[i] /= rhs; return *this; }

    Matrix<T>& operator+= (T rhs)
    { for ( size_t i=0 ; i<size() ; ++i ) _data[i] += rhs; return *this; }

    Matrix<T>& operator-= (T rhs)
    { for ( size_t i=0 ; i<size() ; ++i ) _data[i] -= rhs; return *this; }

    // iterators / pointer
    // -------------------

    const T* data ( void ) const
    { return _data.data(); };

    auto begin ( void )
    { return _data.begin(); }

    auto end ( void )
    { return _data.end(); }

    // return shape array [ndim]
    // -------------------------

    std::vector<size_t> shape ( size_t nd=0 ) const
    {
      if ( nd==0 )
        nd = ndim();

      std::vector<size_t> ret(nd);

      for ( size_t i=0 ; i<nd ; ++i )
        ret[i] = _shape[i];

      return ret;
    };

    // return strides array [ndim]
    // ---------------------------

    std::vector<size_t> strides ( bool bytes=false ) const
    {
      size_t nd = ndim();
      std::vector<size_t> ret(nd);

      for ( size_t i=0 ; i<nd ; ++i )
        ret[i] = _strides[i];

      if ( bytes )
        for ( size_t i=0 ; i<nd ; ++i )
          ret[i] *= sizeof(T);

      return ret;
    };

    // return size
    // -----------

    size_t size ( void ) const
    { return _data.size(); };

    size_t ndim ( void ) const
    { return _shape.size(); };

    // minimum / maximum / mean / sum
    // ------------------------------

    T min ( void ) const
    { return *std::min_element(_data.begin(),_data.end()); };

    T max ( void ) const
    { return *std::max_element(_data.begin(),_data.end()); };

    double mean ( void ) const
    {
      T out = static_cast<T>(0);
      for ( auto i : _data )
        out += i;

      return static_cast<double>(out)/static_cast<double>(size());
    }

    T sum ( void ) const
    {
      T out = static_cast<T>(0);
      for ( auto i : _data )
        out += i;

      return out;
    }

    // initialize to zero/one
    // ----------------------

    void zeros ( void )
    {
      for ( size_t i=0 ; i<size() ; ++i )
        _data[i] = static_cast<T>(0);
    }

    void ones ( void )
    {
      for ( size_t i=0 ; i<size() ; ++i )
        _data[i] = static_cast<T>(1);
    }

    // print to screen
    // ---------------

    void printf(std::string fmt) const
    {
      std::string sep;
      std::string str;

      if ( ndim()==1 ) {
        sep = ",";
        str = fmt+sep;
        for ( size_t h=0 ; h<shape()[0]-1 ; ++h )
          std::printf(str.c_str(),_data[h]);
        sep = "\n";
        str = fmt+sep;
        std::printf(str.c_str(),_data[shape()[0]-1]);
      }
      else if ( ndim()==2 ) {
        for ( size_t h=0 ; h<shape()[0] ; ++h ) {
          sep = ",";
          str = fmt+sep;
          for ( size_t i=0 ; i<shape()[1]-1 ; ++i )
            std::printf(str.c_str(),_data[h*_strides[0]+i*_strides[1]]);
          sep = ";\n";
          str = fmt+sep;
          std::printf(str.c_str(),_data[h*_strides[0]+(shape()[1]-1)*_strides[1]]);
        }
      }
      else if ( ndim()==3 ) {
        for ( size_t h=0 ; h<shape()[0] ; ++h ) {
          for ( size_t i=0 ; i<shape()[1] ; ++i ) {
            sep = ",";
            str = fmt+sep;
            for ( size_t j=0 ; j<shape()[2]-1 ; ++j )
              std::printf(str.c_str(),_data[h*_strides[0]+i*_strides[1]+j*_strides[2]]);
            sep = ";\n";
            str = fmt+sep;
            std::printf(str.c_str(),_data[h*_strides[0]+i*_strides[1]+(shape()[2]-1)*_strides[2]]);
          }
          if ( h<shape()[0]-1 )
            std::printf("\n");
        }
      }
    }

}; // class Matrix

// arithmetic operators
// --------------------

template <class T>
Matrix<T> operator* (const Matrix<T> &A, const Matrix<T> &B)
{ Matrix<T> C = A; return C *= B; }

template <class T>
Matrix<T> operator* (const Matrix<T> &A, T B)
{ Matrix<T> C = A; return C *= B; }

template <class T>
Matrix<T> operator* (T A, const Matrix<T> &B)
{ Matrix<T> C = B; return C *= A; }

template <class T>
Matrix<T> operator/ (const Matrix<T> &A, const Matrix<T> &B)
{ Matrix<T> C = A; return C /= B; }

template <class T>
Matrix<T> operator/ (const Matrix<T> &A, T B)
{ Matrix<T> C = A; return C /= B; }

template <class T>
Matrix<T> operator/ (T A, const Matrix<T> &B)
{ Matrix<T> C = B; return C /= A; }

template <class T>
Matrix<T> operator+ (const Matrix<T> &A, const Matrix<T> &B)
{ Matrix<T> C = A; return C += B; }

template <class T>
Matrix<T> operator+ (const Matrix<T> &A, T B)
{ Matrix<T> C = A; return C += B; }

template <class T>
Matrix<T> operator+ (T A, const Matrix<T> &B)
{ Matrix<T> C = B; return C += A; }

template <class T>
Matrix<T> operator- (const Matrix<T> &A, const Matrix<T> &B)
{ Matrix<T> C = A; return C -= B; }

template <class T>
Matrix<T> operator- (const Matrix<T> &A, T B)
{ Matrix<T> C = A; return C -= B; }

template <class T>
Matrix<T> operator- (T A, const Matrix<T> &B)
{ Matrix<T> C = B; return C -= A; }

// print to "std::cout"
// --------------------

template <class T>
std::ostream& operator<<(std::ostream& out, Matrix<T>& src)
{
  if ( src.ndim()==1 ) {
    for ( size_t i=0 ; i<src.shape()[0]-1 ; ++i )
      out << src(i) << " , ";
    out << src(src.shape()[0]-1) << std::endl;
  }
  else if ( src.ndim()==2 ) {
    for ( size_t i=0 ; i<src.shape()[0] ; ++i ) {
      for ( size_t j=0 ; j<src.shape()[1]-1 ; ++j ) {
        out << src(i,j) << ", ";
      }
      out << src(i,src.shape()[1]-1) << "; " << std::endl;
    }
  }
  else if ( src.ndim()==3 ) {
    for ( size_t h=0 ; h<src.shape()[0] ; ++h ) {
      for ( size_t i=0 ; i<src.shape()[1] ; ++i ) {
        for ( size_t j=0 ; j<src.shape()[2]-1 ; ++j ) {
          out << src(h,i,j) << ", ";
        }
        out << src(i,src.shape()[1]-1) << "; " << std::endl;
      }
      if ( h<src.shape()[0]-1 )
        out << std::endl;
    }
  }

  return out;
}

#endif
