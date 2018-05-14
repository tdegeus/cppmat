
#include <vector>
#include <iostream>

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
class Matrix
{
protected:
  X mData[M*N];

public:
  Matrix() = default;
};

// -------------------------------------------------------------------------------------------------

template<class X, size_t M=2, size_t N=2>
class Matrix2d : public Matrix<X,M,N>
{
public:
  Matrix2d() = default;
};

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
Matrix<X,M,N> operator+ (Matrix<X,M,N> A, const Matrix<X,M,N> &B)
{
  std::cout << "Arbitrary size" << std::endl;

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X>
Matrix<X,2,2> operator+ (Matrix<X,2,2> A, const Matrix<X,2,2> &B)
{
  std::cout << "2-d" << std::endl;

  return A;
}

// -------------------------------------------------------------------------------------------------

int main()
{
  Matrix<double,2,2> A;

  Matrix2d<double> B;

  A + B;

  return 0;
}
