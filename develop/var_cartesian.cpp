
#include <catch/catch.hpp>

#define EQ(a,b) REQUIRE_THAT( (a), Catch::WithinAbs((b), 1.e-10) );

#define CPPMAT_NOCONVERT
// #include <cppmat/cppmat.h>
#include "../src/cppmat/cppmat.h"

#include <Eigen/Eigen>

#include "support.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatD;
typedef Eigen::Matrix<double, Eigen::Dynamic,              1, Eigen::ColMajor> ColD;

static const size_t M = 11;
static const size_t N = 9;
static const size_t O = 6;
static const size_t P = 5;

typedef cppmat::cartesian::tensor4 <double> T4;
typedef cppmat::cartesian::tensor2 <double> T2;
typedef cppmat::cartesian::tensor2s<double> T2s;
typedef cppmat::cartesian::tensor2d<double> T2d;
typedef cppmat::cartesian::vector  <double> V;

// =================================================================================================

TEST_CASE("cppmat::cartesian", "var_cartesian_*.h")
{

// =================================================================================================
// tensor
// =================================================================================================

SECTION( "array += array" )
{
  MatD a = MatD::Random(M,N);
  MatD b = MatD::Random(M,N);

  Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
  Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      a(i,j) += b(i,j);

  A += B;

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      EQ( a(i,j), A(i,j) );
}

// // -------------------------------------------------------------------------------------------------

// SECTION( "array -= array" )
// {
//   MatD a = MatD::Random(M,N);
//   MatD b = MatD::Random(M,N);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
//   Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       a(i,j) -= b(i,j);

//   A -= B;

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( a(i,j), A(i,j) );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "array *= array" )
// {
//   MatD a = MatD::Random(M,N);
//   MatD b = MatD::Random(M,N);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
//   Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       a(i,j) *= b(i,j);

//   A *= B;

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( a(i,j), A(i,j) );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "array /= array" )
// {
//   MatD a = MatD::Random(M,N);
//   MatD b = MatD::Random(M,N) + MatD::Ones(M,N);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
//   Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       a(i,j) /= b(i,j);

//   A /= B;

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( a(i,j), A(i,j) );
// }

// // =================================================================================================
// // arithmetic
// // =================================================================================================

// SECTION( "array += scalar" )
// {
//   MatD a = MatD::Random(M,N);
//   double b = a(0,0);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       a(i,j) += b;

//   A += b;

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( a(i,j), A(i,j) );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "array -= scalar" )
// {
//   MatD a = MatD::Random(M,N);
//   double b = a(0,0);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       a(i,j) -= b;

//   A -= b;

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( a(i,j), A(i,j) );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "array *= scalar" )
// {
//   MatD a = MatD::Random(M,N);
//   double b = a(0,0);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       a(i,j) *= b;

//   A *= b;

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( a(i,j), A(i,j) );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "array /= scalar" )
// {
//   MatD a = MatD::Random(M,N);
//   double b = a(0,0) + 1.0;

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       a(i,j) /= b;

//   A /= b;

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( a(i,j), A(i,j) );
// }

// // =================================================================================================
// // arithmetic
// // =================================================================================================

// SECTION( "array + array" )
// {
//   MatD a = MatD::Random(M,N);
//   MatD b = MatD::Random(M,N);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
//   Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

//   MatD c = MatD::Zero(M,N);

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       c(i,j) = a(i,j) + b(i,j);

//   Arr C = A + B;

//   REQUIRE( static_cast<size_t>(c.size()) == C.size() );

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( c(i,j), C(i,j) );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "array - array" )
// {
//   MatD a = MatD::Random(M,N);
//   MatD b = MatD::Random(M,N);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
//   Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

//   MatD c = MatD::Zero(M,N);

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       c(i,j) = a(i,j) - b(i,j);

//   Arr C = A - B;

//   REQUIRE( static_cast<size_t>(c.size()) == C.size() );

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( c(i,j), C(i,j) );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "array * array" )
// {
//   MatD a = MatD::Random(M,N);
//   MatD b = MatD::Random(M,N);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
//   Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

//   MatD c = MatD::Zero(M,N);

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       c(i,j) = a(i,j) * b(i,j);

//   Arr C = A * B;

//   REQUIRE( static_cast<size_t>(c.size()) == C.size() );

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( c(i,j), C(i,j) );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "array / array" )
// {
//   MatD a = MatD::Random(M,N);
//   MatD b = MatD::Random(M,N) + MatD::Ones(M,N);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
//   Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

//   MatD c = MatD::Zero(M,N);

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       c(i,j) = a(i,j) / b(i,j);

//   Arr C = A / B;

//   REQUIRE( static_cast<size_t>(c.size()) == C.size() );

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( c(i,j), C(i,j) );
// }

// // =================================================================================================
// // arithmetic
// // =================================================================================================

// SECTION( "array + scalar" )
// {
//   MatD a = MatD::Random(M,N);
//   double b = a(0,0);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

//   MatD c = MatD::Zero(M,N);

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       c(i,j) = a(i,j) + b;

//   Arr C = A + b;

//   REQUIRE( static_cast<size_t>(c.size()) == C.size() );

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( c(i,j), C(i,j) );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "array - scalar" )
// {
//   MatD a = MatD::Random(M,N);
//   double b = a(0,0);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

//   MatD c = MatD::Zero(M,N);

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       c(i,j) = a(i,j) - b;

//   Arr C = A - b;

//   REQUIRE( static_cast<size_t>(c.size()) == C.size() );

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( c(i,j), C(i,j) );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "array * scalar" )
// {
//   MatD a = MatD::Random(M,N);
//   double b = a(0,0);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

//   MatD c = MatD::Zero(M,N);

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       c(i,j) = a(i,j) * b;

//   Arr C = A * b;

//   REQUIRE( static_cast<size_t>(c.size()) == C.size() );

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( c(i,j), C(i,j) );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "array / scalar" )
// {
//   MatD a = MatD::Random(M,N);
//   double b = a(0,0) + 1.0;

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

//   MatD c = MatD::Zero(M,N);

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       c(i,j) = a(i,j) / b;

//   Arr C = A / b;

//   REQUIRE( static_cast<size_t>(c.size()) == C.size() );

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( c(i,j), C(i,j) );
// }

// // =================================================================================================
// // arithmetic
// // =================================================================================================

// SECTION( "scalar + array" )
// {
//   MatD b = MatD::Random(M,N);
//   double a = b(0,0);

//   Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

//   MatD c = MatD::Zero(M,N);

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       c(i,j) = a + b(i,j);

//   Arr C = a + B;

//   REQUIRE( static_cast<size_t>(c.size()) == C.size() );

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( c(i,j), C(i,j) );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "scalar - array" )
// {
//   MatD b = MatD::Random(M,N);
//   double a = b(0,0);

//   Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

//   MatD c = MatD::Zero(M,N);

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       c(i,j) = a - b(i,j);

//   Arr C = a - B;

//   REQUIRE( static_cast<size_t>(c.size()) == C.size() );

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( c(i,j), C(i,j) );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "scalar * array" )
// {
//   MatD b = MatD::Random(M,N);
//   double a = b(0,0);

//   Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

//   MatD c = MatD::Zero(M,N);

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       c(i,j) = a * b(i,j);

//   Arr C = a * B;

//   REQUIRE( static_cast<size_t>(c.size()) == C.size() );

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( c(i,j), C(i,j) );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "scalar / array" )
// {
//   MatD b = MatD::Random(M,N) + MatD::Ones(M,N);
//   double a = b(0,0);

//   Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

//   MatD c = MatD::Zero(M,N);

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       c(i,j) = a / b(i,j);

//   Arr C = a / B;

//   REQUIRE( static_cast<size_t>(c.size()) == C.size() );

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( c(i,j), C(i,j) );
// }

// // =================================================================================================
// // algebra - partial
// // =================================================================================================

// SECTION( "min" )
// {
//   ColD a = ColD::Random(M*N*O*P);

//   Arr A = Arr::Copy({M,N,O,P}, a.data(), a.data()+a.size());

//   Arr C = A.min({-1,-2});

//   Arr c = Arr::Constant({M,N}, A.max());

//   for ( size_t i = 0 ; i < M ; i++ )
//     for ( size_t j = 0 ; j < N ; j++ )
//       for ( size_t k = 0 ; k < O ; k++ )
//         for ( size_t l = 0 ; l < P ; l++ )
//           c(i,j) = std::min( c(i,j), A(i,j,k,l) );

//   REQUIRE( c.size() == C.size() );

//   for ( size_t i = 0 ; i < c.size() ; ++i ) EQ( c[i], C[i] );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "max" )
// {
//   ColD a = ColD::Random(M*N*O*P);

//   Arr A = Arr::Copy({M,N,O,P}, a.data(), a.data()+a.size());

//   Arr C = A.max({-1,-2});

//   Arr c = Arr::Constant({M,N}, A.min());

//   for ( size_t i = 0 ; i < M ; i++ )
//     for ( size_t j = 0 ; j < N ; j++ )
//       for ( size_t k = 0 ; k < O ; k++ )
//         for ( size_t l = 0 ; l < P ; l++ )
//           c(i,j) = std::max( c(i,j), A(i,j,k,l) );

//   REQUIRE( c.size() == C.size() );

//   for ( size_t i = 0 ; i < c.size() ; ++i ) EQ( c[i], C[i] );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "sum" )
// {
//   ColD a = ColD::Random(M*N*O*P);

//   Arr A = Arr::Copy({M,N,O,P}, a.data(), a.data()+a.size());

//   Arr C = A.sum({0,1});

//   Arr c = Arr::Zero({O,P});

//   for ( size_t i = 0 ; i < M ; i++ )
//     for ( size_t j = 0 ; j < N ; j++ )
//       for ( size_t k = 0 ; k < O ; k++ )
//         for ( size_t l = 0 ; l < P ; l++ )
//           c(k,l) += A(i,j,k,l);

//   REQUIRE( c.size() == C.size() );

//   for ( size_t i = 0 ; i < c.size() ; ++i ) EQ( c[i], C[i] );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "mean" )
// {
//   ColD a = ColD::Random(M*N*O*P);

//   Arr A = Arr::Copy({M,N,O,P}, a.data(), a.data()+a.size());

//   Arr C = A.mean({0,1});

//   Arr c = Arr::Zero({O,P});

//   for ( size_t i = 0 ; i < M ; i++ )
//     for ( size_t j = 0 ; j < N ; j++ )
//       for ( size_t k = 0 ; k < O ; k++ )
//         for ( size_t l = 0 ; l < P ; l++ )
//           c(k,l) += A(i,j,k,l);

//   c /= static_cast<double>(M*N);

//   REQUIRE( c.size() == C.size() );

//   for ( size_t i = 0 ; i < c.size() ; ++i ) EQ( c[i], C[i] );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "average" )
// {
//   ColD a = ColD::Random(M*N*O*P);
//   ColD b = ColD::Random(M*N*O*P);

//   Arr A = Arr::Copy({M,N,O,P}, a.data(), a.data()+a.size());
//   Arr B = Arr::Copy({M,N,O,P}, b.data(), b.data()+b.size());

//   Arr C = A.average(B,{-2,0});

//   Arr c = Arr::Zero({N,P});
//   Arr d = Arr::Zero({N,P});

//   for ( size_t i = 0 ; i < M ; i++ ) {
//     for ( size_t j = 0 ; j < N ; j++ ) {
//       for ( size_t k = 0 ; k < O ; k++ ) {
//         for ( size_t l = 0 ; l < P ; l++ ) {
//           c(j,l) += B(i,j,k,l) * A(i,j,k,l);
//           d(j,l) += B(i,j,k,l);
//         }
//       }
//     }
//   }

//   c /= d;

//   REQUIRE( c.size() == C.size() );

//   for ( size_t i = 0 ; i < c.size() ; ++i ) EQ( c[i], C[i] );
// }

// // =================================================================================================
// // algebra
// // =================================================================================================

// SECTION( "min" )
// {
//   MatD a = MatD::Random(M,N);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

//   double C = A.min();

//   double c = a.minCoeff();

//   EQ( c, C );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "max" )
// {
//   MatD a = MatD::Random(M,N);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

//   double C = A.max();

//   double c = a.maxCoeff();

//   EQ( c, C );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "sum" )
// {
//   MatD a = MatD::Random(M,N);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

//   double C = A.sum();

//   double c = 0.0;

//   for ( auto i = 0 ; i < a.rows() ; ++i )
//     for ( auto j = 0 ; j < a.cols() ; ++j )
//       c += a(i,j);

//   EQ( c, C );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "mean" )
// {
//   MatD a = MatD::Random(M,N);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

//   double C = A.mean();

//   double c = 0.0;

//   for ( auto i = 0 ; i < a.rows() ; ++i )
//     for ( auto j = 0 ; j < a.cols() ; ++j )
//       c += a(i,j);

//   c /= static_cast<double>(M*N);

//   EQ( c, C );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "average" )
// {
//   MatD a = MatD::Random(M,N);
//   MatD b = MatD::Random(M,N);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());
//   Arr B = Arr::Copy({M,N}, b.data(), b.data()+b.size());

//   double C = A.average(B);

//   double c = 0.0;
//   double d = 0.0;

//   for ( auto i = 0 ; i < a.rows() ; ++i ) {
//     for ( auto j = 0 ; j < a.cols() ; ++j ) {
//       c += b(i,j) * a(i,j);
//       d += b(i,j);
//     }
//   }

//   c /= d;

//   EQ( c, C );
// }

// // =================================================================================================
// // absolute value
// // =================================================================================================

// SECTION( "abs" )
// {
//   MatD a = MatD::Random(M,N) - MatD::Constant(M,N,.5);

//   Arr A = Arr::Copy({M,N}, a.data(), a.data()+a.size());

//   Arr C = A.abs();

//   MatD c = a.cwiseAbs();

//   for ( size_t i = 0 ; i < M ; ++i )
//     for ( size_t j = 0 ; j < N ; ++j )
//       EQ( c(i,j), C(i,j) );
// }

// // =================================================================================================
// // index operators
// // =================================================================================================

// SECTION( "at" )
// {
//   ColD a = ColD::Random(M*N*O*P);

//   Arr A = Arr::Copy({M,N,O,P}, a.data(), a.data()+a.size());

//   std::vector<size_t> idx = {1,2,3,4};

//   EQ( A.at(idx.begin(), idx.end()), A(1,2,3,4) );
// }

// // -------------------------------------------------------------------------------------------------

// SECTION( "decompress" )
// {
//   ColD a = ColD::Random(M*N*O*P);

//   Arr A = Arr::Copy({M,N,O,P}, a.data(), a.data()+a.size());

//   std::vector<size_t> idx = A.decompress(A.compress(1,2,3,4));
//   std::vector<size_t> jdx = {1,2,3,4};

//   REQUIRE( idx == jdx );
// }


// =================================================================================================

}
