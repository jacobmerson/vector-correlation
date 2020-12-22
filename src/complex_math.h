#ifndef VC_COMPLEX_MATH_H__
#define VC_COMPLEX_MATH_H__
#include "common_types.h"
#include "Kokkos_Core.hpp"
namespace VectorCorrelation
{
  //KOKKOS_FUNCTION
  //auto conjugate_times(ComplexScalarType a, ComplexScalarType b) -> ComplexScalarType;
  template<typename T>
  KOKKOS_INLINE_FUNCTION
  auto conjugate_times(T a, T b) -> T { return a*b; }

  template<>
  KOKKOS_INLINE_FUNCTION
  auto conjugate_times (ComplexScalarType a, ComplexScalarType b) ->ComplexScalarType
  {
    return ComplexScalarType{a.real() * b.real() + a.imag() * b.imag(),
               a.real() * b.imag() - b.real() * a.imag()};
  }

  template <typename T>
  KOKKOS_INLINE_FUNCTION
  ScalarType conjugate_square(T val) { return val*val; }


  template<>
  KOKKOS_INLINE_FUNCTION
  ScalarType conjugate_square(ComplexScalarType val) { return conjugate_times(val,val).real(); }

  template <typename T>
  KOKKOS_INLINE_FUNCTION
  ScalarType magnitude(T val) { return fabs(val); }

  template <>
  KOKKOS_INLINE_FUNCTION
  ScalarType magnitude(ComplexScalarType val) { return sqrt(conjugate_square(val)); }
}
#endif
