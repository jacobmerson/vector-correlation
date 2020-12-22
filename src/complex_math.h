#ifndef VC_COMPLEX_MATH_H__
#define VC_COMPLEX_MATH_H__
#include "common_types.h"
#include "Kokkos_Core.hpp"
namespace VectorCorrelation
{
  KOKKOS_FUNCTION
  auto conjugate_times(CST a, CST b) -> CST;
}
#endif
