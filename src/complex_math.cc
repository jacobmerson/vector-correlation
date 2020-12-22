#include "complex_math.h"
namespace VectorCorrelation
{
  KOKKOS_FUNCTION
  auto conjugate_times(CST a, CST b) -> CST
  {
    return CST{a.real() * b.real() + a.imag() * b.imag(),
               a.real() * b.imag() - b.real() * a.imag()};
  }
}
