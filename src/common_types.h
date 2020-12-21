#ifndef COMMON_TYPES_H__
#define COMMON_TYPES_H__
#include "Kokkos_Core.hpp"
using OrdinalType = int;
using ScalarType = double;
using ComplexScalarType = Kokkos::complex<ScalarType>;
// short alias for complex scalar
using CST = ComplexScalarType;
#endif
