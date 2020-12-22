#ifndef COMMON_TYPES_H__
#define COMMON_TYPES_H__
#include "Kokkos_Core.hpp"

namespace VectorCorrelation
{
using OrdinalType = int;
using ScalarType = double;
using ComplexScalarType = Kokkos::complex<ScalarType>;

// images are 2d views of data
using ComplexImageType = Kokkos::View<ComplexScalarType **>;
using ScalarImageType = Kokkos::View<ScalarType **>;
using OrdinalImageType = Kokkos::View<OrdinalType **>;
using ConstComplexImageType = typename ComplexImageType::const_type;
using ConstScalarImageType = typename ScalarImageType::const_type;
using ConstOrdinalImageType = typename OrdinalImageType::const_type;

// redefine Kokkos backend types, so end user doesn't need to know about Kokkos
using Serial = Kokkos::Serial;
using OpenMP = Kokkos::OpenMP;
}
#endif
