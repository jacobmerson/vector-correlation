#include "image_magnitude.h"
#include "complex_math.h"
#include <cmath>

namespace VectorCorrelation
{

template <typename Backend, typename ImageType>
void ImageMagnitude<Backend, ImageType>::Run()
{
  if(image_.extent(0) <1 || image_.extent(1) <1)
  {
    throw std::runtime_error("ImageMean filter must have an image assigned.");
  }
  // all new kokkos arrays are zero initialized
  magnitude_ = ScalarImageType{"magnitude", image_.extent(0), image_.extent(1)};
  Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy({0,0}, {image_.extent(0), image_.extent(1)});
  Kokkos::parallel_for("compute magnitude", policy, KOKKOS_LAMBDA(const int i, const int j)
  {
      magnitude_(i,j) = magnitude(image_(i,j));
  });
}

template class ImageMagnitude<Serial, ComplexImageType>;
}
