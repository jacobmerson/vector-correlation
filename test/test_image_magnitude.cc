#include "catch2/catch_test_macros.hpp"
#include "catch2/catch_approx.hpp"
#include "common_types.h"
#include "image_magnitude.h"
#include "complex_math.h"

using Catch::Approx;
using CST = VectorCorrelation::ComplexScalarType;
using VectorCorrelation::ComplexImageType;
using VectorCorrelation::ImageMagnitude;
using VectorCorrelation::Serial;
using VectorCorrelation::magnitude;

static void initialize_complex_image(ComplexImageType image)
{
  for (int i=0; i<image.extent(0); ++i)
  {
    for(int j=0; j<image.extent(1); ++j)
    {
      double val = i*image.extent(0)+j;
      image(i,j) = CST{val, val};
    }
  }
}

TEST_CASE("image_magnitude", "[image]")
{
  ComplexImageType image{"image mag test", 10, 10};
  initialize_complex_image(image);
  ImageMagnitude<Serial,ComplexImageType> image_magnitude;
  image_magnitude.SetInput(image);
  image_magnitude.Run();
  auto magnitude_image = image_magnitude.GetOutput();
  REQUIRE(magnitude_image.extent(0) == image.extent(0));
  REQUIRE(magnitude_image.extent(1) == image.extent(1));
  double result = 0.0;
  for(int i=0; i<image.extent(0); ++i)
  {
    for(int j=0; j<image.extent(1); ++j)
    {
      result += (magnitude(image(i,j)) - image_magnitude.GetOutput()(i,j));
    }
  }
  REQUIRE(result == Approx(0));
}
