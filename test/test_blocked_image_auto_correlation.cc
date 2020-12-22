#include "catch2/catch_test_macros.hpp"
#include "catch2/catch_approx.hpp"
#include "common_types.h"
#include "blocked_image_auto_correlation.h"
#include "blocked_image_mean.h"
#include "complex_math.h"

using Catch::Approx;
using VectorCorrelation::ComplexScalarType;
using VectorCorrelation::ComplexImageType;
using VectorCorrelation::ScalarImageType;
using VectorCorrelation::BlockedImageMean;
using VectorCorrelation::BlockedImageAutoCorrelation;
using VectorCorrelation::Serial;
using VectorCorrelation::conjugate_square;


static void initialize_complex_image(ComplexImageType image)
{
  for (int i=0; i<image.extent(0); ++i)
  {
    for(int j=0; j<image.extent(1); ++j)
    {
      double val = 3*((i*image.extent(0)+j)%3);
      image(i,j) = ComplexScalarType{val, val};
    }
  }
}
TEST_CASE("auto_correlation", "[image]")
{
  ComplexImageType image{"auto correlation", 11, 11};
  initialize_complex_image(image);
  int block_size = 1;
  BlockedImageMean<Serial, ComplexImageType> blocked_image_mean;
  BlockedImageAutoCorrelation<Serial, ComplexImageType> blocked_image_auto_correlation;
  blocked_image_mean.SetInput(image);

  blocked_image_mean.SetBlockSize(block_size);
  blocked_image_mean.Run();
  auto mean = blocked_image_mean.GetOutput();

  blocked_image_auto_correlation.SetInput(image, mean);
  blocked_image_auto_correlation.SetBlockSize(block_size);
  blocked_image_auto_correlation.Run();
  auto auto_correlation = blocked_image_auto_correlation.GetOutput();

  REQUIRE(auto_correlation.extent(0) == image.extent(0));
  REQUIRE(auto_correlation.extent(1) == image.extent(1));
  double result = 0.0;
  for(int i=0; i<auto_correlation.extent(0); ++i)
  {
    for(int j=0; j<auto_correlation.extent(1); ++j)
    {
      result += abs(auto_correlation(i,j));
    }
  }
  REQUIRE(result == Approx(0));

  block_size = 3;
  int buffer = block_size/2;
  blocked_image_mean.SetBlockSize(block_size);
  blocked_image_mean.Run();
  mean = blocked_image_mean.GetOutput();
  blocked_image_auto_correlation.SetInput(image, mean);
  blocked_image_auto_correlation.SetBlockSize(block_size);
  blocked_image_auto_correlation.Run();
  auto_correlation = blocked_image_auto_correlation.GetOutput();

  for(int i=0; i<auto_correlation.extent(0); ++i)
  {
    for(int j=0; j<auto_correlation.extent(1); ++j)
    {
      if (i<buffer || j<buffer ||
          i > auto_correlation.extent(0)-(1+buffer) ||
          j > auto_correlation.extent(1)-(1+buffer))
      {
        result += abs(auto_correlation(i,j));
      }
      else
      {
        result += abs(auto_correlation(i,j)-12.0);
      }
    }
  }
  REQUIRE(result == Approx(0));
}
