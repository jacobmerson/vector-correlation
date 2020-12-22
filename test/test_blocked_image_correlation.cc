#include "catch2/catch_test_macros.hpp"
#include "catch2/catch_approx.hpp"
#include "common_types.h"
#include "blocked_image_auto_correlation.h"
#include "blocked_image_correlation.h"
#include "blocked_image_mean.h"

using Catch::Approx;
using VectorCorrelation::ComplexScalarType;
using VectorCorrelation::ComplexImageType;
using VectorCorrelation::ScalarImageType;
using VectorCorrelation::BlockedImageMean;
using VectorCorrelation::BlockedImageAutoCorrelation;
using VectorCorrelation::BlockedImageCorrelation;
using VectorCorrelation::Serial;


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
TEST_CASE("correlation", "[image]")
{
  ComplexImageType image{"auto correlation", 11, 11};
  initialize_complex_image(image);
  int block_size = 1;
  BlockedImageMean<Serial, ComplexImageType> blocked_image_mean;
  BlockedImageAutoCorrelation<Serial, ComplexImageType> blocked_image_auto_correlation;
  BlockedImageCorrelation<Serial, ComplexImageType> blocked_image_correlation;


  block_size = 3;
  int buffer = block_size/2;
  blocked_image_mean.SetInput(image);
  blocked_image_mean.SetBlockSize(block_size);
  blocked_image_mean.Run();
  auto mean = blocked_image_mean.GetOutput();
  blocked_image_auto_correlation.SetInput(image, mean);
  blocked_image_auto_correlation.SetBlockSize(block_size);
  blocked_image_auto_correlation.Run();
  auto auto_correlation = blocked_image_auto_correlation.GetOutput();

  blocked_image_correlation.SetInput(image,mean,auto_correlation,image,mean,auto_correlation);
  blocked_image_correlation.SetBlockSize(block_size);
  blocked_image_correlation.Run();
  auto correlation = blocked_image_correlation.GetOutput();
  ComplexScalarType result;
  int sum = 0;
  for(int i=0; i<auto_correlation.extent(0); ++i)
  {
    for(int j=0; j<auto_correlation.extent(1); ++j)
    {
      if (i<buffer || j<buffer ||
          i > auto_correlation.extent(0)-(1+buffer) ||
          j > auto_correlation.extent(1)-(1+buffer))
      {
        result += abs(correlation(i,j));
      }
      else
      {
        result += abs(correlation(i,j));
        sum += 1;
      }
    }
  }
  REQUIRE(result.real() == Approx(sum));
  REQUIRE(result.imag() == Approx(0.0));
}
