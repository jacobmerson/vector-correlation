#include "catch2/catch_test_macros.hpp"
#include "catch2/catch_approx.hpp"
#include "blocked_image_mean.h"

using Serial = VectorCorrelation::Serial;
using ComplexImageType = VectorCorrelation::ComplexImageType;
using CST = VectorCorrelation::ComplexScalarType;
TEST_CASE("image mean fails with no image", "[blocked image,filter]")
{
  VectorCorrelation::BlockedImageMean<Serial,ComplexImageType> blocked_image_mean{};
  REQUIRE_THROWS(blocked_image_mean.Run());
  ComplexImageType image{"my image", 10,10};
  blocked_image_mean.SetInput(image);
  REQUIRE_NOTHROW(blocked_image_mean.Run());
  // blocked image tasks only work with odd sized blocks
  REQUIRE_THROWS(blocked_image_mean.SetBlockSize(4));
  REQUIRE_NOTHROW(blocked_image_mean.SetBlockSize(3));

  // blocked image cannot run when block size > image size
  blocked_image_mean.SetBlockSize(11);
  REQUIRE_THROWS(blocked_image_mean.Run());
}

static void initialize_image(ComplexImageType image)
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

TEST_CASE("blocked image mean returns correct values", "[blocked image, filter]")
{
  VectorCorrelation::BlockedImageMean<Serial,ComplexImageType> blocked_image_mean{};
  ComplexImageType image{"my image", 10,10};
  initialize_image(image);
  blocked_image_mean.SetInput(image);
  // block result with block size 1 is the same as the input
  blocked_image_mean.SetBlockSize(1);
  blocked_image_mean.Run();
  auto mean = blocked_image_mean.GetOutput();
  REQUIRE(mean.extent(0) == image.extent(0));
  REQUIRE(mean.extent(1) == image.extent(1));

  int result = 0;
  for (int i=0; i<mean.extent(0); ++i)
  {
    for(int j=0; j<mean.extent(1); ++j)
    {
      double val = i*image.extent(0)+j;
      result += (mean(i,j) != CST{val,val});
    }
  }
  REQUIRE(result == 0);

  // result with block size same as image size is 
  // the mean in the central pixel with other pixels
  // zero
  ComplexImageType image2{"my image", 11,11};
  initialize_image(image2);
  blocked_image_mean.SetInput(image2);
  blocked_image_mean.SetBlockSize(11);
  blocked_image_mean.Run();
  mean = blocked_image_mean.GetOutput();
  REQUIRE(mean.extent(0) == image2.extent(0));
  REQUIRE(mean.extent(1) == image2.extent(1));
  REQUIRE(mean(5,5).real() == Catch::Approx(60));
  REQUIRE(mean(5,5).imag() == Catch::Approx(60));

  // compute mean with a block size of 3
  blocked_image_mean.SetBlockSize(3);
  blocked_image_mean.Run();
  mean = blocked_image_mean.GetOutput();
  result = 0;
  for (int i=0; i<mean.extent(0); ++i)
  {
    for(int j=0; j<mean.extent(1); ++j)
    {
      if(i==0 || j==0 || i==mean.extent(0)-1 || j == mean.extent(1)-1)
      {
        // the boundary of the image outside of the mask is supposed
        // to be zero
        result += (mean(i,j) != 0);
      }
      else
      {
        CST local_mean{0,0};
        for(int l=0; l<blocked_image_mean.GetBlockSize(); ++l)
        {
          for(int m=0; m<blocked_image_mean.GetBlockSize(); ++m)
          {
            int block_row = i-blocked_image_mean.GetBlockSize()/2+l;
            int block_col = j-blocked_image_mean.GetBlockSize()/2+m;
            // within a column each row is idential
            local_mean += image2(block_row, block_col);
          }
        }
        local_mean/= blocked_image_mean.GetBlockSize()*blocked_image_mean.GetBlockSize();
        result+= (mean(i,j).real() != Catch::Approx(local_mean.real()));
        result+= (mean(i,j).imag() != Catch::Approx(local_mean.imag()));
      }
    }
  }
}
