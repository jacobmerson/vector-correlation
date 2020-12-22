#include "blocked_image_correlation.h"
#include "complex_math.h"

namespace VectorCorrelation
{
  template <typename Backend, typename ImageType>
  void BlockedImageCorrelation<Backend, ImageType>::Run()
  {
    if(image1_.extent(0) <1 || image1_.extent(1) <1)
    {
      throw std::runtime_error("BlockedImageCorrelation filter must have an image assigned.");
    }
    if(block_size_ > image1_.extent(0) || block_size_ > image1_.extent(1))
    {
      throw std::runtime_error("Block Size must be smaller than image size.");
    }
    correlation_image_ = ImageType{"correlation", image1_.extent(0), image1_.extent(1)};
    // loop over  every pixel and every subgrid within the pixels
    OrdinalType buffer = block_size_ / 2;
    OrdinalType block_size_sq = block_size_*block_size_;
    Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy(
        {buffer, buffer, 0, 0}, {static_cast<long>(image1_.extent(0) - buffer),
                                 static_cast<long>(image1_.extent(1) - buffer), block_size_, block_size_});
    Kokkos::parallel_for("compute correlation", policy, KOKKOS_LAMBDA(const int i, const int j, const int l, const int m)
    {
        auto block_row = i-buffer+l;
        auto block_col = j-buffer+m;
        auto sigma = conjugate_times(image1_(block_row,block_col)-mean_image1_(i,j), image2_(block_row,block_col)-mean_image2_(i,j));
        sigma /= (block_size_sq*sqrt(sigma_image1_(i,j))*sqrt(sigma_image2_(i,j)));
        Kokkos::atomic_add(&correlation_image_(i,j),sigma);
    });
  }
  template<typename Backend, typename ImageType>
  void BlockedImageCorrelation<Backend,ImageType>::SetBlockSize(OrdinalType bs)
  {
    if(bs%2 == 0)
    {
      throw std::runtime_error("Block size must be an odd number");
    }
    block_size_ = bs;
  }
  template<typename Backend, typename ImageType>
  void BlockedImageCorrelation<Backend,ImageType>::SetInput(ImageType image1, ImageType mean1,
      ScalarImageType sigma1, ImageType image2, ImageType mean2, ScalarImageType sigma2) {
      if((image1.extent(0) != mean1.extent(0)) || (image1.extent(1) != mean1.extent(1)) ||
        (image1.extent(0) != sigma1.extent(0)) || (image1.extent(1) != sigma1.extent(1)) ||
        (image1.extent(0) != image2.extent(0)) || (image1.extent(1) != image2.extent(1)) ||
        (image2.extent(0) != mean2.extent(0)) || (image2.extent(1) != mean2.extent(1)) ||
        (image2.extent(0) != sigma2.extent(0)) || (image2.extent(1) != sigma2.extent(1)))
      {
        throw std::runtime_error("All input images must be the same size.");
      }
        image1_ = image1;
        mean_image1_ = mean1;
        sigma_image1_ = sigma1;
        image2_ = image2;
        mean_image2_ = mean2;
        sigma_image2_ = sigma2;
      }
  template class BlockedImageCorrelation<Serial, ComplexImageType>;
}
