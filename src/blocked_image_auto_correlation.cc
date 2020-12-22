#include "blocked_image_auto_correlation.h"
#include "complex_math.h"

namespace VectorCorrelation
{
  template <typename Backend, typename ImageType>
  void BlockedImageAutoCorrelation<Backend,ImageType>::Run()
  {
    if(image_.extent(0) <1 || image_.extent(1) <1)
    {
      throw std::runtime_error("BlockedImageAutoCorrelation filter must have an image assigned.");
    }
    if(block_size_ > image_.extent(0) || block_size_ > image_.extent(1))
    {
      throw std::runtime_error("Block Size must be smaller than image size.");
    }
    correlation_image_ = ScalarImageType{"auto correlation", image_.extent(0), image_.extent(1)};

    OrdinalType buffer = block_size_ / 2;
    OrdinalType block_size_sq = block_size_*block_size_;
    Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy(
        {buffer, buffer, 0, 0}, {static_cast<long>(image_.extent(0) - buffer),
                                 static_cast<long>(image_.extent(1) - buffer), block_size_, block_size_});
    Kokkos::parallel_for("compute autocorrelation", policy, KOKKOS_LAMBDA(const int i, const int j, const int l, const int m)
    {
        auto block_row = i-buffer+l;
        auto block_col = j-buffer+m;
        auto val = image_(block_row,block_col)-mean_image_(i,j);
        auto sigma = conjugate_square(val);
        // there should be no imaginary part of sigma
        Kokkos::atomic_add(&correlation_image_(i,j),sigma/block_size_sq);
    });
  }
  template<typename Backend, typename ImageType>
  void BlockedImageAutoCorrelation<Backend,ImageType>::SetBlockSize(OrdinalType bs)
  {
    if(bs%2 == 0)
    {
      throw std::runtime_error("Block size must be an odd number");
    }
    block_size_ = bs;
  }
  template<typename Backend, typename ImageType>
  void BlockedImageAutoCorrelation<Backend,ImageType>::SetInput(ImageType image, ImageType mean) {
      if((image.extent(0) != mean.extent(0)) || (image.extent(1) != mean.extent(1)))
      {
        throw std::runtime_error("Image and MeanImage must be the same size.");
      }
        image_ = image;
        mean_image_ = mean;
      }
  template class BlockedImageAutoCorrelation<Serial,ComplexImageType>;
}
