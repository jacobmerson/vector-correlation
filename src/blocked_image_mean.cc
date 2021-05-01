#include "blocked_image_mean.h"
#include "common_types.h"
#include <exception>
#include "fmt/format.h"
#include "Kokkos_Core.hpp"

namespace VectorCorrelation
{
  template<typename Backend, typename ImageType>
  void BlockedImageMean<Backend, ImageType>::Run()
  {

    if(image_.extent(0) <1 || image_.extent(1) <1)
    {
      throw std::runtime_error("BlockedImageMean filter must have an image assigned.");
    }
    if(block_size_ > image_.extent(0) || block_size_ > image_.extent(1))
    {
      throw std::runtime_error("Block Size must be smaller than image size.");
    }
    mean_image_ = ImageType{"mean", image_.extent(0), image_.extent(1)};
    OrdinalType buffer = block_size_ / 2;
    OrdinalType block_size_sq = block_size_*block_size_;
    Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy(
        {buffer, buffer, 0, 0}, {static_cast<long>(image_.extent(0) - buffer),
                                 static_cast<long>(image_.extent(1) - buffer), block_size_, block_size_});
    Kokkos::parallel_for("compute mean", policy, KOKKOS_LAMBDA(const int i, const int j, const int l, const int m)
        {
        auto block_row = i-buffer+l;
        auto block_col = j-buffer+m;
        // FIXME: replace with a scatter view
        Kokkos::atomic_add(&mean_image_(i,j),image_(block_row,block_col)/block_size_sq);
        });
  }
  template<typename Backend, typename ImageType>
  void BlockedImageMean<Backend,ImageType>::SetBlockSize(OrdinalType bs)
  {
    if(bs%2 == 0)
    {
      throw std::runtime_error("Block size must be an odd number");
    }
    block_size_ = bs;
  }
  template class BlockedImageMean<Serial, ComplexImageType>;
  template class BlockedImageMean<Serial, ScalarImageType>;
}
