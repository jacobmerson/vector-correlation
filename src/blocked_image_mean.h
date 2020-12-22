#ifndef BLOCKED_IMAGE_MEAN_H__
#define BLOCKED_IMAGE_MEAN_H__
#include "common_types.h"
#include "Kokkos_Core.hpp"

namespace VectorCorrelation
{
/*
 * Compute the mean of a small block within the image and
 * set the mean of that block to the central pixel of the block
 * Backend is the Kokkos backend to use for the computation
 */
template <typename Backend, typename ImageType>
class BlockedImageMean
{
  static_assert(Kokkos::SpaceAccessibility<Backend, Kokkos::HostSpace>::accessible, "currently backend must be able to use CPU memory");
  public:
    BlockedImageMean() :mean_image_(ImageType{}), block_size_(1) {};
    void Run();
    ImageType GetOutput() const {return mean_image_;}
    void SetInput(ImageType image) {image_ = image;}
    void SetBlockSize(OrdinalType bs);
    const OrdinalType GetBlockSize() const {return block_size_;}

  private:
    ImageType image_;
    ImageType mean_image_;
    OrdinalType block_size_;
};
}
#endif
