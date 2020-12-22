#ifndef BLOCKED_IMAGE_AUTO_CORRELATION_H__
#define BLOCKED_IMAGE_AUTO_CORRELATION_H__
#include "common_types.h"

namespace VectorCorrelation
{
  template <typename Backend, typename ImageType>
  class BlockedImageAutoCorrelation
  {
    public:
      BlockedImageAutoCorrelation() : block_size_(1) {};
      void Run();
      ScalarImageType GetOutput() const {return correlation_image_;}
      void SetInput(ImageType image, ImageType mean);
      void SetBlockSize(OrdinalType bs);
      const OrdinalType GetBlockSize() const {return block_size_;}
    private:
      ImageType image_;
      ImageType mean_image_;
      ScalarImageType correlation_image_;
      OrdinalType block_size_;
  };
}
#endif
