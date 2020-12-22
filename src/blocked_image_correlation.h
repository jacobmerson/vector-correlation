#ifndef BLOCKED_IMAGE_CORRELATION_H__
#define BLOCKED_IMAGE_CORRELATION_H__
#include "common_types.h"
namespace VectorCorrelation
{
  template <typename Backend, typename ImageType>
  class BlockedImageCorrelation
  {
    public:
      BlockedImageCorrelation() : block_size_(1) {};
      void Run();
      auto GetOutput() const {return correlation_image_;}
      void SetInput(ImageType image1, ImageType mean1, ScalarImageType sigma1, ImageType image2, ImageType mean2, ScalarImageType sigma2);
      void SetBlockSize(OrdinalType bs);
      const OrdinalType GetBlockSize() const {return block_size_;}
    private:
      ImageType image1_;
      ImageType image2_;
      ImageType mean_image1_;
      ImageType mean_image2_;
      ScalarImageType sigma_image1_;
      ScalarImageType sigma_image2_;
      ImageType correlation_image_;
      OrdinalType block_size_;
  };
}
#endif
