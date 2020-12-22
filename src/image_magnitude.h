#ifndef IMAGE_MAGNITUDE_H__
#define IMAGE_MAGNITUDE_H__
#include "common_types.h"
#include <type_traits>

namespace VectorCorrelation
{
  template <typename Backend, typename ImageType>
  class ImageMagnitude
  {
    public:
      void Run();
      void SetInput(ImageType image) {image_ = image;}
      ScalarImageType GetOutput() const {return magnitude_;}
    private:
      ImageType image_;
      ScalarImageType magnitude_;
  };
}
//#include "image_magnitude_impl.h"
#endif
