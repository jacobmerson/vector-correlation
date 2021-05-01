#include "vector_correlation_analysis.h"
#include "image_magnitude.h"
#include "blocked_image_mean.h"
#include "blocked_image_correlation.h"
#include "blocked_image_auto_correlation.h"

namespace VectorCorrelation
{
  template <typename Backend>
  void VectorCorrelationAnalysis<Backend>::AddFrame(ComplexImageType orientation_vector,
                                            OrdinalImageType mask)
  {
    spdlog::debug("Adding frame by orienation vector and mask");
    // FIXME if frame sizes are not correct...throw
    if (orientation_vector.extent(0) != NX_ ||
        orientation_vector.extent(1) != NY_ || mask.extent(0) != NX_ ||
        mask.extent(1) != NY_)
    {
      throw std::length_error(fmt::format(
          "orientation_vector({},{}) or mask({},{}) does not have correct "
          "size({},{})",
          orientation_vector.extent(0), orientation_vector.extent(1),
          mask.extent(0), mask.extent(1), NX_, NY_));
    }
    orientation_vector_frames_.push_back(orientation_vector);
    mask_frames_.push_back(mask);
  }
  template <typename Backend>
  void VectorCorrelationAnalysis<Backend>::AddFrame(ComplexImageType orientation_vector)
  {
    OrdinalImageType mask("mask", orientation_vector.extent(0), orientation_vector.extent(1));
    Kokkos::deep_copy(mask,1);
    AddFrame(orientation_vector,mask);
  }
  template <typename Backend>
  void VectorCorrelationAnalysis<Backend>::Run()
  {
    spdlog::debug("Run Correlation Analysis");
    vector_correlation_frames_.reserve(orientation_vector_frames_.size() - 1);
    // initialize filters we will use
    ImageMagnitude<Serial, ComplexImageType> image_magnitude;
    BlockedImageMean<Serial, ComplexImageType> blocked_image_mean;
    BlockedImageAutoCorrelation<Serial, ComplexImageType> blocked_image_auto_correlation;
    BlockedImageCorrelation<Serial, ComplexImageType> blocked_image_correlation;
    blocked_image_mean.SetBlockSize(number_correlation_pixels_);
    blocked_image_auto_correlation.SetBlockSize(number_correlation_pixels_);
    blocked_image_correlation.SetBlockSize(number_correlation_pixels_);

    auto image1 = orientation_vector_frames_[0];
    blocked_image_mean.SetInput(image1);
    blocked_image_mean.Run();
    //auto mask1 = mask_frames_[0];
    auto mean1 = blocked_image_mean.GetOutput();
    blocked_image_auto_correlation.SetInput(image1, mean1);
    blocked_image_auto_correlation.Run();
    auto sigma1 = blocked_image_auto_correlation.GetOutput();
    for (int i = 1; i < orientation_vector_frames_.size(); ++i)
    {
      auto image2 = orientation_vector_frames_[i];
      //auto mask2 = mask_frames_[i];
      blocked_image_mean.SetInput(image2);
      blocked_image_mean.Run();
      //auto mask1 = mask_frames_[0];
      auto mean2 = blocked_image_mean.GetOutput();
      blocked_image_auto_correlation.SetInput(image2, mean2);
      blocked_image_auto_correlation.Run();
      auto sigma2 = blocked_image_auto_correlation.GetOutput();
      blocked_image_correlation.SetInput(image1, mean1, sigma1, image2, mean2, sigma2);
      blocked_image_correlation.Run();
      auto vector_correlation = blocked_image_correlation.GetOutput();
      image_magnitude.SetInput(vector_correlation);
      image_magnitude.Run();
      vector_correlation_frames_.push_back(image_magnitude.GetOutput());
      // swap image 2 to image 1 so we don't need to recompute the values
      image1 = image2;
      //mask1 = mask2;
      mean1 = mean2;
      sigma1 = sigma2;
    }
  }
  // explicitly instantiate backend types
  template class VectorCorrelationAnalysis<Serial>;
}  // namespace VectorCorrelation
