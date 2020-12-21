#include "vector_correlation.h"
#include <exception>
namespace VectorCorrelation
{
  KOKKOS_INLINE_FUNCTION
  auto conjugate_times(CST a, CST b) -> CST
  {
    return CST{a.real() * b.real() + a.imag() * b.imag(),
               a.real() * b.imag() - b.real() * a.imag()};
  }

  template <typename Backend>
  auto compute_mean_image(ComplexImageType image,
                        OrdinalImageType mask,
                        OrdinalType ncorrpx)  // ->ScalarImageType
  {
    spdlog::debug("Computing mean image");
    spdlog::critical("Currently ignores mask!");
    ComplexImageType mean("mean", image.extent(0),
                                image.extent(1));
    OrdinalType buffer = ncorrpx / 2;
    OrdinalType ncorrpxsq = ncorrpx*ncorrpx;
    Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy(
        {buffer, buffer, 0, 0}, {image.extent(0) - buffer,
                                 image.extent(1) - buffer, ncorrpx, ncorrpx});
    Kokkos::parallel_for("compute mean", policy, KOKKOS_LAMBDA(const int i, const int j, const int l, const int m)
        {
        // FIXME: replace with a scatter view
        Kokkos::atomic_add(&mean(i,j),image(i-buffer+l,j-buffer+m)/ncorrpxsq);
        });
    return mean;
  }
  
  template<typename Backend>
  auto complex_image_magnitude(ComplexImageType image)
  {
    ScalarImageType magnitude("magnitude", image.extent(0), image.extent(1));
    Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy({0,0}, {image.extent(0), image.extent(1)});
    Kokkos::parallel_for("compute magnitude", policy, KOKKOS_LAMBDA(const int i, const int j)
    {
        magnitude(i,j) = sqrt(conjugate_times(image(i,j), image(i,j)).real());
     });
    return magnitude;
  }
                                
  /*
   * computes the magnitude of vector correlation between image1 and image2
   * \returns scalar image with the magnitude of the vector correlation
   */
  template <typename Backend>
  auto correlate_images(ComplexImageType image1,
                        OrdinalImageType mask1,
                        ComplexImageType mean1,
                        ScalarImageType sigma1,
                        ComplexImageType image2,
                        OrdinalImageType mask2,
                        ComplexImageType mean2,
                        ScalarImageType sigma2,
                        OrdinalType ncorrpx)  // ->ScalarImageType
  {
    spdlog::debug("Computing vector correlation");
    spdlog::critical("Current implementation ignores image masks");
    ComplexImageType correlation("correlation", image1.extent(0),
                                image1.extent(1));
    // loop over  every pixel and every subgrid within the pixels
    OrdinalType buffer = ncorrpx / 2;
    OrdinalType ncorrpxsq = ncorrpx*ncorrpx;
    Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy(
        {buffer, buffer, 0, 0}, {image1.extent(0) - buffer,
                                 image1.extent(1) - buffer, ncorrpx, ncorrpx});
    Kokkos::parallel_for("compute correlation", policy, KOKKOS_LAMBDA(const int i, const int j, const int l, const int m)
    {
        int subimage_i = i-buffer+l;
        int subimage_j = j-buffer+m;
        CST sigma = conjugate_times(image1(subimage_i,subimage_j)-mean1(i,j), image2(subimage_i,subimage_j)-mean2(i,j));
        sigma /= (ncorrpx*sigma1(i,j)*sigma2(i,j));
        Kokkos::atomic_add(&correlation(i,j),sigma/ncorrpxsq);
    });
    auto correlation_magnitude = complex_image_magnitude<Backend>(correlation);
    return correlation_magnitude;
  }
  /*
   * compute the auto correlation between two images
   */
  template <typename Backend>
  auto compute_auto_correlation(ComplexImageType image,
                                OrdinalImageType mask,
                                ComplexImageType mean,
                                OrdinalType ncorrpx)
  {
    spdlog::debug("Computing auto correlation");
    spdlog::critical("Current implementation ignores image masks");
    ScalarImageType correlation("auto correlation", image.extent(0),
                                image.extent(1));
    // loop over  every pixel and every subgrid within the pixels
    OrdinalType buffer = ncorrpx / 2;
    OrdinalType ncorrpxsq = ncorrpx*ncorrpx;
    Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy(
        {buffer, buffer, 0, 0}, {image.extent(0) - buffer,
                                 image.extent(1) - buffer, ncorrpx, ncorrpx});
    Kokkos::parallel_for("compute autocorrelation", policy, KOKKOS_LAMBDA(const int i, const int j, const int l, const int m)
    {
        int subimage_i = i-buffer+l;
        int subimage_j = j-buffer+m;
        CST sigma = conjugate_times(image(subimage_i,subimage_j)-mean(i,j), image(subimage_i,subimage_j)-mean(i,j));
        // there should be no imaginary part of sigma
        Kokkos::atomic_add(&correlation(i,j),sigma.real()/ncorrpxsq);
    });
    return correlation;
  }
  template <typename Backend>
  void VectorCorrelation<Backend>::AddFrame(ComplexImageType orientation_vector,
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
  void VectorCorrelation<Backend>::Run()
  {
    spdlog::debug("Run Correlation Analysis");
    vector_correlation_frames_.resize(orientation_vector_frames_.size() - 1);
    auto image1 = orientation_vector_frames_[0];
    auto mask1 = mask_frames_[0];
    auto mean1 = compute_mean_image<Backend>(image1,mask1,number_correlation_pixels_);
    auto sigma1 = compute_auto_correlation<Backend>(image1, mask1, mean1, number_correlation_pixels_);
    for (int i = 1; i < orientation_vector_frames_.size(); ++i)
    {
      auto image2 = orientation_vector_frames_[i];
      auto mask2 = mask_frames_[i];
      auto mean2 = compute_mean_image<Backend>(image2,mask2,number_correlation_pixels_);
      auto sigma2 = compute_auto_correlation<Backend>(image2, mask2, mean2, number_correlation_pixels_);
      auto correlated = correlate_images<Backend>(
          image1, mask1, mean1, sigma1,
          image2, mask1, mean2, sigma2,
          number_correlation_pixels_);
      // swap image 2 to image 1 so we don't need to recompute the values
      image1 = image2;
      mask1 = mask2;
      mean1 = mean2;
      sigma1 = sigma2;
    }
  }
}  // namespace VectorCorrelation
