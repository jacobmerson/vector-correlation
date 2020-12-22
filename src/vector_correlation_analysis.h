#ifndef VECTOR_CORRELATION_H__
#define VECTOR_CORRELATION_H__
#include <istream>
#include <memory>
#include "Kokkos_Core.hpp"
#include "common_types.h"
#include "spdlog/spdlog.h"
namespace VectorCorrelation
{
  class ScopeGuard
  {
    public:
    // initialize any libraries that need initialization
    ScopeGuard(int & narg, char * arg[])
    {
      spdlog::debug("Initializing Kokkos");
      // initialize kokkos scope guard
      kokkos_guard_ = std::make_unique<Kokkos::ScopeGuard>(narg, arg);
    }

    private:
    std::unique_ptr<Kokkos::ScopeGuard> kokkos_guard_;
  };
  // main class that holds vector correlation data
  template <typename Backend>
  class VectorCorrelationAnalysis
  {
    public:
    using BackendType = Backend;
    /*
     *
     */
    VectorCorrelationAnalysis(OrdinalType NX, OrdinalType NY)
        : NX_(NX), NY_(NY), number_correlation_pixels_(5){};
    void Run();
    /**
     * Add frame from input streams
     */
    // void AddFrame(std::istream& orientation, std::istream& magnitude,
    // std::istream& mask, delimiter='\t');
    /**
     * Add frame from orientation and magnitude images
     */
    // void AddFrame(ScalarImageType orientation, ScalarImageType magnitude,
    // OrdinalImageType mask);
    /*
     * add frame directly from orientation vector and mask
     */
    void AddFrame(ComplexImageType orientation_vector, OrdinalImageType mask);
    const auto& GetOutput() const {return vector_correlation_frames_;}
    private:
    /**
     * vector frames over which the vector correlation will be computed
     * these vectors are a combination of the orientation and magnitude
     * of alignment data
     */
    std::vector<ComplexImageType> orientation_vector_frames_;
    /**
     * the masks for each vector correlation frame.
     * values greater than 0 correspond to data which will be included
     * in the correlation. Values less than or equal to 0 will be ignored
     * in the vector correlation.
     */
    std::vector<OrdinalImageType> mask_frames_;
    /**
     * The output vector correlatio magnitude image
     * This image will have a halo of number_correlation_pixels_/2
     * rows and columns with a zero value.
     */
    std::vector<ScalarImageType> vector_correlation_frames_;
    OrdinalType NX_;
    OrdinalType NY_;
    OrdinalType number_correlation_pixels_;
    // bool ignore_max_vector
    // bool ignore_vector_magnitude
  };
}  // namespace VectorCorrelation
#endif
