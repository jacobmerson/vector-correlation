#include <cstdlib>  // for EXIT_FAILURE
#include <filesystem>
#include "common_types.h"
#include "spdlog/spdlog.h"
#include "vector_correlation_analysis.h"
int main(int argc, char ** argv)
{
  spdlog::set_level(spdlog::level::debug);
  spdlog::info("Vector Correlation Code");
  VectorCorrelation::ScopeGuard scope_guard(argc, argv);
  constexpr int NX = 1000;
  constexpr int NY = 1000;
  VectorCorrelation::ComplexImageType image1("image1", NX, NY);
  VectorCorrelation::ComplexImageType image2("image2", NX, NY);
  VectorCorrelation::OrdinalImageType mask("mask", NX, NY);
  for (int i = 0; i < NX; ++i)
  {
    for (int j = 0; j < NY; ++j)
    {
      image1(i, j) = CST{1, 1};
      image2(i, j) = CST{1, 1};
      mask(i, j) = 1;
    }
  }
  VectorCorrelation::VectorCorrelationAnalysis<VectorCorrelation::Serial>
      vector_correlation{NX, NY};
  vector_correlation.AddFrame(image1, mask);
  vector_correlation.AddFrame(image2, mask);
  vector_correlation.Run();
  return 0;
}
