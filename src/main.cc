#include <cstdlib>  // for EXIT_FAILURE
#include "common_types.h"
#include "spdlog/spdlog.h"
#include "vector_correlation_analysis.h"
#include "fmt/format.h"
using CST = VectorCorrelation::ComplexScalarType;
int main(int argc, char ** argv)
{
  spdlog::set_level(spdlog::level::debug);
  spdlog::info("Vector Correlation Code");
  VectorCorrelation::ScopeGuard scope_guard(argc, argv);
  constexpr int NX = 10;
  constexpr int NY = 10;
  VectorCorrelation::ComplexImageType image1("image1", NX, NY);
  VectorCorrelation::ComplexImageType image2("image2", NX, NY);
  VectorCorrelation::OrdinalImageType mask("mask", NX, NY);
  for (int i = 0; i < NX; ++i)
  {
    for (int j = 0; j < NY; ++j)
    {
      if(i<5)
        image1(i, j) = CST{1, j*0.5};
      else
        image1(i, j) = CST{i, j*0.5};
      image2(i, j) = CST{j*0.5,1};
      mask(i, j) = 1;
    }
  }
  VectorCorrelation::VectorCorrelationAnalysis<VectorCorrelation::Serial>
      vector_correlation{NX, NY};
  vector_correlation.AddFrame(image1, mask);
  vector_correlation.AddFrame(image2, mask);
  vector_correlation.Run();
  auto correlations = vector_correlation.GetOutput();
  for(auto& frame: correlations)
  {
    for(int i=0; i<frame.extent(0); ++i)
    {
      for(int j=0; j<frame.extent(1); ++j)
      {
        fmt::print("{:2.2f} ",frame(i,j));
      }
      fmt::print("\n");
    }
    fmt::print("\n");
  }
  return 0;
}
