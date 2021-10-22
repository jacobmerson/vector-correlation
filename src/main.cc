#include <cstdlib>  // for EXIT_FAILURE
#include "common_types.h"
#include "spdlog/spdlog.h"
#include "vector_correlation_analysis.h"
#include "fmt/format.h"
#include <fstream>
#include <filesystem>
#include "blocked_image_mean.h"

using CST = VectorCorrelation::ComplexScalarType;
using VectorCorrelation::ComplexImageType;
using VectorCorrelation::OrdinalImageType;
using VectorCorrelation::ScalarImageType;
using VectorCorrelation::VectorCorrelationAnalysis;
using VectorCorrelation::Serial;
using VectorCorrelation::BlockedImageMean;
namespace fs = std::filesystem;

struct ImageFrame {
  ImageFrame(fs::path a,fs::path r) : angle(std::move(a)), retardance(std::move(r)) {
    if (!fs::exists(angle)) {
      spdlog::error("Angle Path: {} doesn't exist!\n", angle.c_str());
      std::abort();
    }
    if (!fs::exists(retardance)) {
      spdlog::error("Retardance Path: {} doesn't exist!\n", retardance.c_str());
      std::abort();
    }
  }
  fs::path angle;
  fs::path retardance;
};

// Read image data from fstream
std::vector<std::vector<double>> read_data(std::ifstream& ifile);
ComplexImageType read_complex_image(ImageFrame frame);

int main(int argc, char ** argv)
{
  spdlog::set_level(spdlog::level::warn);
  spdlog::info("Vector Correlation Code");
  VectorCorrelation::ScopeGuard scope_guard(argc, argv);
  fs::path base_path{"/lore/mersoj/vector-correlation/AR_alignments"};
  auto frame1 = read_complex_image(ImageFrame{base_path / "align9_0_angle.dat",base_path / "align9_0_alignment.dat"});
  VectorCorrelationAnalysis<Serial> vector_correlation(frame1.extent(0), frame1.extent(1));
  vector_correlation.AddFrame(frame1);
  //vector_correlation.AddFrame(read_complex_image(ImageFrame{base_path / fmt::format("align9_{}_angle.dat",1),base_path / fmt::format("align9_{}_alignment.dat",1)}));
  for(int i=10;i<50; i+=10)
  {
    vector_correlation.AddFrame(read_complex_image(ImageFrame{base_path / fmt::format("align9_{}_angle.dat",i),base_path / fmt::format("align9_{}_alignment.dat",i)}));
  }
  /*
  //fs::path base_path{"/lore/mersoj/vector-correlation/R54-csv"};
  fs::path base_path{"/lore/mersoj/vector-correlation/R55"};

  auto frame1 = read_complex_image(ImageFrame{base_path / "min100_Ang.csv",base_path / "min100_Ret.csv"});
  VectorCorrelationAnalysis<Serial> vector_correlation(frame1.extent(0), frame1.extent(1));
  vector_correlation.AddFrame(frame1);
  vector_correlation.AddFrame(read_complex_image(ImageFrame{base_path / "min80_Ang.csv",base_path / "min80_Ret.csv"}));
  vector_correlation.AddFrame(read_complex_image(ImageFrame{base_path / "min60_Ang.csv",base_path / "min60_Ret.csv"}));
  vector_correlation.AddFrame(read_complex_image(ImageFrame{base_path / "min40_Ang.csv",base_path / "min40_Ret.csv"}));
  vector_correlation.AddFrame(read_complex_image(ImageFrame{base_path / "min20_Ang.csv",base_path / "min20_Ret.csv"}));
  vector_correlation.AddFrame(read_complex_image(ImageFrame{base_path / "AR_Ang.csv",base_path / "AR_Ret.csv"}));
  */

  vector_correlation.Run();
  auto correlations = vector_correlation.GetOutput();
  BlockedImageMean<Serial, ScalarImageType> mean;
  mean.SetBlockSize(3);
  for(auto& frame: correlations)
  {
    mean.SetInput(frame);
    mean.Run();
    auto mean_image = mean.GetOutput();
    fmt::print("# Frame\n");
    for(int i=0; i<mean_image.extent(0); ++i)
    {
      for(int j=0; j<mean_image.extent(1); ++j)
      {
        fmt::print("{:2.2f} ",mean_image(i,j));
      }
      fmt::print("\n");
    }
  }

  return 0;
}
std::vector<std::vector<double>> read_data(std::ifstream& ifile)
{
  std::vector<std::vector<double>> data;
    for(std::string line; std::getline(ifile, line);)
    {
      data.emplace_back();
      size_t end=0;
      size_t start = 0;
      do {
        //end = line.find(',',start);
        end = line.find(',',start);
        double val = std::stod(line.substr(start,end));
        data.back().push_back(val);
        if(end == std::string::npos)
        {
          break;
        }
        start = end+1;
      }while(true);
      if(data[0].size() != data.back().size())
      {
        spdlog::error("all rows of input data must have the same dimension");
        exit(1);
      }
    }

  spdlog::debug("NX: {}, NY: {}", data[0].size(), data.size());
  return data;
}
ComplexImageType read_complex_image(ImageFrame frame) {
  std::ifstream iangle(frame.angle);
  std::ifstream iretardance(frame.retardance);
  auto angle = read_data(iangle);
  auto retardance = read_data(iretardance);
  if(angle.size() != retardance.size() ||
      angle[0].size() != retardance[0].size()) {
    spdlog::error("angle and retardance must be the same size!");
    std::abort();
  }

  auto nrows = angle.size();
  auto ncols = angle[0].size();
  ComplexImageType image("image", nrows, ncols);
  for(size_t r=0; r<nrows; ++r) {
    for(size_t c=0; c<ncols; ++c) {
      //image(r,c) = CST{std::cos(angle[r][c])*retardance[r][c], std::sin(angle[r][c])*retardance[r][c]};
      image(r,c) = CST{std::cos(angle[r][c]), std::sin(angle[r][c])};
      //image(r,c) = CST{std::cos(2*angle[r][c])*std::pow(std::sin(retardance[r][c]),2), std::sin(2*angle[r][c])*std::pow(std::sin(retardance[r][c]),2)};
    }
  }
  return image;
}
