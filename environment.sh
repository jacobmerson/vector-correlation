## CMAKE 3.19.1
#spack load gcc@6.5.0 # needed for gcc 10 symbols
#spack load gcc@10.1.0
#spack load /msm7qpr
#spack load spdlog
#spack load fmt
### Load OMP Kokkos
#spack load /3sbuwnf
##spack load cmake@3.19.1 spdlog fmt /3sbuwnf

spack load llvm@11.0.0
spack load /zj6wmta # cmake
spack load spdlog%clang@11.0.0
spack load fmt%clang@11.0.0
spack load kokkos%clang@11.0.0
spack load catch2%clang@11.0.0

