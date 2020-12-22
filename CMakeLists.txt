cmake_minimum_required(VERSION 3.19.1)
project(VectorCorrelation CXX)

include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# This block of code is not needed when gcc/clang is properly loaded
# Not entirely sure why cmake isn't finding threading library
# https://github.com/alicevision/geogram/issues/2
#set(CMAKE_THREAD_LIBS_INIT "-lpthread")
#set(CMAKE_HAVE_THREADS_LIBRARY 1)
#set(CMAKE_USE_WIN32_THREADS_INIT 0)
#set(CMAKE_USE_PTHREADS_INIT 1)
#set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package(fmt REQUIRED)
find_package(spdlog REQUIRED)
find_package(Kokkos REQUIRED)

add_subdirectory(src)

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
