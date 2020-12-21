
add_library(vector_correlation vector_correlation.cc)
target_link_libraries( vector_correlation PUBLIC Kokkos::kokkos fmt::fmt)
target_compile_definitions(vector_correlation PUBLIC SPDLOG_FMT_EXTERNAL)

add_executable(main main.cc)
target_link_libraries(main PUBLIC spdlog vector_correlation)
