add_library(vector_correlation 
            blocked_image_auto_correlation.cc
            blocked_image_correlation.cc
            blocked_image_mean.cc
            complex_math.cc 
            image_magnitude.cc
            vector_correlation_analysis.cc)
target_link_libraries( vector_correlation PUBLIC Kokkos::kokkos fmt::fmt)
target_compile_definitions(vector_correlation PUBLIC SPDLOG_FMT_EXTERNAL)
target_include_directories(vector_correlation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(main main.cc)
target_link_libraries(main PUBLIC spdlog vector_correlation)
