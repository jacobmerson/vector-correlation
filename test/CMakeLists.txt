
find_package(Catch2 3 REQUIRED)

add_library(vector_correlation_main vector_correlation_main.cc)
target_link_libraries(vector_correlation_main PUBLIC Catch2::Catch2 vector_correlation)

add_executable(test_correlation 
               test_blocked_image_auto_correlation.cc
               test_blocked_image_correlation.cc
               test_blocked_image_mean.cc
               test_complex_math.cc
               test_image_magnitude.cc
               test_vector_correlation.cc)
#target_link_libraries(test_correlation Catch2::Catch2WithMain)
target_link_libraries(test_correlation vector_correlation_main)

