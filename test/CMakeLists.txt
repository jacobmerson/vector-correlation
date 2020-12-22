
find_package(Catch2 3 REQUIRED)

add_library(vector_correlation_main vector_correlation_main.cc)
target_link_libraries(vector_correlation_main PUBLIC Catch2::Catch2 vector_correlation)

add_executable(test_correlation test_complex_math.cc test_vector_correlation.cc)
#target_link_libraries(test_correlation Catch2::Catch2WithMain)
target_link_libraries(test_correlation vector_correlation_main)

