#include "catch2/catch_test_macros.hpp"
#include "common_types.h"
#include "complex_math.h"

TEST_CASE("conjugate_times", "[complex]")
{
  CST num1{7,2};
  CST num2{1,10};
  auto result = VectorCorrelation::conjugate_times(num1, num1);
  REQUIRE(result.real() == num1.real()*num1.real()+num1.imag()*num1.imag());
  REQUIRE(result.imag() == 0);

  result = VectorCorrelation::conjugate_times(num1, num2);
  REQUIRE(result.real() == 27);
  REQUIRE(result.imag() == 68);

}
