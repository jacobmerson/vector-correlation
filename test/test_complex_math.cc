#include "catch2/catch_test_macros.hpp"
#include "catch2/catch_approx.hpp"
#include "common_types.h"
#include "complex_math.h"

using Catch::Approx;
using VectorCorrelation::conjugate_times;
using VectorCorrelation::conjugate_square;
using VectorCorrelation::magnitude;

using CST = VectorCorrelation::ComplexScalarType;
TEST_CASE("conjugate_times", "[complex]")
{
  CST num1{7,2};
  CST num2{1,10};
  auto result = conjugate_times(num1, num1);
  REQUIRE(result.real() == num1.real()*num1.real()+num1.imag()*num1.imag());
  REQUIRE(result.imag() == Approx(0));

  result = conjugate_times(num2, num2);
  REQUIRE(result.real() == num2.real()*num2.real()+num2.imag()*num2.imag());
  REQUIRE(result.imag() == Approx(0));

  result = conjugate_times(num1, num2);
  REQUIRE(result.real() == Approx(27));
  REQUIRE(result.imag() == Approx(68));

  double real1 = 3.0;
  double real2 = 7.0;
  REQUIRE(conjugate_times(real1,real2) == Approx(21));
}
TEST_CASE("conjugate_square", "[complex]")
{
  CST complex1{7,2};
  CST complex2{1,10};
  double real1 = 3.0;
  double real2 = 7.0;
  // test complex values
  REQUIRE(conjugate_square(complex1) == Approx(conjugate_times(complex1,complex1).real()));
  REQUIRE(conjugate_square(complex2) == Approx(conjugate_times(complex2,complex2).real()));
  // test real values
  REQUIRE(conjugate_square(real1) == Approx(real1*real1));
  REQUIRE(conjugate_square(real2) == Approx(real2*real2));
}
TEST_CASE("magnitude", "[complex]")
{
  CST complex1{7,2};
  CST complex2{1,10};
  double real1 = 3.0;
  double real2 = 7.0;
  // test complex values
  REQUIRE(magnitude(complex1) == Approx(sqrt(conjugate_times(complex1,complex1).real())));
  REQUIRE(magnitude(complex2) == Approx(sqrt(conjugate_times(complex2,complex2).real())));
  // test real values
  REQUIRE(magnitude(real1) == Approx(sqrt(real1*real1)));
  REQUIRE(magnitude(real2) == Approx(sqrt(real2*real2)));
}
