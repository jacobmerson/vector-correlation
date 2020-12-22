#include "catch2/catch_session.hpp"
#include "vector_correlation.h"

int main(int argc, char*argv[]) 
{
  VectorCorrelation::ScopeGuard scope_guard(argc, argv);
  int result = Catch::Session().run(argc, argv);
  return result;
}
