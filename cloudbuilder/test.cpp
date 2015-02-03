#define CATCH_CONFIG_MAIN
#include "../catch/single_include/catch.hpp"
#include "./build_clouds.hpp"

TEST_CASE( "Factorials are computed", "[factorial]" ) {
	REQUIRE( pclouds::Factorial(0) == 1 );
	REQUIRE( pclouds::Factorial(1) == 1 );
	REQUIRE( pclouds::Factorial(2) == 2 );
	REQUIRE( pclouds::Factorial(3) == 6 );
}
