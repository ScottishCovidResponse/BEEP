#include "../catch.hpp"

//#include <iostream>
//#include <iomanip>
//#include <string>

#include "../data.hh"

//#include "../consts.hh"

using namespace std;

///////////////////////////////////////////////////////////////////////////
// getint
///////////////////////////////////////////////////////////////////////////

const char* suite = "[data.cc:getint()]";
TEST_CASE("getint reads an integer",suite) {
	unsigned int i = getint("42", "anon", UNSET);
	REQUIRE(i == 42);
}

TEST_CASE("getint reads 'NA' as UNKNOWN",suite) {
	unsigned int i = getint("NA", "anon", UNSET);
	REQUIRE(i == UNKNOWN);
}

TEST_CASE("getint reads '*' as THRESH if threshold is given",suite) {
	unsigned int i = getint("*", "anon", 5);
	REQUIRE(i == THRESH);
}

TEST_CASE("getint fails on junk",suite) {
	unsigned int i;
	REQUIRE_THROWS_AS(i = getint("xx", "anon", UNSET),std::runtime_error);
}

TEST_CASE("getint fails on trailing junk",suite) {
	unsigned int i;
	REQUIRE_THROWS_AS(i = getint("xx", "anon", UNSET),std::runtime_error);
}
