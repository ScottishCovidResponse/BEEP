#include "../catch.hpp"

#include <iostream>
#include <iomanip>

#include "../utils.hh"

#include "../consts.hh"

using namespace std;

/*
TEST_CASE("test_ran") {

	sran(0);
	double r = ran();

//	cout << "r == " << std::setprecision(19) << r << endl;

  CHECK(r == 0.5928446159645536229);
}

TEST_CASE("split") {
	vector<string> parts = split("ab*defgh**p*", '*');

	CHECK(parts.at(0) == "ab");
	CHECK(parts.at(1) == "defgh");
	CHECK(parts.at(2) == "");
	CHECK(parts.at(3) == "p");
	CHECK(parts.size() == 4);

	// ianhinder: I would have expected the following to work, but the last part
	// is not returned
	// CHECK(parts.at(4) == "");
	// CHECK(parts.size() == 5);
}
*/

///////////////////////////////////////////////////////////////////////////
// getint
///////////////////////////////////////////////////////////////////////////

const char* suite = "[data.cc:getint()]";
TEST_CASE("getint reads an integer",suite) {
	unsigned int i = getint("42", UNSET);
	REQUIRE(i == 42);
}

TEST_CASE("getint reads 'NA' as UNKNOWN",suite) {
	unsigned int i = getint("NA", UNSET);
	REQUIRE(i == UNKNOWN);
}

TEST_CASE("getint reads '*' as THRESH if threshold is given",suite) {
	unsigned int i = getint("*", 5);
	REQUIRE(i == THRESH);
}

TEST_CASE("getint fails on '*' if no threshold is given",suite) {
	unsigned int i;
	REQUIRE_THROWS_AS(i = getint("*", UNSET),std::runtime_error);
}

TEST_CASE("getint fails on junk",suite) {
	unsigned int i;
	REQUIRE_THROWS_AS(i = getint("xx", UNSET),std::runtime_error);
}

TEST_CASE("getint fails on trailing junk",suite) {
	unsigned int i;
	REQUIRE_THROWS_AS(i = getint("xx", UNSET),std::runtime_error);
}
