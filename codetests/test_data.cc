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

TEST_CASE("getint reads UNKNOWN for NA",suite) {
	unsigned int i = getint("NA", "anon", UNSET);
	REQUIRE(i == UNKNOWN);
}

/*
TEST_CASE("DATA::getint *") {
	DATA d;
	unsigned int i = d.getint("*", "anon");
	REQUIRE(i == THRESH);
}
*/
