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

TEST_CASE("Reads an integer","[data.cc:getint()]") {
	unsigned int i = getint("42", "anon", UNSET);
	REQUIRE(i == 42);
}

/*
TEST_CASE("DATA::getint NA") {
	DATA d;
	unsigned int i = d.getint("NA", "anon");
	CHECK(i == UNKNOWN);
}

TEST_CASE("DATA::getint *") {
	DATA d;
	unsigned int i = d.getint("*", "anon");
	CHECK(i == THRESH);
}
*/
