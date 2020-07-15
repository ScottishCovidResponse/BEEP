#include "../catch.hpp"

#include <iostream>
#include <iomanip>
#include <string>

#define private public
#include "../data.hh"

#include "../consts.hh"

using namespace std;

///////////////////////////////////////////////////////////////////////////
// getint
///////////////////////////////////////////////////////////////////////////

TEST_CASE("DATA::getint 45") {
	DATA d;
	unsigned int i = d.getint("45", "anon");
	CHECK(i == 45);
}

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
