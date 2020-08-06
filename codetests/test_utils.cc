#include "../catch.hpp"

#include <iostream>
#include <iomanip>

#include "../utils.hh"

#define private public // Comments welcome :)
//#include "../data.hh"

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

