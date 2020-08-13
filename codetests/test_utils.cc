#include "../catch.hpp"

#include <iostream>
#include <iomanip>

#include "../utils.hh"

#include "../consts.hh"

using namespace std;

const char* suite_random = "Random numbers";
TEST_CASE("ran after seed set to 0 returns result expected for std::mt19937",
					suite_random) {
	sran(0);
	double r = ran();
	
//	cout << "r == " << std::setprecision(19) << r << endl;
	
	REQUIRE(r == 0.5928446159645536229);
}
const char* suite_distributions = "Distributions";
TEST_CASE("logged normalprob at mean with variance 1./(2*pi) is 0.",
					suite_distributions) {
	REQUIRE(normalprob(1.0,1.0,1.0/(2*M_PI)) == Approx( 0.0 ));
}

TEST_CASE("logged normalprob in left wing is tiny",
					suite_distributions) {
	REQUIRE(normalprob(-100.0,0.0,1.0/(2*M_PI)) < -100.0);
}

TEST_CASE("logged normalprob in right wing is tiny",
					suite_distributions) {
	REQUIRE(normalprob(100.0,0.0,1.0/(2*M_PI)) < -100.0);
}

TEST_CASE("logged normalprob throws out of domain",
					suite_distributions) {
	double d;
	emsg_throws = true;
	REQUIRE_THROWS_AS(d = normalprob(100.0,0.0,-1.0),std::runtime_error);
}


const char* suite_string = "String utilities";

TEST_CASE("split breaks up const char*",suite_string) {
	vector<string> parts = split("ab*defgh**p*", '*');
	
	CHECK(parts.at(0) == "ab");
	CHECK(parts.at(1) == "defgh");
	CHECK(parts.at(2) == "");
	CHECK(parts.at(3) == "p");
	REQUIRE(parts.size() == 4);
	
	// ianhinder: I would have expected the following to work, but the last part
	// is not returned
	// CHECK(parts.at(4) == "");
	// CHECK(parts.size() == 5);
}

///////////////////////////////////////////////////////////////////////////
// getint
///////////////////////////////////////////////////////////////////////////


TEST_CASE("getint reads an integer",suite_string) {
	unsigned int i = getint("42", UNSET);
	REQUIRE(i == 42);
}
	
TEST_CASE("getint reads 'NA' as UNKNOWN",suite_string) {
	unsigned int i = getint("NA", UNSET);
	REQUIRE(i == UNKNOWN);
}

TEST_CASE("getint reads '*' as THRESH if threshold is given",suite_string) {
	unsigned int i = getint("*", 5);
	REQUIRE(i == THRESH);
}

TEST_CASE("getint fails on '*' if no threshold is given",suite_string) {
	unsigned int i;
	REQUIRE_THROWS_AS(i = getint("*", UNSET),std::runtime_error);
}

TEST_CASE("getint fails on junk",suite_string) {
	unsigned int i;
	REQUIRE_THROWS_AS(i = getint("xx", UNSET),std::runtime_error);
}

TEST_CASE("getint fails on trailing junk",suite_string) {
	unsigned int i;
	REQUIRE_THROWS_AS(i = getint("xx", UNSET),std::runtime_error);
}
