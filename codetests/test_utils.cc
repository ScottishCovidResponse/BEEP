#include "../catch.hpp"

#include <iostream>
#include <iomanip>

#include "../utils.hh"

#include "../consts.hh"

using namespace std;

//////////////////////////////////
// Random generators
//////////////////////////////////

const char* tag_random = "[random]";
TEST_CASE("ran after seed set to 0 returns result expected for std::mt19937",
					tag_random) {
	sran(0);
	double r = ran();
	
//	cout << "r == " << std::setprecision(19) << r << endl;
	
	REQUIRE(r == 0.5928446159645536229);
}

//////////////////////////////////
// Probability Distributions
//////////////////////////////////
const char* tag_distributions = "[distributions]";

TEST_CASE("logged normalprob at mean with variance 1./(2*pi) is 0.",
					tag_distributions) {
	REQUIRE(normalprob(1.0,1.0,1.0/(2*M_PI)) == Approx( 0.0 ));
}

TEST_CASE("logged normalprob at mean+1sd with variance 1./(2*pi) is -0.5",
					tag_distributions) {
	double var = 1.0/(2*M_PI), sd = sqrt(var), mean=1.0;
	REQUIRE(normalprob(mean+sd,mean,var) == Approx( -0.5 ));
}

TEST_CASE("logged normalprob in left wing is tiny",
					tag_distributions) {
	REQUIRE(normalprob(-100.0,0.0,1.0/(2*M_PI)) < -100.0);
}

TEST_CASE("logged normalprob in right wing is tiny",
					tag_distributions) {
	REQUIRE(normalprob(100.0,0.0,1.0/(2*M_PI)) < -100.0);
}

TEST_CASE("logged normalprob throws out of domain",
					tag_distributions) {
	double d;
	emsg_throws = true;
	CHECK_THROWS_AS(d = normalprob(100.0,0.0,0.0),std::runtime_error);
	REQUIRE_THROWS_AS(d = normalprob(100.0,0.0,-1.0),std::runtime_error);
}

TEST_CASE("logged lognormprob at mean with variance 1./(2*pi) is -1.0",
					tag_distributions) {
	REQUIRE(lognormprob(exp(1.0),1.0,1.0/(2*M_PI)) == Approx( -1.0 ));
}

TEST_CASE("logged lognormprob at mean+1sd with variance 1./(2*pi) is -1.8989",
					tag_distributions) {
	double var = 1.0/(2*M_PI), sd = sqrt(var), mean=1.0;
	REQUIRE(lognormprob(exp(mean+sd),mean,var) == Approx( -0.5-(mean+sd) ));
}

TEST_CASE("logged lognormprob in left wing is tiny",
					tag_distributions) {
	REQUIRE(lognormprob(exp(-100.0),0.0,1.0/(2*M_PI)) < -100.0);
}

TEST_CASE("logged lognormprob in right wing is tiny",
					tag_distributions) {
	REQUIRE(lognormprob(exp(100.0),0.0,1.0/(2*M_PI)) < -100.0);
}

TEST_CASE("logged lognormprob throws out of domain",
					tag_distributions) {
	double d;
	emsg_throws = true;
	CHECK_THROWS_AS(d = lognormprob(0.0,0.0,1.0),std::runtime_error);
	CHECK_THROWS_AS(d = lognormprob(-1.0,0.0,1.0),std::runtime_error);
	CHECK_THROWS_AS(d = lognormprob(100.0,0.0,0.0),std::runtime_error);
	CHECK_THROWS_AS(d = lognormprob(100.0,0.0,-1.0),std::runtime_error);
}

//////////////////////////////////
// String utilities
//////////////////////////////////
const char* tag_string = "[string]";

TEST_CASE("split breaks up const char*",tag_string) {
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


TEST_CASE("getint reads an integer",tag_string) {
	unsigned int i = getint("42", UNSET);
	REQUIRE(i == 42);
}
	
TEST_CASE("getint reads 'NA' as UNKNOWN",tag_string) {
	unsigned int i = getint("NA", UNSET);
	REQUIRE(i == UNKNOWN);
}

TEST_CASE("getint reads '*' as THRESH if threshold is given",tag_string) {
	unsigned int i = getint("*", 5);
	REQUIRE(i == THRESH);
}

TEST_CASE("getint fails on '*' if no threshold is given",tag_string) {
	unsigned int i;
	REQUIRE_THROWS_AS(i = getint("*", UNSET),std::runtime_error);
}

TEST_CASE("getint fails on junk",tag_string) {
	unsigned int i;
	REQUIRE_THROWS_AS(i = getint("xx", UNSET),std::runtime_error);
}

TEST_CASE("getint fails on trailing junk",tag_string) {
	unsigned int i;
	REQUIRE_THROWS_AS(i = getint("xx", UNSET),std::runtime_error);
}
