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

TEST_CASE("normal after seed set to 0 returns result expected for std::mt19937",
					tag_random) {
	sran(0);
	double r = normal(1.,1.);
	
	REQUIRE(r == Approx(1.5708607051));
}

TEST_CASE("gammasamp with invalid bounds throws",
					tag_random) {
	sran(0);
	double r;
	emsg_throws = true;
	
	CHECK_THROWS_AS(r = gammasamp(-1.,1.),std::runtime_error);
	CHECK_NOTHROW(r = gammasamp(0.,1.));
	CHECK_THROWS_AS(r = gammasamp(1.,-1.),std::runtime_error);
	CHECK_NOTHROW(r = gammasamp(1.,0.));
}

TEST_CASE("gammasamp after seed set to 0 returns result expected for small a with std::mt19937",
					tag_random) {
	sran(0);
	double r = gammasamp(0.5,1.);
	
	REQUIRE(r == Approx(0.5648378411));
}

TEST_CASE("gammasamp after seed set to 0 returns result expected for large a with std::mt19937",
					tag_random) {
	sran(0);
	double r = gammasamp(2.,1.);
	
	REQUIRE(r == Approx(2.5176090456));
}
TEST_CASE("gammasamp after seed set to 0 returns result expected for a at boundary with std::mt19937",
					tag_random) {
	sran(0);
	double r = gammasamp(1.,1.);
	
	REQUIRE(r == Approx(1.2498384326));
}

//////////////////////////////////
// Probability Distributions
//////////////////////////////////
const char* tag_distributions = "[distributions]";

TEST_CASE("logged normalprob throws out of domain",
					tag_distributions) {
	double d;
	emsg_throws = true;
	CHECK_THROWS_AS(d = normalprob(100.0,0.0,0.0),std::runtime_error);
	REQUIRE_THROWS_AS(d = normalprob(100.0,0.0,-1.0),std::runtime_error);
}

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

TEST_CASE("logged lognormprob throws out of domain",
					tag_distributions) {
	double d;
	emsg_throws = true;
	CHECK_THROWS_AS(d = lognormprob(0.0,0.0,1.0),std::runtime_error);
	CHECK_THROWS_AS(d = lognormprob(-1.0,0.0,1.0),std::runtime_error);
	CHECK_THROWS_AS(d = lognormprob(100.0,0.0,0.0),std::runtime_error);
	CHECK_THROWS_AS(d = lognormprob(100.0,0.0,-1.0),std::runtime_error);
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

TEST_CASE("logged gammaprob throws out of domain",
					tag_distributions) {
	double d;
	emsg_throws = true;
	CHECK_NOTHROW(d = gammaprob(1.0,1.0,1.0));
	CHECK_THROWS_AS(d = gammaprob(0.0,1.0,1.0),std::runtime_error);
	CHECK_THROWS_AS(d = gammaprob(-1.0,1.0,1.0),std::runtime_error);
	CHECK_THROWS_AS(d = gammaprob(1.0,-1.0,1.0),std::runtime_error);
	CHECK_THROWS_AS(d = gammaprob(1.0,0.0,1.0),std::runtime_error);
	CHECK_THROWS_AS(d = gammaprob(1.0,1.0,-1.0),std::runtime_error);
	CHECK_THROWS_AS(d = gammaprob(1.0,1.0,0.0),std::runtime_error);
}

TEST_CASE("logged gammaprob at (x=1,a=1,b=1) is -1.",
					tag_distributions) {
	REQUIRE(gammaprob(1.0,1.0,1.0) == Approx( -1.0 ));
}

//>>> import scipy.stats as stats
//>>> import math
//>>> x=3;a=2;b=1;stats.gamma(a).logpdf(x*b)+math.log(b)
//-1.9013877113318902
//(%i) load(distrib)$
//(%i) log(pdf_gamma(3.0,2.0,1.0));
//(%o)    - 1.90138771133189
TEST_CASE("logged gammaprob at (x=3,a=2,b=1) is -1.90138771133189",
					tag_distributions) {
	REQUIRE(gammaprob(3.0,2.0,1.0) == Approx( -1.90138771133189 ));
}

// Note different convention for b parameter in maxima
//>>> x=1;a=1;b=2;stats.gamma(a).logpdf(x*b)+math.log(b)
//-1.3068528194400546
//(%i) load(distrib)$
//(%i) log(pdf_gamma(1.0,1.0,1./2.0));
//(%o)    - 1.306852819440054
TEST_CASE("logged gammaprob at (x=1,a=1,b=2) is -1.3068528194400546",
					tag_distributions) {
	REQUIRE(gammaprob(1.0,1.0,2.0) == Approx( -1.3068528194400546 ));
}

TEST_CASE("logged gammaprob in right wing is tiny",
					tag_distributions) {
	REQUIRE(gammaprob(200.0,4.0,1.0) < -100.0);
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

TEST_CASE("filebasename with no subdirectory finds same name",tag_string) {
  REQUIRE(filebasename("file") == "file");
}

TEST_CASE("filebasename with one subdirectory finds last name",tag_string) {
  REQUIRE(filebasename("dir/file") == "file");
}

TEST_CASE("filebasename with two subdirectories finds last name",tag_string) {
  REQUIRE(filebasename("dir1/dir2/file") == "file");
}

TEST_CASE("stringhasending with short string doesn't match",tag_string) {
  REQUIRE(!stringhasending("test.img","lotsandlotsamdlots"));
}

TEST_CASE("stringhasending with non-matching string doesn't match",tag_string) {
  REQUIRE(!stringhasending("test.img","lots"));
}
TEST_CASE("stringhasending with matching string does match",tag_string) {
  REQUIRE(stringhasending("test.img","img"));
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

TEST_CASE("getint fails on negative",tag_string) {
	unsigned int i;
	REQUIRE_THROWS_AS(i = getint("-42", UNSET),std::runtime_error);
}
TEST_CASE("getint fails on huge",tag_string) {
	unsigned int i;
	REQUIRE_THROWS_AS(i = getint("55555555555555555555", UNSET),
										std::runtime_error);
}

TEST_CASE("getint fails on junk",tag_string) {
	unsigned int i;
	REQUIRE_THROWS_AS(i = getint("xx", UNSET),std::runtime_error);
}

TEST_CASE("getint fails on trailing junk",tag_string) {
	unsigned int i;
	REQUIRE_THROWS_AS(i = getint("xx", UNSET),std::runtime_error);
}
