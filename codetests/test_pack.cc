#include "../catch.hpp"

#include "../pack.hh"

const char* tag_pack = "[pack]";
TEST_CASE("Pack initializes to zero size", tag_pack) {
	packinit(0);
	REQUIRE(packsize() == 0);
}
TEST_CASE("Pack can store and read back an unsigned int", tag_pack) {
	packinit(0);
	pack(1u);
	REQUIRE(packsize() == 1);
	packinit(1);
	CHECK(packsize() == 0);
	unsigned int i=0;
	unpack(i);	
	REQUIRE(i == 1u);
	REQUIRE(packsize() == 1);
}
TEST_CASE("Pack can store and read back two unsigned ints", tag_pack) {
	packinit(0);
	pack(1u);
	pack(2u);
	REQUIRE(packsize() == 2);
	// Unpack
	packinit(2);
	unsigned int i=0;
	unpack(i);	
	REQUIRE(i == 1u);
	unpack(i);	
	REQUIRE(i == 2u);
	REQUIRE(packsize() == 2);
}
TEST_CASE("Pack can store and read back large unsigned ints", tag_pack) {
	packinit(0);
	pack(std::numeric_limits<unsigned int>::max());
	pack(std::numeric_limits<unsigned int>::max()-1u);
	pack(std::numeric_limits<unsigned int>::max()-2u);
	// Unpack
	packinit(3);
	unsigned int i=0;
	unpack(i);	
	REQUIRE(i == std::numeric_limits<unsigned int>::max());
	unpack(i);	
	REQUIRE(i == std::numeric_limits<unsigned int>::max()-1u);
	unpack(i);	
	REQUIRE(i == std::numeric_limits<unsigned int>::max()-2u);
}

TEST_CASE("Pack can store and read back an unsigned short", tag_pack) {
	packinit(0);
	pack((unsigned short) 1);
	REQUIRE(packsize() == 1);
	packinit(1);
	CHECK(packsize() == 0);
	unsigned short i=0;
	unpack(i);	
	REQUIRE(i == (unsigned short) 1);
	REQUIRE(packsize() == 1);
}
TEST_CASE("Pack can store and read back two unsigned shorts", tag_pack) {
	packinit(0);
	pack((unsigned short) 1);
	pack((unsigned short) 2);
	REQUIRE(packsize() == 2);
	// Unpack
	packinit(2);
	unsigned short i=0;
	unpack(i);	
	REQUIRE(i == (unsigned short) 1);
	unpack(i);	
	REQUIRE(i == (unsigned short) 2);
	REQUIRE(packsize() == 2);
}
TEST_CASE("Pack can store and read back large unsigned shorts", tag_pack) {
	packinit(0);
	pack(std::numeric_limits<unsigned short>::max());
	pack((unsigned short) (std::numeric_limits<unsigned short>::max()-1));
	pack((unsigned short) (std::numeric_limits<unsigned short>::max()-2));
	// Unpack
	packinit(3);
	unsigned short i=0;
	unpack(i);	
	REQUIRE(i == std::numeric_limits<unsigned short>::max());
	unpack(i);	
	REQUIRE(i == std::numeric_limits<unsigned short>::max()-(unsigned short) 1);
	unpack(i);	
	REQUIRE(i == std::numeric_limits<unsigned short>::max()-(unsigned short) 2);
}

TEST_CASE("Pack can store and read back an double", tag_pack) {
	packinit(0);
	pack(1.);
	REQUIRE(packsize() == 1);
	packinit(1);
	CHECK(packsize() == 0);
	double i=0;
	unpack(i);	
	REQUIRE(i == 1.);
	REQUIRE(packsize() == 1);
}
TEST_CASE("Pack can store and read back two doubles", tag_pack) {
	packinit(0);
	pack(1.);
	pack(2.);
	REQUIRE(packsize() == 2);
	// Unpack
	packinit(2);
	double i=0;
	unpack(i);	
	REQUIRE(i == 1.);
	unpack(i);	
	REQUIRE(i == 2.);
	REQUIRE(packsize() == 2);
}
TEST_CASE("Pack can store and read back an vector of unsigned int", tag_pack) {
	std::vector<unsigned int> vec;	
	vec.push_back(1u);
	vec.push_back(2u);
	packinit(0);
	pack(vec);
	REQUIRE(packsize() == 3);
	packinit(3);
	CHECK(packsize() == 0);
	vec.resize(0);
	unpack(vec);	
	REQUIRE(packsize() == 3);
	REQUIRE(vec.size() == 2u);
	REQUIRE(vec[0] == 1u);
	REQUIRE(vec[1] == 2u);
}
TEST_CASE("Pack can store and read back an vector of int", tag_pack) {
	std::vector<int> vec;	
	vec.push_back(1);
	vec.push_back(2);
	vec.push_back(-2);
	packinit(0);
	pack(vec);
	REQUIRE(packsize() == 4);
	packinit(4);
	CHECK(packsize() == 0);
	vec.resize(0);
	unpack(vec);	
	REQUIRE(packsize() == 4);
	REQUIRE(vec.size() == 3u);
	REQUIRE(vec[0] == 1);
	REQUIRE(vec[1] == 2);
	REQUIRE(vec[2] == -2);
}
TEST_CASE("Pack can store and read back an vector of unsigned short", tag_pack) {
	std::vector<unsigned short> vec;	
	vec.push_back(1u);
	vec.push_back(2u);
	packinit(0);
	pack(vec);
	REQUIRE(packsize() == 3);
	packinit(3);
	CHECK(packsize() == 0);
	vec.resize(0);
	unpack(vec);	
	REQUIRE(packsize() == 3);
	REQUIRE(vec.size() == 2u);
	REQUIRE(vec[0] == 1u);
	REQUIRE(vec[1] == 2u);
}
TEST_CASE("Pack can store and read back an vector of double", tag_pack) {
	std::vector<double> vec;	
	vec.push_back(1.);
	vec.push_back(2.);
	packinit(0);
	pack(vec);
	REQUIRE(packsize() == 3);
	packinit(3);
	CHECK(packsize() == 0);
	vec.resize(0);
	unpack(vec);	
	REQUIRE(packsize() == 3);
	REQUIRE(vec.size() == 2u);
	REQUIRE(vec[0] == 1.);
	REQUIRE(vec[1] == 2.);
}
TEST_CASE("Pack can store and read back an vector of vector of unsigned int",
					tag_pack) {
	std::vector<std::vector<unsigned int>> vec;	
	vec.resize(2);
	vec[0].push_back(1u);
	vec[0].push_back(2u);
	vec[1].push_back(3u);
	packinit(0);
	pack(vec);
	REQUIRE(packsize() == 6);
	packinit(6);
	CHECK(packsize() == 0);
	vec.resize(0);
	unpack(vec);	
	REQUIRE(packsize() == 6);
	REQUIRE(vec.size() == 2u);
	REQUIRE(vec[0].size() == 2u);
	REQUIRE(vec[1].size() == 1u);
	REQUIRE(vec[0][0] == 1u);
	REQUIRE(vec[0][1] == 2u);
	REQUIRE(vec[1][0] == 3u);
}
TEST_CASE("Pack can store and read back an vector of vector of double",
					tag_pack) {
	std::vector<std::vector<double>> vec;	
	vec.resize(2);
	vec[0].push_back(1.);
	vec[0].push_back(2.);
	vec[1].push_back(-3.);
	packinit(0);
	pack(vec);
	REQUIRE(packsize() == 6);
	packinit(6);
	CHECK(packsize() == 0);
	vec.resize(0);
	unpack(vec);	
	REQUIRE(packsize() == 6);
	REQUIRE(vec.size() == 2u);
	REQUIRE(vec[0].size() == 2u);
	REQUIRE(vec[1].size() == 1u);
	REQUIRE(vec[0][0] == 1.);
	REQUIRE(vec[0][1] == 2.);
	REQUIRE(vec[1][0] == -3.);
}

#if 0
TEST_CASE("Unpack handles buffer overflow", tag_pack) {
	packinit();
	emsg_throws = true;
	for (long i=0; i<MAX_NUMBERS; ++i)
		pack(1u);
	REQUIRE(packsize() == MAX_NUMBERS);
	packinit();
	unsigned int n;
	for (long i=0; i<MAX_NUMBERS; ++i)
		unpack(n);
	REQUIRE(packsize() == MAX_NUMBERS);
	CHECK_THROWS_AS(unpack(n),std::runtime_error);
}

// Other cases to address
void pack(const string &vec);
void pack(const vector <unsigned short> &vec);
void pack(const vector <int> &vec);
void pack(const vector <double> &vec);
void pack(const vector< vector <unsigned int> > &vec);
void pack(const vector< vector <double> > &vec);
void pack(const vector< vector <float> > &vec);
void pack(const vector< vector< vector <double> > > &vec);
void pack(const vector <string> &vec);
void pack(const vector< vector <FEV> > &vec, unsigned int fedivmin, unsigned int fedivmax);
void pack(const vector <AREA> &vec);
void pack(const vector <REGION> &vec);
void pack(const vector <DEMOCAT> &vec);
void pack(const vector <vector <EVREF> > &vec);
void pack(const unsigned short *vec, unsigned int imax);
void pack(float **vec, unsigned int imax, unsigned int jmax);
void pack(const vector < vector <vector <unsigned int> > > &vec);

#endif
