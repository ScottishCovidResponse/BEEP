#include "../catch.hpp"

#include "../pack.hh"

const char* tag_pack = "[pack]";
TEST_CASE("Pack initializes to zero size", tag_pack) {
	packinit();
	REQUIRE(packsize() == 0);
}
TEST_CASE("Pack can store and read back an unsigned int", tag_pack) {
	packinit();
	pack(1u);
	REQUIRE(packsize() == 1);
	packinit();
	CHECK(packsize() == 0);
	unsigned int i=0;
	unpack(i);	
	REQUIRE(i == 1u);
	REQUIRE(packsize() == 1);
}
TEST_CASE("Pack can store and read back two unsigned ints", tag_pack) {
	packinit();
	pack(1u);
	pack(2u);
	REQUIRE(packsize() == 2);
	packinit();
	unsigned int i=0;
	unpack(i);	
	REQUIRE(i == 1u);
	unpack(i);	
	REQUIRE(i == 2u);
	REQUIRE(packsize() == 2);
}
TEST_CASE("Pack can store and read back large unsigned ints", tag_pack) {
	packinit();
	pack(std::numeric_limits<unsigned int>::max());
	pack(std::numeric_limits<unsigned int>::max()-1u);
	pack(std::numeric_limits<unsigned int>::max()-2u);
	packinit();
	unsigned int i=0;
	unpack(i);	
	REQUIRE(i == std::numeric_limits<unsigned int>::max());
	unpack(i);	
	REQUIRE(i == std::numeric_limits<unsigned int>::max()-1u);
	unpack(i);	
	REQUIRE(i == std::numeric_limits<unsigned int>::max()-2u);
}

TEST_CASE("Pack handles buffer overflow", tag_pack) {
	packinit();
	emsg_throws = true;
	for (long i=0; i<MAX_NUMBERS; ++i)
		pack(1u);
	REQUIRE(packsize() == MAX_NUMBERS);
	CHECK_THROWS_AS(pack(2u),std::runtime_error);
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
void pack(const unsigned short num);
void pack(const double num);
void pack(const string &vec);
void pack(const vector <unsigned int> &vec);
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
