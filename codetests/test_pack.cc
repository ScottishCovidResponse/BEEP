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
