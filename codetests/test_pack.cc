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

TEST_CASE("Pack can store and read back a double", tag_pack) {
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
TEST_CASE("Pack can store and read back a string", tag_pack) {
	packinit(0);
	pack(std::string("Hello, world"));
	REQUIRE(packsize() == 13);
	packinit(1);
	CHECK(packsize() == 0);
	std::string s;
	unpack(s);	
	REQUIRE(s == "Hello, world");
	REQUIRE(packsize() == 13);
}


TEST_CASE("Pack can store and read back a vector of unsigned int", tag_pack) {
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
TEST_CASE("Pack can store and read back a vector of int", tag_pack) {
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
TEST_CASE("Pack can store and read back a vector of unsigned short", tag_pack) {
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
TEST_CASE("Pack can store and read back a vector of double", tag_pack) {
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
TEST_CASE("Pack can store and read back a vector of vector of unsigned int",
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

TEST_CASE("Pack can store and read back a vector of vector of double",
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

TEST_CASE("Pack can store and read back a vector of vector of float",
					tag_pack) {
	std::vector<std::vector<float>> vec;	
	vec.resize(2);
	vec[0].push_back(1.f);
	vec[0].push_back(2.f);
	vec[1].push_back(-3.f);
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
	REQUIRE(vec[0][0] == 1.f);
	REQUIRE(vec[0][1] == 2.f);
	REQUIRE(vec[1][0] == -3.f);
}

TEST_CASE("Pack can store and read back a vector^3 of unsigned int",
					tag_pack) {
	std::vector<std::vector<std::vector<unsigned int>>> vec;	
	vec.resize(2);
	vec[0].resize(1);
	vec[1].resize(3);
	vec[0][0].push_back(1u);
	vec[1][0].push_back(2u);
	vec[1][1].push_back(3u);
	vec[1][2].push_back(4u);
	packinit(0);
	pack(vec);
	REQUIRE(packsize() == 11);
	packinit(6);
	CHECK(packsize() == 0);
	vec.resize(0);
	unpack(vec);	
	REQUIRE(packsize() == 11);
	REQUIRE(vec.size() == 2u);
	REQUIRE(vec[0].size() == 1u);
	REQUIRE(vec[1].size() == 3u);
	REQUIRE(vec[0][0].size() == 1u);
	REQUIRE(vec[1][0].size() == 1u);
	REQUIRE(vec[1][1].size() == 1u);
	REQUIRE(vec[1][2].size() == 1u);
	REQUIRE(vec[0][0][0] == 1u);
	REQUIRE(vec[1][0][0] == 2u);
	REQUIRE(vec[1][1][0] == 3u);
	REQUIRE(vec[1][2][0] == 4u);
}
TEST_CASE("Pack can store and read back a vector^3 of double",
					tag_pack) {
	std::vector<std::vector<std::vector<double>>> vec;	
	vec.resize(2);
	vec[0].resize(1);
	vec[1].resize(3);
	vec[0][0].push_back(1.);
	vec[1][0].push_back(2.);
	vec[1][1].push_back(-3.);
	vec[1][2].push_back(4.);
	packinit(0);
	pack(vec);
	REQUIRE(packsize() == 11);
	packinit(6);
	CHECK(packsize() == 0);
	vec.resize(0);
	unpack(vec);	
	REQUIRE(packsize() == 11);
	REQUIRE(vec.size() == 2u);
	REQUIRE(vec[0].size() == 1u);
	REQUIRE(vec[1].size() == 3u);
	REQUIRE(vec[0][0].size() == 1u);
	REQUIRE(vec[1][0].size() == 1u);
	REQUIRE(vec[1][1].size() == 1u);
	REQUIRE(vec[1][2].size() == 1u);
	REQUIRE(vec[0][0][0] == 1.);
	REQUIRE(vec[1][0][0] == 2.);
	REQUIRE(vec[1][1][0] == -3.);
	REQUIRE(vec[1][2][0] == 4.);
}

TEST_CASE("Pack can store and read back a string vector", tag_pack) {
	packinit(0);
	std::vector<std::string> s;
	s.push_back("Hello");
	s.push_back("world");
	pack(s);
	REQUIRE(packsize() == 13);
	packinit(13);
	CHECK(packsize() == 0);
	s.resize(0);
	unpack(s);	
	REQUIRE(s.size() == 2u);
	REQUIRE(s[0] == "Hello");
	REQUIRE(s[1] == "world");
	REQUIRE(packsize() == 13);
}
TEST_CASE("Pack can store and read back an AREA vector", tag_pack) {
	packinit(0);
	std::vector<AREA> av;
	AREA a;
	a.code = "code";
	a.region = 1u;
	a.x = 1.;
	a.y = 2.;
	a.agepop = vector<unsigned int>{15,45,75};
	a.pop = vector<unsigned int>{2,8,4};
	a.covar = vector<double>{1.,2.,3.};
	a.ind = vector<vector<unsigned int>>{
		vector<unsigned int>{1,2,3},
		vector<unsigned int>{4,5}};
	av.push_back(a);
	a.code = "recode";
	a.region = 2u;
	a.x = 3.;
	a.y = 4.;
	a.agepop = vector<unsigned int>{18,44,78};
	a.pop = vector<unsigned int>{1,7,3};
	a.covar = vector<double>{-1.,-2.,-3.};
	a.ind = vector<vector<unsigned int>>{
		vector<unsigned int>{6,7,8},
		vector<unsigned int>{9,10}};	
	av.push_back(a);
	pack(av);
	size_t size = packsize();
	packinit(size);
	av.resize(0);
	unpack(av);
	REQUIRE(av.size() == 2u);
	REQUIRE(av[0].code == "code");
	REQUIRE(av[1].code == "recode");
	REQUIRE(av[1].ind.size() == 2);
	REQUIRE(av[1].ind[1][1] == 10);
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
void pack(const vector< vector <FEV> > &vec, unsigned int fedivmin, unsigned int fedivmax);
void pack(const vector <AREA> &vec);
void pack(const vector <REGION> &vec);
void pack(const vector <DEMOCAT> &vec);
void pack(const vector <vector <EVREF> > &vec);

#endif
