#include "../catch.hpp"

#include "../pack.hh"

const char* tag_pack = "[pack]";
TEST_CASE("Pack initializes to zero size", tag_pack) {
	pack_init(0);
	REQUIRE(packsize() == 0);
}
TEST_CASE("Pack can store and read back an unsigned int", tag_pack) {
	pack_init(0);
	pack(1u);
	REQUIRE(packsize() == 1);
	pack_init(1);
	CHECK(packsize() == 0);
	unsigned int i=0;
	unpack(i);	
	REQUIRE(i == 1u);
	REQUIRE(packsize() == 1);
}
TEST_CASE("Pack can store and read back two unsigned ints", tag_pack) {
	pack_init(0);
	pack(1u);
	pack(2u);
	REQUIRE(packsize() == 2);
	// Unpack
	pack_init(2);
	unsigned int i=0;
	unpack(i);	
	REQUIRE(i == 1u);
	unpack(i);	
	REQUIRE(i == 2u);
	REQUIRE(packsize() == 2);
}
TEST_CASE("Pack can store and read back LARGE unsigned ints", tag_pack) {
	pack_init(0);
	pack(std::numeric_limits<unsigned int>::max());
	pack(std::numeric_limits<unsigned int>::max()-1u);
	pack(std::numeric_limits<unsigned int>::max()-2u);
	// Unpack
	pack_init(3);
	unsigned int i=0;
	unpack(i);	
	REQUIRE(i == std::numeric_limits<unsigned int>::max());
	unpack(i);	
	REQUIRE(i == std::numeric_limits<unsigned int>::max()-1u);
	unpack(i);	
	REQUIRE(i == std::numeric_limits<unsigned int>::max()-2u);
}

TEST_CASE("Pack can store and read back an unsigned short", tag_pack) {
	pack_init(0);
	pack((unsigned short) 1);
	REQUIRE(packsize() == 1);
	pack_init(1);
	CHECK(packsize() == 0);
	unsigned short i=0;
	unpack(i);	
	REQUIRE(i == (unsigned short) 1);
	REQUIRE(packsize() == 1);
}
TEST_CASE("Pack can store and read back two unsigned shorts", tag_pack) {
	pack_init(0);
	pack((unsigned short) 1);
	pack((unsigned short) 2);
	REQUIRE(packsize() == 2);
	// Unpack
	pack_init(2);
	unsigned short i=0;
	unpack(i);	
	REQUIRE(i == (unsigned short) 1);
	unpack(i);	
	REQUIRE(i == (unsigned short) 2);
	REQUIRE(packsize() == 2);
}
TEST_CASE("Pack can store and read back LARGE unsigned shorts", tag_pack) {
	pack_init(0);
	pack(std::numeric_limits<unsigned short>::max());
	pack((unsigned short) (std::numeric_limits<unsigned short>::max()-1));
	pack((unsigned short) (std::numeric_limits<unsigned short>::max()-2));
	// Unpack
	pack_init(3);
	unsigned short i=0;
	unpack(i);	
	REQUIRE(i == std::numeric_limits<unsigned short>::max());
	unpack(i);	
	REQUIRE(i == std::numeric_limits<unsigned short>::max()-(unsigned short) 1);
	unpack(i);	
	REQUIRE(i == std::numeric_limits<unsigned short>::max()-(unsigned short) 2);
}

TEST_CASE("Pack can store and read back a double", tag_pack) {
	pack_init(0);
	pack(1.);
	REQUIRE(packsize() == 1);
	pack_init(1);
	CHECK(packsize() == 0);
	double i=0;
	unpack(i);	
	REQUIRE(i == 1.);
	REQUIRE(packsize() == 1);
}
TEST_CASE("Pack can store and read back two doubles", tag_pack) {
	pack_init(0);
	pack(1.);
	pack(2.);
	REQUIRE(packsize() == 2);
	// Unpack
	pack_init(2);
	double i=0;
	unpack(i);	
	REQUIRE(i == 1.);
	unpack(i);	
	REQUIRE(i == 2.);
	REQUIRE(packsize() == 2);
}
TEST_CASE("Pack can store and read back a string", tag_pack) {
	pack_init(0);
	pack(std::string("Hello, world"));
	REQUIRE(packsize() == 13);
	pack_init(1);
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
	pack_init(0);
	pack(vec);
	REQUIRE(packsize() == 3);
	pack_init(3);
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
	pack_init(0);
	pack(vec);
	REQUIRE(packsize() == 4);
	pack_init(4);
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
	pack_init(0);
	pack(vec);
	REQUIRE(packsize() == 3);
	pack_init(3);
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
	pack_init(0);
	pack(vec);
	REQUIRE(packsize() == 3);
	pack_init(3);
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
	pack_init(0);
	pack(vec);
	REQUIRE(packsize() == 6);
	pack_init(6);
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
	pack_init(0);
	pack(vec);
	REQUIRE(packsize() == 6);
	pack_init(6);
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
	pack_init(0);
	pack(vec);
	REQUIRE(packsize() == 6);
	pack_init(6);
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
	pack_init(0);
	pack(vec);
	REQUIRE(packsize() == 11);
	pack_init(6);
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
	pack_init(0);
	pack(vec);
	REQUIRE(packsize() == 11);
	pack_init(6);
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
	pack_init(0);
	std::vector<std::string> s;
	s.push_back("Hello");
	s.push_back("world");
	pack(s);
	REQUIRE(packsize() == 13);
	pack_init(13);
	CHECK(packsize() == 0);
	s.resize(0);
	unpack(s);	
	REQUIRE(s.size() == 2u);
	REQUIRE(s[0] == "Hello");
	REQUIRE(s[1] == "world");
	REQUIRE(packsize() == 13);
}
TEST_CASE("Pack can store and read back an AREA vector", tag_pack) {
	pack_init(0);
	std::vector<Area> av;
	Area a;
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
	pack_init(size);
	av.resize(0);
	unpack(av);
	REQUIRE(av.size() == 2u);
	REQUIRE(av[0].code == "code");
	REQUIRE(av[1].code == "recode");
	REQUIRE(av[1].ind.size() == 2);
	REQUIRE(av[1].ind[1][1] == 10);
}
TEST_CASE("Pack can store and read back an REGION vector", tag_pack) {
	pack_init(0);
	std::vector<DataRegion> rv;
	DataRegion r;
	r.name = "name1";
	r.code = "code1";
	rv.push_back(r);
	r.name = "name2";
	r.code = "code2";
	rv.push_back(r);
	pack(rv);
	size_t size = packsize();
	pack_init(size);
	rv.resize(0);
	unpack(rv);
	REQUIRE(rv.size() == 2u);
	REQUIRE(rv[0].code == "code1");
	REQUIRE(rv[1].code == "code2");
	REQUIRE(rv[0].name == "name1");
	REQUIRE(rv[1].name == "name2");
}
TEST_CASE("Pack can store and read back an DEMOCAT vector", tag_pack) {
	pack_init(0);
	std::vector<DemographicCategory> rv;
	DemographicCategory r;
	r.name = "name1";
	r.value = std::vector<std::string>{"value1a","value1b"};
	rv.push_back(r);
	r.name = "name2";
	r.value = std::vector<std::string>{"value2a","value2b"};
	rv.push_back(r);
	pack(rv);
	size_t size = packsize();
	pack_init(size);
	rv.resize(0);
	unpack(rv);
	REQUIRE(rv.size() == 2u);
	REQUIRE(rv[0].name == "name1");
	REQUIRE(rv[1].name == "name2");
	REQUIRE(rv[0].value[0] == "value1a");
	REQUIRE(rv[0].value[1] == "value1b");
	REQUIRE(rv[1].value[0] == "value2a");
	REQUIRE(rv[1].value[1] == "value2b");
}
TEST_CASE("Pack can store and read back an EVREF vector of vector", tag_pack) {
	pack_init(0);
	std::vector<std::vector<EventRef>> rv;
	rv.resize(2);
	EventRef r;
	r.ind = 1;
	r.e = 2;
	rv[0].push_back(r);
	r.ind = 3;
	r.e = 4;
	rv[1].push_back(r);
	pack(rv);
	size_t size = packsize();
	pack_init(size);
	rv.resize(0);
	unpack(rv);
	REQUIRE(rv.size() == 2u);
	REQUIRE(rv[0][0].ind == 1);
	REQUIRE(rv[0][0].e == 2);
	REQUIRE(rv[1][0].ind == 3);
	REQUIRE(rv[1][0].e == 4);
}

