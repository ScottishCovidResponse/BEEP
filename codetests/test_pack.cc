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
#if 0
TEST_CASE("Pack handles buffer overflow", tag_pack) {
	packinit();
	emsg_throws = true;
	for (long i=0; i<MAX_NUMBERS; ++i)
		pack(1u);
	REQUIRE(packsize() == MAX_NUMBERS);
	CHECK_THROWS_AS(pack(2u),std::runtime_error);
}
TEST_CASE("Unpack handles buffer overflow", tag_pack) {
	packinit();
	emsg_throws = true;
	for (long i=0; i<MAX_NUMBERS; ++i)
		pack(1u);
	REQUIRE(packsize() == MAX_NUMBERS);
	packinit();
	emsg_throws = true;
	unsigned int n;
	for (long i=0; i<MAX_NUMBERS; ++i)
		unpack(n);
	REQUIRE(packsize() == MAX_NUMBERS);
	CHECK_THROWS_AS(unpack(n),std::runtime_error);
}
#endif
