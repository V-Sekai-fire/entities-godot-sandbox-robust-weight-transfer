// Host build of the self tests: no engine, no RISC-V toolchain.
#include "../src/property_test.hpp"

int main() {
	const int failed = rwt::run_all_tests();
	std::cout << (failed == 0 ? "All tests passed!" : "Some tests failed.") << std::endl;
	return failed == 0 ? 0 : 1;
}
