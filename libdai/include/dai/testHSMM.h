
#ifndef __defined_libdai_testhsmm_h
#define __defined_libdai_testhsmm_h


#include <vector>
#include <dai/alldai.h>

namespace dai{

class TestHSMM{
	std::vector<std::vector<std::pair<size_t, size_t> > > test_data;

public:
	TestHSMM(const char* filename);
	void test();
};

}

#endif
