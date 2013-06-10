
#ifndef __defined_libdai_testhmm_h
#define __defined_libdai_testhmm_h


#include <vector>
#include <dai/alldai.h>

namespace dai{

class TestHMM{
	std::vector<std::vector<std::pair<size_t, size_t> > > test_data;
public:
	TestHMM(const char* filename);
	void test_loglik(const char* filename, size_t ID);
	void test_marginal(const char* filename, size_t ID);
};

}

#endif
