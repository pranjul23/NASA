
#ifndef __defined_libdai_testhsmm_h
#define __defined_libdai_testhsmm_h


#include <vector>
#include <string>
#include <dai/alldai.h>

namespace dai{

class TestHSMM{
	std::vector<std::vector<std::pair<size_t, size_t> > > test_data;

public:
	TestHSMM(const char* filename);
	void test_loglik(const char* filename, size_t ID, std::string type);

	//same as above except that we use different HSMM model
	void test_loglik(const char* filename, size_t ID, int dummy);

	//same as above but with specification of type of testing
	//last parameter is for special hsmm testing (where training was run specified # of iter)
	void test_loglik(const char* filename, size_t ID, std::string type, int dummy, int ID2 =- 1);

	void test_marginal(const char* filename, size_t ID);
	void test_marginal_cut(const char* filename, size_t ID);
};

}

#endif
