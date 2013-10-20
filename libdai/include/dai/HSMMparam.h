#ifndef __defined_libdai_hsmmparam_h
#define __defined_libdai_hsmmparam_h


#include <vector>
#include <dai/alldai.h>


namespace dai{

class HSMMparam{
private:
	std::vector<Var> trans, dur, durat, obs;
	Permute permDur, permDurat, permTrans, permObs;

public:

	std::vector<Factor> init, dist;

	HSMMparam(const char* filename);

	void printHSMMparam(const char* filename);

	void saveHSMMparam(const char* filename);

	// version of the above functions for different HSMM model
	HSMMparam(const char* filename, int dummy);
	void printHSMMparam(const char* filename, int dummy);
	void saveHSMMparam(const char* filename, int dummy);
};

}

#endif
