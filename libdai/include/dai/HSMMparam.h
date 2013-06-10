#ifndef __defined_libdai_hsmmparam_h
#define __defined_libdai_hsmmparam_h


#include <vector>
#include <dai/alldai.h>


namespace dai{

class HSMMparam{
private:
	std::vector<Var> trans, durat, obs;
	Permute permDurat, permTrans, permObs;

public:

	std::vector<Factor> init, dist;

	HSMMparam(const char* filename);

	void saveHSMMparam(const char* filename);

	void printHSMMparam(const char* filename);
};

}

#endif
