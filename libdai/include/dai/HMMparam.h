#ifndef __defined_libdai_hmmparam_h
#define __defined_libdai_hmmparam_h

#include <vector>
#include <dai/alldai.h>


namespace dai{

class HMMparam{
private:
	std::vector<Var> trans, obs;
	Permute permTrans, permObs;

public:

	std::vector<Factor> init, dist;

	HMMparam(const char* filename);

	void printHMMparam(const char* filename);

	void saveHMMparam(const char* filename);
};

}

#endif
