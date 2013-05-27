#ifndef __defined_libdai_trainhmm_h
#define __defined_libdai_trainhmm_h


#include <vector>
#include <dai/alldai.h>


namespace dai{

class TrainHMM{
	std::vector<std::vector<std::pair<size_t, size_t> > > data;
public:
	TrainHMM(const char* filename);
	void train();
};

}
#endif
