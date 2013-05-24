#ifndef __defined_libdai_trainhsmm_h
#define __defined_libdai_trainhsmm_h


#include <vector>
#include <dai/alldai.h>


namespace dai{

class TrainHSMM{
	std::vector<std::vector<std::pair<size_t, size_t> > > data;

public:
	TrainHSMM(const char* filename);
	void train();
};

}
#endif
