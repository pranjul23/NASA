#ifndef __defined_libdai_trainhsmm_h
#define __defined_libdai_trainhsmm_h


#include <vector>
#include <dai/alldai.h>


namespace dai{

class TrainHSMM{
	std::vector<std::vector<std::pair<size_t, size_t> > > data;

public:
	TrainHSMM(const char* filename);
	void train(const char* filename, size_t ID, size_t max_num_iter);

	//version of the above function for different HSMM model
	//has only prior for hidden state and at the last step does not draw new duration variable
	//Ntrain is the number of examples to use for training (-1 means use all of them)
	//Ltrain is the max length of sequences to use during training (-1 means use whole length)
	void train(const char* filename, size_t ID, size_t max_num_iter, int dummy, int Ntrain = -1, int Ltrain = -1);
};

}
#endif
