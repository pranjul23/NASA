
#ifndef __defined_libdai_observation_h
#define __defined_libdai_observation_h

#include <vector>

class Observation{
private:
	std::vector<std::vector<std::pair<size_t, size_t> > > data;

public:
	Observation(const char* filename);

	std::vector<std::vector<std::pair<size_t, size_t> > >& getData(){
		return data;
	}
};

#endif

