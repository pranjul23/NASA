#include <dai/observation.h>
#include <fstream>
#include <iostream>

Observation::Observation(const char* filename){

	std::ifstream file(filename);
	size_t var, value, numObs, dataLen;

	file >> dataLen;
	data.reserve(dataLen);

	for(size_t i=0; i<dataLen; i++){
		file >> numObs;
		data.push_back(std::vector<std::pair<size_t, size_t> >(numObs, std::make_pair(0, 0)));

		for(size_t j=0; j<numObs; j++){
			file >> var;
			data.back().at(j).first = var;
		}

		for(size_t j=0; j<numObs; j++){
			file >> value;
			data.back().at(j).second = value;
		}

	}
/*
		std::cout << " data.size = " << data.size() << "\n";
		for(int i=0; i<data.size();i++){
			std::cout << " data[i].size = " << data[i].size() << "\n";
			for(int j=0; j<data[i].size(); j++){
				std::cout << "data[i][j].first = " << data[i][j].first  << "\n";
				std::cout << "data[i][j].second = " << data[i][j].second  << "\n";
			}
		}
		*/
}


