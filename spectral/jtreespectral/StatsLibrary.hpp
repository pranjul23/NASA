#pragma once

#include <vector>
#include <sstream>
#include <cmath>

using namespace std;

class StatsLibrary
{
public:
	static int sample_multinomial(vector<double>& probs);
};



template <class T>
string to_string (T& t)
{
	std::stringstream ss;
	ss << t;
	return ss.str();
}
