#pragma once

#include "StatsLibrary.hpp"
#include <assert.h>
#include <stdlib.h> 

using namespace std;

int StatsLibrary::sample_multinomial(vector<double>& probs) {
	int K = probs.size();
    double sum_so_far = 0;
    int sample = -1;
    double rand_num = ((double)rand()) / RAND_MAX;
    for (int k = 0; k < K; ++k) {
        sum_so_far += probs[k];
        if (rand_num <= sum_so_far || k == K - 1) {
            sample = k;
            break;
        }
    }
    assert(sample != -1);
    return sample;
}
