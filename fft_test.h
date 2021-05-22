#ifndef FFT_TEST_H
#define FFT_TEST_H

#include "fft_utils.h"

void generateSine(std::vector<float>& data, int frequency, int noOfSamples);
void runFftTest(std::vector<float>& data);

#endif