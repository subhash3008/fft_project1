#ifndef FFT_UTILS_H
#define FFT_UTILS_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <stdint.h>

void swap (std::vector<float>& data, unsigned long first, unsigned long second);
void generatePlot(std::vector<float>& data, const std::string& fileName, bool isComplexPlot = false);
void addImaginaryComponents(std::vector<float>& data);
void padZeroTillPowerOfTwo(std::vector<float>& data);
void erasePaddedZeroes(std::vector<float>& data, int originalSize);
void applyTransmissionEqn(std::vector<float>& data);
void FFT (std::vector<float>& data, unsigned long number_of_complex_samples, int isForwardFFT);
void logData(std::vector<float>& data);

#endif