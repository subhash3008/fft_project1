#include "fft_test.h"

const double PI = 2 * acos(0.0);

// x = sin(2 * PI * frequency)
void generateSine(std::vector<float>& data, int frequency, int noOfSamples) {
  for (int i = 0; i < noOfSamples; ++i) {
    // float sampleVal = sin(2 * PI * frequency * i / static_cast<float>(samplesPerSecond));
    float sampleVal = 100 * sin((2 * PI / 180) * i * frequency);
    data.push_back(sampleVal);
  }
}

void runFftTest(std::vector<float>& data) {
  int noOfSamples = 1000;
  int frequency = 10;

  generateSine(data, frequency, noOfSamples);

  addImaginaryComponents(data);

  generatePlot(data, "testInputSignal.dat");

  padZeroTillPowerOfTwo(data);

  FFT(data, data.size() / 2 , 1);

  generatePlot(data, "testFftSignal.dat", true);

  FFT(data, data.size() / 2, -1);

  std::cout << "Took inverse of the data." << std::endl;

  erasePaddedZeroes(data, noOfSamples);

  generatePlot(data, "testReverseFftSignal.dat");

}