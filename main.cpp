// This code transforms the data stored in "alice.h" into a WAV-file
// and saves it as "alice.wav".

#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <stdint.h>
// #include "alice.h"
#include "abcd.h"

using namespace std;

#define FILE_SIZE 1000000

typedef struct WAV_HEADER
{
  /* RIFF Chunk Descriptor */
  uint8_t RIFF[4] = {'R', 'I', 'F', 'F'}; // RIFF Header Magic header
  uint32_t ChunkSize;                     // RIFF Chunk Size
  uint8_t WAVE[4] = {'W', 'A', 'V', 'E'}; // WAVE Header
  /* "fmt" sub-chunk */
  uint8_t fmt[4] = {'f', 'm', 't', ' '};    // FMT header
  uint32_t Subchunk1Size = 16;              // Size of the fmt chunk
  uint16_t AudioFormat = 1;                 // Audio format 1=PCM,6=mulaw,7=alaw,     257=IBM
                                            // Mu-Law, 258=IBM A-Law, 259=ADPCM
  uint16_t NumOfChan = 1;                   // Number of channels 1=Mono 2=Sterio
  uint32_t SamplesPerSec = 20000;           // Sampling Frequency in Hz
  uint32_t bytesPerSec = SamplesPerSec * 2; // bytes per second
  uint16_t blockAlign = 2;                  // 2=16-bit mono, 4=16-bit stereo
  uint16_t bitsPerSample = 16;              // Number of bits per sample
  /* "data" sub-chunk */
  uint8_t Subchunk2ID[4] = {'d', 'a', 't', 'a'}; // "data"  string
  uint32_t Subchunk2Size;                        // Sampled data length
} wav_hdr;

// Global Variables
std::vector<float> paddedData;

const double PI = 2*acos(0.0);

void swap (vector<float>& data, unsigned long first, unsigned long second) {
  double temp = data.at(first);
  data.at(first) = data.at(second);
  data.at(second) = temp;
}

// Generate dat file for gnu plot
void generatePlot(const std::string& fileName) {
  std::ofstream plot { fileName };
  for (int x = 0; x < paddedData.size(); x += 2) {
    plot << (x / 2) << " " << std::sqrt(std::pow(paddedData.at(x), 2) + pow(paddedData.at(x + 1), 2)) << std::endl;
  }
  plot.close();
}

/*
  isForwardFFT = 1 for FFT and -1 for Inverse FFT
*/
void FFT (vector<float>& data, unsigned long number_of_complex_samples, int isForwardFFT) {
  //variables for trigonometric recurrences
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;
  
  n = number_of_complex_samples * 2;

  //binary inversion (note that the indexes
  //start from 0 witch means that the
  //real part of the complex is on the even-indexes
  //and the complex part is on the odd-indexes
  j = 0;
  for (i = 0; i < n; i += 2) {
    if (j > i) {
      //swap the real part
      swap(data, j, i);
      //swap the complex part
      swap(data, j + 1, i + 1);
      // checks if the changes occurs in the first half
      // and use the mirrored effect on the second half
      // if ((j / 2) < (n / 4)) {
      //   //swap the real part
      //   swap(data, n - (i + 2), n - (j + 2));
      //   //swap the complex part
      //   swap(data, (n - (i + 2)) + 1, (n - (j + 2)) + 1);
      // }
    }
    m = n / 2;
    while (m >= 2 && j >= m) {
      j -= m;
      m = m / 2;
    }
    j += m;
  }

  //Danielson-Lanzcos routine
  mmax=2;
  //external loop
  while (n > mmax) {
    istep = mmax << 1;
    theta = isForwardFFT * (2 * PI / mmax);
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    //internal loops
    for (m = 1; m < mmax; m += 2) {
      for (i = m; i <= n; i += istep) {
          j = i + mmax;
          tempr = wr * data[j - 1] - wi * data[j];
          tempi = wr * data[j] + wi * data[j - 1];
          data[j - 1] = data[i - 1] - tempr;
          data[j] = data[i] - tempi;
          data[i - 1] += tempr;
          data[i] += tempi;
      }
      wr = (wtemp = wr) * wpr - wi * wpi + wr;
      wi = wi * wpr + wtemp * wpi + wi;
    }
    mmax = istep;
  }
  if (isForwardFFT == -1) {
    for (int x = 0; x < n; ++x) {
      data.at(x) /= number_of_complex_samples;
    }
  }
}

int main() {
  static_assert(sizeof(wav_hdr) == 44, "");

  int i = sizeof(rawData) / sizeof(int16_t);
  int fsize = 2 * i;
  cout << "Lines in the input file: " << i << '\n';
  printf("file size: %u\n", fsize);

  // Pad the data with zeroes
  int nearestPowerOf2 = 1;
  long nearestPowerOf2Num = 1;
  while (nearestPowerOf2Num < i) {
    nearestPowerOf2Num = nearestPowerOf2Num << 1;
    ++nearestPowerOf2;
  }

  std::cout << nearestPowerOf2 << " " << nearestPowerOf2Num << std::endl;

  // copy data into vector
  for (int x = 0; x < i; ++x) {
    paddedData.push_back(static_cast<float>(rawData[x]));
  }
  std::cout << "Size :: " << paddedData.size() << std::endl;

  std::cout << (nearestPowerOf2Num - i) << std::endl;

  for (int x = 0; x < (nearestPowerOf2Num - i); ++x) {
    paddedData.push_back(0.0);
  }

  // for (int x = 0; x < paddedData.size(); ++x) {
  //   std::cout << "data : " << x << " " << paddedData.at(x) << std::endl;
  // }
  std::cout << "Size :: " << paddedData.size() << std::endl;

  // Adding imaginary components
  std::cout << "added alternate imaginary values : " << paddedData.size() << endl;
  int paddedDataSize = paddedData.size() * 2;
  for (int x = 0; x < paddedDataSize; x += 2) {
    paddedData.insert(paddedData.begin() + x + 1, 0.0);
    // std::cout << x << std::endl;
  }
  std::cout << "added alternate imaginary values : " << paddedData.size() << endl;

  generatePlot("inputSignal.dat");

  FFT(paddedData, i, 1);

  generatePlot("fftSignal.dat");

  for (int x = 0; x < paddedData.size(); ++x) {
    cout << paddedData.at(x) << ", ";
    if (x % 2 != 0) {
      cout << endl;
    }
  }

  FFT(paddedData, i, -1);

  generatePlot("reverseFftSignal.dat");

  std::cout << "Took inverse of the data : \n";
    for (int x = 0; x < paddedData.size(); ++x) {
    cout << paddedData.at(x) << ", ";
    if (x % 2 != 0) {
      cout << endl;
    }
  }

  int16_t processedData[FILE_SIZE];
  // std::memcpy(processedData, rawData, FILE_SIZE);
  ////////////////////////////////////////////////////////////////////
  // DO YOUR PROCESSING HERE, ON "processedData"
  ////////////////////////////////////////////////////////////////////

  // CONVERTS "processedData" TO AN AUDIO FILE
  // wav_hdr wav;
  // wav.ChunkSize = fsize + sizeof(wav_hdr) - 8;
  // wav.Subchunk2Size = fsize + sizeof(wav_hdr) - 44;

  // FILE-NAME: "alice.wav", SAVE LOCATION: same as this cpp-file
  // std::ofstream out("alice.wav", std::ios::binary);
  // out.write(reinterpret_cast<const char *>(&wav), sizeof(wav));
  // out.write(reinterpret_cast<char *>(processedData), fsize);

  // out.close();

  return 0;
}