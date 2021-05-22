#include "fft_utils.h"

const double PI = 2 * acos(0.0);

void swap (std::vector<float>& data, unsigned long first, unsigned long second) {
  double temp = data.at(first);
  data.at(first) = data.at(second);
  data.at(second) = temp;
}

// Adding imaginary components
void addImaginaryComponents(std::vector<float>& data) {
  int dataSize = data.size() * 2;
  for (int x = 0; x < dataSize; x += 2) {
    data.insert(data.begin() + x + 1, 0.0);
  }
  std::cout << "Added alternate imaginary values : " << data.size() << std::endl;
}

// Pad the data with zeroes
void padZeroTillPowerOfTwo(std::vector<float>& data) {
  int nearestPowerOf2 = 1;
  long nearestPowerOf2Num = 1;
  while (nearestPowerOf2Num < data.size()) {
    nearestPowerOf2Num = nearestPowerOf2Num << 1;
    ++nearestPowerOf2;
  }

  std::cout << "Nearest Power Of Two : count : " << nearestPowerOf2 << " num : " << nearestPowerOf2Num << std::endl;

  int toBeFilled = nearestPowerOf2Num - data.size();
  
  std::cout << "To be filled places :: " << toBeFilled << std::endl;

  for (int x = 0; x < toBeFilled; ++x) {
    data.push_back(0.0);
  }
  std::cout << "Size after padding with 0 :: " << data.size() << std::endl;
}

void erasePaddedZeroes(std::vector<float>& data, int originalSize) {
  std::cout << "Erasing the padded Zeroes." << std::endl;
  std::cout << "Original size with imaginary : " << originalSize * 2 << std::endl;
  std::cout << "Padded size with imaginary : " << data.size() << std::endl;
  data.erase(data.begin() + originalSize * 2, data.begin() + data.size());

  std::cout << "FINAL SIZE : " << data.size() << std::endl;
}

// T(ν) = (1 − δ) exp(−ν2/σ2) + δ
void applyTransmissionEqn(std::vector<float>& data) {
  std::ofstream tvPlot { "tvPlot.dat" };
  float sigma = 1010;
  float delta = .0001;
  long paddedSize = data.size();
  for (int x = 0; x < paddedSize; x += 2) {
    float nu1 = (x / 2);
    float nu2 = (x / 2) - (paddedSize / 2);
    float tv = ((1 - delta) * (
      std::exp(- std::pow(nu1, 2) / std::pow(sigma, 2)) +
      std::exp(- std::pow(nu2, 2) / std::pow(sigma, 2))
    )) + delta;
    tvPlot << static_cast<int>(x/2) << " " << tv << std::endl;
    data.at(x) = data.at(x) / tv;
    data.at(x + 1) = data.at(x + 1) / tv;
  }
  tvPlot.close();
}

// Generate dat file for gnu plot
void generatePlot(std::vector<float>& data, const std::string& fileName, bool isComplexPlot) {
  std::ofstream plot { fileName };
  for (int x = 0; x < data.size(); x += 2) {
    float val = !isComplexPlot ?
      data.at(x) :
      std::sqrt(std::pow(data.at(x), 2) + std::pow(data.at(x + 1), 2));
    plot << (x / 2) << " " << val << std::endl;
  }
  plot.close();
}

/*
  isForwardFFT = 1 for FFT and -1 for Inverse FFT
*/
void FFT (std::vector<float>& data, unsigned long number_of_complex_samples, int isForwardFFT) {
  //variables for trigonometric recurrences
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;
  
  n = number_of_complex_samples * 2;

  //real part of the complex is on the even-indexes
  //and the complex part is on the odd-indexes
  j = 0;
  for (i = 0; i < n; i += 2) {
    if (j > i) {
      //swap the real part
      swap(data, j, i);
      //swap the complex part
      swap(data, j + 1, i + 1);
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
          tempr = wr * data[j - 1] + wi * data[j];
          tempi = wr * data[j] - wi * data[j - 1];
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
