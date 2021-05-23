// This code transforms the data stored in "alice.h" into a WAV-file
// and saves it as "alice.wav".

#include "alice.h"
// #include "test_eight_bit.h"
#include "fft_utils.h"
#include "fft_test.h"

using namespace std;

#define FILE_SIZE 100000 

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

int main() {

  // // TASK 2
  // std::vector<float> testData;
  // runFftTest(testData);
  // return 0;
  // // END TASK 2

  static_assert(sizeof(wav_hdr) == 44, "");

  int i = sizeof(rawData) / sizeof(int16_t);
  int fsize = 2 * i;
  cout << "Lines in the input file: " << i << '\n';
  printf("file size: %u\n", fsize);

  ////////////////////////////////////////////////////////////////////
  // DO YOUR PROCESSING HERE
  ////////////////////////////////////////////////////////////////////

  // copy data into vector
  std::vector<float> paddedData;
  for (int x = 0; x < i; ++x) {
    paddedData.push_back(static_cast<float>(rawData[x]));
  }
  std::cout << "Initial Size :: " << paddedData.size() << std::endl;

  addImaginaryComponents(paddedData);
  
  generatePlot(paddedData, "inputSignal.dat");

  padZeroTillPowerOfTwo(paddedData);

  FFT(paddedData, paddedData.size() / 2 , 1);

  generatePlot(paddedData, "fftSignal.dat", true);

  applyTransmissionEqn(paddedData);

  generatePlot(paddedData, "modifiedFFT.dat", true);

  FFT(paddedData, paddedData.size() / 2, -1);

  std::cout << "Took inverse of the data." << std::endl;

  erasePaddedZeroes(paddedData, i);

  generatePlot(paddedData, "reverseFftSignal.dat");

  int16_t processedData[FILE_SIZE];
  // Populating processedData
  for (int x = 0; x < paddedData.size(); x += 2) {
    // processedData[x / 2] = static_cast<int16_t>(std::sqrt(
    //   std::pow(paddedData.at(x), 2) +
    //   std::pow(paddedData.at(x + 1), 2)
    // ));
    processedData[x / 2] = static_cast<int16_t>(paddedData.at(x));
  }

  /////////////////////////
  // PROCESSING ENDS     //
  /////////////////////////
  // int16_t processedData[FILE_SIZE];
  // std::memcpy(processedData, rawData, FILE_SIZE);

  // CONVERTS "processedData" TO AN AUDIO FILE
  std::cout << "Generating audio file\n.";
  wav_hdr wav;
  wav.ChunkSize = fsize + sizeof(wav_hdr) - 8;
  wav.Subchunk2Size = fsize + sizeof(wav_hdr) - 44;

  // FILE-NAME: "alice.wav", SAVE LOCATION: same as this cpp-file
  std::ofstream out("alice.wav", std::ios::binary);
  out.write(reinterpret_cast<const char *>(&wav), sizeof(wav));
  out.write(reinterpret_cast<char *>(processedData), fsize);

  out.close();

  return 0;
}
