# Fourier eavesdropping.
Eve is listening though a wall to a conversation coming from the neighboring room between
Alice and Bob. Since not all frequencies are transmitted equally well through the wall, the voices
sound muted, and Eve cannot understand a word of what is being said. The file "alice.h"
contains the waveform of one such muted sentence, spoken by Alice. Try to reconstruct the
original signal, to help Eve in her eavesdropping.

1> Write your own program to compute the Fast Fourier Transform (FFT) - due to the
large number of samples the DFT would take too long. You may use and adapt the code
provided in Numerical Recipes.

2> Test your program first on an artificially generated signal, consisting of only a few known
frequencies, and check if they appear in the Fourier spectrum at the expected positions.

3> The sentence spoken by Alice was recorded at the rate of 20000 samples per second. The
supplied data will require padding with zeroes, to the closest power of two.

4> You may assume that the low frequencies have been transmitted better that the high
frequencies. Use the following transmission function:
  T(v) = (1 - delta) * e^(-v^2/sigma^2) + delta

5> Fourier transform the signal, reconstruct it in the Fourier space by dividing by T, and
transform it back to the time space.

6> Use the provided code, or an alternative, to convert your data to an audio signal. What
was the sentence spoken by Alice?

## File Details
* alice.h: contains input data for audio signal
* alice input.wav: input audio file
* main.cpp: contains logic to create audio file from given data and main calling function for processing the data to generate the same
* fft_utils: header and cpp file for processing of the given data (also includes code to generate plot file)
* fft_test: header and cpp file for task 2 (generating a sine wave and testing fft logic)
* test_eight_bit.h: input 8 bit data to test the fft code
* alice.wav : Output audio file

On running the main.exe after compiling with below command(Please update project location on local system), many *.dat files would be generated in the format for 2d plotting. Any plotting tool can be used to plot the same (We have used gnuplot for this purpose.)
```
C:\MinGW\bin\g++.exe -std=c++17 -g C:\Users\lenovo\Desktop\audio_denoise\*.cpp C:\Users\lenovo\Desktop\audio_denoise\*.h -o C:\Users\lenovo\Desktop\audio_denoise\main.exe
```

## Images
* audio-denoise.png: flow chart for the evesdropping code
* input.png: Plot for input audio signal on time scale
* fft.png: Frequency domain plot for given audio signal after taking fourier transform
* tvfunction.png: Transmission function plot
* final-fft.png: FFT plot after applying the transmission function i.e. FFT / T(v)
* ifft.png: Denoised audio signal plot after taking inverse fft
* sineInput.png: Task 2 signal generated using code in fft_test.cpp
* sineFft.png: Frequency plot after taking the fft of automatically generated sine wave
* sineIfft.png: Plot after taking inverse fft of the generate signal