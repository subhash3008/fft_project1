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

