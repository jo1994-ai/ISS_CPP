fft_example_q15: 16-bit fixed-point STFT
========================================

Q15 fixed-point version of `fft_example`: the same STFT / iSTFT approach as
the original piva C++ code (periodic Hann analysis window, dual synthesis
window, same framing and overlap-add), but with 16-bit integer data and
integer arithmetic, the way it is done on fixed-point DSP chips.

No FFTW, no libraries: only a C++17 compiler.

Q15 in short
------------

An `int16_t` value `v` represents `v / 32768`, i.e. numbers in [-1, 1).

| Data                        | Stored as | Format        |
|-----------------------------|-----------|---------------|
| audio samples               | int16     | Q15           |
| analysis window (Hann)      | int16     | Q15           |
| synthesis window (dual)     | int16     | Q14 (peak > 1)|
| cos / sin tables            | int16     | Q15           |
| STFT spectrum (re, im)      | int16     | Q15, scaled 1/N |
| products, sums, overlap-add | int32 / int64 intermediates | |

Floating point is used only to create the tables at startup and to print
the accuracy report; the STFT / iSTFT themselves are integer-only.

Folder contents
---------------

    q15_fft.hpp            Q15 FFT, windows, STFT, iSTFT
    wav16.hpp              16-bit WAV reading / writing as int16
    fft_example_q15.cpp    demo + accuracy report
    audio/mixture.wav      2-channel speech mixture, 16 kHz, 16-bit
    build.sh               macOS / Linux
    build.bat, CMakeLists.txt   Windows (Visual Studio)

Build and run
-------------

macOS / Linux:

    cd ~/Desktop/fft_example_q15
    ./build.sh

Windows: double-click `build.bat` (needs Visual Studio with
"Desktop development with C++").

Output:

    Input audio/mixture.wav: 2 channels, 72764 samples at 16000 Hz (16-bit)
    Analysis window: Q15, synthesis window: Q14
    STFT shape (n_freq, n_chan, n_frames): (513, 2, 144), stored as int16
    iSTFT output: 73216 samples x 2 channels (input: 72764 samples)
    Wrote output/reconstructed_q15.wav

    Reconstruction vs input:
      SNR                  45.7 dB
      max error            113 LSB (of 32768)
      samples changed      142942 of 145528
    Q15 spectrum vs exact double FFT:
      SNR                  38.0 dB

Other settings:

    ./fft_example_q15 <input.wav> <output.wav> [n_fft=1024] [n_shift=n_fft/2]

`n_fft` must be a power of two, input must be 16-bit PCM WAV.

Accuracy
--------

The spectrum is scaled by 1/N so that 16-bit bins cannot overflow. Every
FFT stage halves the values, so each doubling of N loses one more bit and
costs about 3 dB:

| n_fft | Reconstruction SNR | Spectrum SNR |
|-------|--------------------|--------------|
| 64    | 57.9 dB            | 51.1 dB      |
| 256   | 51.7 dB            | 44.2 dB      |
| 512   | 48.7 dB            | 41.1 dB      |
| 1024  | 45.7 dB            | 38.0 dB      |
| 4096  | 39.6 dB            | 31.7 dB      |

For comparison: double precision reconstructs the input exactly (~300 dB),
and 16-bit audio itself has ~96 dB dynamic range. The Q15 noise floor is
about -46 dB below the signal at N = 1024: audible in quiet passages with
headphones, fine for many detection / feature tasks.

Ways to get more accuracy on fixed-point hardware:

* smaller frames or more overlap (75 % overlap gives +3 dB)
* keep spectra in 32 bits (Q31) instead of 16
* scale by 1/sqrt(N) instead of 1/N, using the headroom the signal actually has


How it works:

Q15 format: an int16 value v stands for v/32768.
Everything stored is int16: samples, windows, sine/cosine tables and the spectra.
Wider integers only for intermediate math: products and sums use int32/int64 and are shifted back, as fixed-point DSP libraries such as CMSIS-DSP do.
Floating point: used only to build the tables at startup and for the accuracy report. The STFT and iSTFT themselves are integer-only.
Same approach as piva: periodic Hann analysis window, dual synthesis window, same framing and overlap-add. The synthesis window peaks above 1, so it's stored in Q14 (range ±2).
The FFT: a radix-2 fixed-point FFT written from scratch. It halves the values at every stage so 16-bit spectra can't overflow, giving a forward result scaled by 1/N. The inverse applies no scaling.
Is the loss a bug or just 16-bit limits? It's the limits: SNR drops about 3 dB every time the frame size doubles, the expected pattern for rounding noise.

n_fft	64	256	512	1024	4096
Reconstruction SNR	57.9 dB	51.7 dB	48.7 dB	45.7 dB	39.6 dB
Each extra FFT stage throws away one more bit. For comparison, double precision reconstructs the input exactly, and 16-bit audio itself has about 96 dB of range.

In practice:

Sound: noise about 46 dB below the speech. It's audible in quiet parts with headphones, so compare output/reconstructed_q15.wav with audio/mixture.wav.
Improving it:
use smaller frames, or 75% overlap (+3 dB),
keep spectra in 32 bits (Q31),
scale by 1/√N instead of 1/N.
Online IVA in Q15: far harder than the FFT, because of the covariance updates, divisions and square roots. Most real products use float32 instead.
Tested: built from a clean shell with build.sh and also through CMake (the Windows build path), with no warnings. The program uses only the system C and C++ libraries. build.bat itself hasn't been run on Windows.




Same in both
The formula: both compute the discrete Fourier transform $$X[k] = \sum_{n=0}^{N-1} x[n], e^{-j2\pi kn/N}$$
Frequency bins: N/2 + 1 bins, from 0 Hz to fs/2.
Window: periodic Hann for analysis, its dual for synthesis.
Framing, hop and overlap-add: the same as the original piva code.
Inverse sign: e^{+j…}.
Result: a full round trip gives the signal back.
What's different
Double version (FFTW)	Q15 version
Numbers	64-bit floating point	16-bit integers
FFT algorithm	FFTW's optimized routines	simple radix-2 FFT
Where the 1/N goes	on the inverse	on the forward
Rounding	none that matters (~1e-16)	at every multiply and every halving
Result	exact	exact spectrum plus rounding noise
The two spectra relate like this:
$$X_{\text{double}}[k] \approx N \cdot \frac{X_{\text{Q15}}[k]}{32768}$$
That's the same spectrum with a different scale factor, plus noise.

Why the 1/N moved to the forward transform
A 16-bit value can't go above 32767. Without scaling, the forward FFT adds up N samples, so a bin can reach N × 32767 (about 33 million for N = 1024), which overflows immediately. Halving at each of the log₂N stages keeps every value within 16 bits. In double precision, overflow isn't a concern, so FFTW can leave the forward transform unscaled.

Where the 45.7 dB loss comes from
The operation is unchanged; the loss comes from rounding:

Windowing: sample × window is rounded back to 16 bits.
Every FFT stage: the rotation multiply is rounded, and the halving drops the lowest bit. With N = 1024 that's 10 stages, so 10 bits of precision are lost.
Storage: the spectrum is stored in 16 bits.
The measurement shows this. The Q15 spectrum matches the exact double-precision FFT with 38 dB SNR: same values, plus rounding noise. With no rounding, as in the double version, the output is identical to the input.

In short, the FFT is the same mathematical operation in both versions. The Q15 version does it with 16-bit integers, so it scales earlier to avoid overflow and rounds at each step, which adds noise.

is Q15 different library?
No, Q15 isn't a library. It's a number format: a convention for storing fractions in an int16_t.

What Q15 means
Q: fixed-point format.
15: 15 bits after the binary point, plus 1 sign bit, for 16 bits total.
Value: real value = integer / 2^15 = integer / 32768.
int16 stored	Q15 value
32767	0.99997
16384	0.5
0	0.0
-16384	-0.5
-32768	-1.0
Multiplying two Q15 numbers gives a Q30 result, so you shift right by 15 to get back to Q15. That shift is the round_shift(a * b, 15) you see all through the code. Other Q formats follow the same idea: Q14 has range ±2 (used for your synthesis window), and Q31 is a 32-bit version with more precision.

In your fft_example_q15
No library is used. I wrote the Q15 FFT by hand in q15_fft.hpp. The only things the program needs are the C++ compiler and its standard library.

Folder	FFT comes from
fft_example, fft_example_windows	FFTW library (double)
fft_example_q15	your own code (q15_fft.hpp), no library
Libraries that provide Q15 FFTs
FFTW supports only floating point. For fixed-point hardware there are dedicated libraries:

Library	Q15 FFT functions	Used on	License
CMSIS-DSP (Arm)	arm_cfft_q15, arm_rfft_q15	ARM Cortex-M/A microcontrollers	Apache 2.0
TI DSPLIB	DSP_fft16x16	Texas Instruments C6000/C5000 DSPs	TI license
KissFFT	compiled with FIXED_POINT=16	any CPU, portable C	BSD
Cadence NatureDSP	fft_cplx16x16	Tensilica HiFi audio DSPs	commercial
ADI libraries	vendor-specific	Analog Devices SHARC/Blackfin	vendor
CMSIS-DSP: the most common choice on ARM microcontrollers. Its arm_cfft_q15 also scales by 1/N, like your code, and uses the same approach of int16 storage with 32-bit intermediate math.
KissFFT: the easiest to try on a PC. It's small, portable C, and becomes a Q15 FFT when compiled with FIXED_POINT=16.
The chip vendor's library: use this when targeting that DSP chip, because it's hand-optimized in assembly.
If you want to compare against a library, I can add a version of fft_example_q15 that uses KissFFT or CMSIS-DSP next to the hand-written FFT and compares their outputs.

so there is dofference:	Double version (FFTW)	Q15 version
Numbers	64-bit floating point	16-bit integers
FFT algorithm	FFTW's optimized routines	simple radix-2 FFT
Where the 1/N goes	on the inverse	on the forward
Rounding	none that matters (~1e-16)	at every multiply and every halving
Result	exact	exact spectrum plus rounding noise
Yes, those differences are real, but they aren't equally important. Only one of them changes the sound: the number format.

Difference	Changes the output?	Why
Numbers: 64-bit float vs 16-bit int	Yes, the only real cause	16 bits force rounding, which adds noise (45.7 dB SNR)
Algorithm: FFTW vs radix-2	No	Both compute the same DFT. Run in double precision, a radix-2 FFT gives the same numbers as FFTW (to ~1e-16); FFTW is just faster
Where 1/N goes	No	It only moves a constant factor. Forward × inverse is 1/N either way, so the round trip matches
Rounding	Yes	This is a consequence of the 16-bit format, not a separate choice
Result		exact vs. exact + noise, which follows from the rows above
In the Q15 program, the "exact double FFT" used for the accuracy report is the same radix-2 algorithm in double precision, and it matches the exact transform. So the algorithm isn't the problem.

Only the number format was really chosen; the rest follows from it:


16-bit integers
   ├─→ values must stay ≤ 32767  →  scale by 1/N on the forward transform
   ├─→ products must be shifted back  →  rounding at every step  →  noise
   └─→ FFTW has no integer version  →  own radix-2 FFT (or CMSIS-DSP, KissFFT)
Put simply, both versions compute the same FFT. The Q15 version uses 16-bit integers, which forces it to scale earlier, round at each step and use its own FFT code. The rounding is the only thing you can hear.