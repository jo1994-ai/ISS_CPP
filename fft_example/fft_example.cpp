// STFT / iSTFT with FFTW, using the same analysis and synthesis approach as
// the original piva C++ code (include/piva/stft.hpp):
//
//   analysis window   periodic Hann            (piva::windows::hann)
//   synthesis window  dual of the Hann window  (piva::windows::make_dual)
//   stft              piva::impl::stft,  output (n_freq, n_chan, n_frames)
//   istft             piva::impl::istft, weighted overlap-add
//
// Framing, zero-padding, frame count, output length and scaling (inverse FFT
// divided by n_fft) all follow the original, so the results are identical.
// The only differences: FFTW is called directly instead of through
// xtensor-fftw, the FFTW plans are made once instead of for every frame, and
// the loops are single-threaded instead of TBB.
//
// The demo reads a WAV file (by default the same 2-mic speech mixture used
// for the online IVA algorithms), goes to the STFT domain and back, and
// writes the reconstructed audio.
//
// Usage:
//   fft_example [input.wav=audio/mixture.wav] [output.wav=output/reconstructed.wav]
//               [n_analysis=1024] [n_shift=n_analysis/2]
//
// Dependencies: C++17 compiler + FFTW3.
//
// Build and run:  ./build.sh   (or: make run)
#define _USE_MATH_DEFINES  // M_PI on MSVC
#include <fftw3.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

#include "wav_io.hpp"

using cplx = std::complex<double>;

// Periodic Hann window, same as piva::windows::hann
std::vector<double> hann(long length) {
  std::vector<double> window(length);
  double length_inv = 1. / length;
  for (long n = 0; n < length; n++)
    window[n] = 0.5 * (1. - std::cos(2 * M_PI * length_inv * n));
  return window;
}

// Synthesis window for perfect reconstruction with the given analysis
// window and shift, same as piva::windows::make_dual
std::vector<double> make_dual(const std::vector<double>& awindow, long n_shift) {
  long length = long(awindow.size());
  std::vector<double> norm(length, 0.0);

  long n = 0;

  // move the window back as far as possible while still overlapping
  while (n - n_shift > -length) n -= n_shift;

  while (n < length) {
    if (n == 0)
      for (long m = 0; m < length; m++) norm[m] += awindow[m] * awindow[m];
    else if (n < 0)
      for (long m = 0; m < n + length; m++) norm[m] += awindow[m - n] * awindow[m - n];
    else
      for (long m = n; m < length; m++) norm[m] += awindow[m - n] * awindow[m - n];
    n += n_shift;
  }

  std::vector<double> swindow(length);
  for (long m = 0; m < length; m++) swindow[m] = awindow[m] / norm[m];
  return swindow;
}

// STFT tensor (n_freq, n_chan, n_frames), stored row-major
struct Stft {
  long n_freq = 0, n_chan = 0, n_frames = 0;
  std::vector<cplx> data;
  cplx& operator()(long f, long c, long t) { return data[(f * n_chan + c) * n_frames + t]; }
};

// Same as piva::impl::stft.
// td_signals: (n_samples, n_chan), row-major
Stft stft(const std::vector<double>& td_signals, long n_samples, long n_chan,
          const std::vector<double>& window, long n_shift, long n_zeropad_front,
          long n_zeropad_back) {
  long n_analysis = long(window.size());

  long n_frames = (n_samples + n_analysis - n_shift) / n_shift;
  if (n_frames * n_shift - n_analysis + n_shift < n_samples) n_frames++;

  long n_fft = n_zeropad_front + n_analysis + n_zeropad_back;
  long n_freq = n_fft / 2 + 1;

  Stft out;
  out.n_freq = n_freq;
  out.n_chan = n_chan;
  out.n_frames = n_frames;
  out.data.assign(n_freq * n_chan * n_frames, 0.0);

  double* input = fftw_alloc_real(n_fft);
  fftw_complex* spectrum = fftw_alloc_complex(n_freq);
  fftw_plan plan = fftw_plan_dft_r2c_1d(int(n_fft), input, spectrum, FFTW_ESTIMATE);

  for (long frame = 0; frame < n_frames; frame++) {
    long n = frame * n_shift;   // current offset in the input buffer
    long s2 = n + n_shift;      // end of current analysis window
    long s1 = s2 - n_analysis;  // beginning of analysis window

    long i1 = n_zeropad_front;
    long i2 = i1 + n_analysis;

    // handle the beginning
    if (s1 < 0) {
      i1 -= s1;
      s1 = 0;
    }

    // handle the end of the buffer
    if (s2 > n_samples) {
      i2 = i1 + n_samples - s1;
      s2 = n_samples;
    }

    for (long chan = 0; chan < n_chan; chan++) {
      // zeros outside the signal, samples inside, then the analysis window
      std::fill(input, input + n_fft, 0.0);
      for (long i = i1, s = s1; i < i2; i++, s++) input[i] = td_signals[s * n_chan + chan];
      for (long i = 0; i < n_analysis; i++) input[n_zeropad_front + i] *= window[i];

      fftw_execute(plan);
      for (long f = 0; f < n_freq; f++)
        out(f, chan, frame) = cplx(spectrum[f][0], spectrum[f][1]);
    }
  }

  fftw_destroy_plan(plan);
  fftw_free(input);
  fftw_free(spectrum);
  return out;
}

// Same as piva::impl::istft.
// Returns (n_samples, n_chan), row-major, with
// n_samples = n_frames * n_shift - (n_analysis - n_shift).
std::vector<double> istft(Stft& stft_signals, const std::vector<double>& window,
                          long n_shift, long n_zeropad_front, long n_zeropad_back,
                          long& n_samples_out) {
  long n_freq = stft_signals.n_freq;
  long n_chan = stft_signals.n_chan;
  long n_frames = stft_signals.n_frames;
  long n_analysis = long(window.size());
  long n_fft = n_zeropad_front + n_analysis + n_zeropad_back;

  if (n_fft / 2 + 1 != n_freq) {
    std::fprintf(stderr, "istft: window and zero-padding do not match the STFT size\n");
    n_samples_out = 0;
    return {};
  }

  long n_samples = n_frames * n_shift + (n_shift - n_analysis);
  n_samples_out = n_samples;
  std::vector<double> td_signals(n_samples * n_chan, 0.0);

  fftw_complex* spectrum = fftw_alloc_complex(n_freq);
  double* output = fftw_alloc_real(n_fft);
  fftw_plan plan = fftw_plan_dft_c2r_1d(int(n_fft), spectrum, output, FFTW_ESTIMATE);
  double inv_n = 1.0 / double(n_fft);  // FFTW's inverse is unnormalized

  // do the overlap add thing
  for (long frame = 0; frame < n_frames; frame++) {
    long n = frame * n_shift;                 // current offset in the output buffer
    long s2 = n + n_shift + n_zeropad_back;  // end of current fft
    long s1 = s2 - n_fft;                     // beginning of fft

    long i1 = 0;  // offset in the current frame buffer
    long i2 = n_fft;

    // handle the beginning
    if (s1 < 0) {
      i1 -= s1;
      s1 = 0;
    }

    // handle the end of the buffer
    if (s2 > n_samples) {
      i2 = i1 + n_samples - s1;
      s2 = n_samples;
    }

    for (long chan = 0; chan < n_chan; chan++) {
      // inverse fft (a c2r plan may overwrite its input, so refill every time)
      for (long f = 0; f < n_freq; f++) {
        spectrum[f][0] = stft_signals(f, chan, frame).real();
        spectrum[f][1] = stft_signals(f, chan, frame).imag();
      }
      fftw_execute(plan);
      for (long i = 0; i < n_fft; i++) output[i] *= inv_n;

      // apply synthesis window
      for (long i = 0; i < n_analysis; i++) output[n_zeropad_front + i] *= window[i];

      // overlap and add
      for (long i = i1, s = s1; i < i2; i++, s++) td_signals[s * n_chan + chan] += output[i];
    }
  }

  fftw_destroy_plan(plan);
  fftw_free(spectrum);
  fftw_free(output);
  return td_signals;
}

int main(int argc, char** argv) {
  std::string input_filename = argc > 1 ? argv[1] : "audio/mixture.wav";
  std::string output_filename = argc > 2 ? argv[2] : "output/reconstructed.wav";

  // STFT parameters, as in examples-cpp/src/auxiva.cpp
  const long n_analysis = argc > 3 ? std::stol(argv[3]) : 1024;
  const long n_shift = argc > 4 ? std::stol(argv[4]) : n_analysis / 2;
  const long n_zeropad_front = 0, n_zeropad_back = 0;

  // read the speech signal
  piva_online::WavData wav;
  try {
    wav = piva_online::read_wav(input_filename);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
  const long n_samples = long(wav.n_samples), n_chan = long(wav.n_chan);
  const std::vector<double>& x = wav.samples;
  std::printf("Input %s: %ld channels, %ld samples at %d Hz\n", input_filename.c_str(), n_chan,
              n_samples, wav.samplerate);

  // create the windows for the STFT
  auto analysis_window = hann(n_analysis);
  auto synthesis_window = make_dual(analysis_window, n_shift);

  // go to time-frequency domain
  Stft X = stft(x, n_samples, n_chan, analysis_window, n_shift, n_zeropad_front, n_zeropad_back);
  std::printf("STFT shape (n_freq, n_chan, n_frames): (%ld, %ld, %ld)\n", X.n_freq, X.n_chan,
              X.n_frames);

  // ... process X here (e.g. separation) ...

  // go back to time domain
  long n_out = 0;
  auto y = istft(X, synthesis_window, n_shift, n_zeropad_front, n_zeropad_back, n_out);
  std::printf("iSTFT output: %ld samples x %ld channels (input: %ld samples)\n", n_out, n_chan,
              n_samples);

  // output sample s lines up with input sample s
  double max_err = 0.0;
  for (long s = 0; s < std::min(n_samples, n_out); s++)
    for (long c = 0; c < n_chan; c++)
      max_err = std::max(max_err, std::abs(y[s * n_chan + c] - x[s * n_chan + c]));
  std::printf("Reconstruction max error: %.2e\n", max_err);

  // write the reconstruction, trimmed to the input length, at the original
  // level (no normalization) so it can be compared sample by sample
  y.resize(n_samples * n_chan);
  try {
    piva_online::write_wav(output_filename, y, size_t(n_chan), wav.samplerate);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s (does the output folder exist?)\n", e.what());
    return 1;
  }
  std::printf("Wrote %s\n", output_filename.c_str());
}
