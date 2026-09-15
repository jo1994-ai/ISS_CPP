// 16-bit fixed-point (Q15) version of fft_example.
//
// Reads a 16-bit speech recording, computes the Q15 STFT and iSTFT with
// integer arithmetic (q15_fft.hpp), writes the reconstructed audio, and
// reports how close it is to the input and to an exact double-precision FFT.
//
// Usage:
//   fft_example_q15 [input.wav=audio/mixture.wav]
//                   [output.wav=output/reconstructed_q15.wav]
//                   [n_fft=1024] [n_shift=n_fft/2]
//
// Dependencies: a C++17 compiler only.
#define _USE_MATH_DEFINES  // M_PI on MSVC
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "q15_fft.hpp"
#include "wav16.hpp"

// ---------------------------------------------------------------------------
// Double-precision reference, used only to measure the Q15 accuracy.
// Returns numpy.fft.rfft(x) / N, the same scaling as q15::FFT::forward.
static std::vector<std::complex<double>> reference_rfft(const std::vector<double>& x) {
  size_t n = x.size();
  std::vector<std::complex<double>> a(x.begin(), x.end());
  for (size_t i = 1, j = 0; i < n; i++) {
    size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(a[i], a[j]);
  }
  for (size_t len = 2; len <= n; len <<= 1) {
    std::complex<double> wlen = std::polar(1.0, -2.0 * M_PI / double(len));
    for (size_t i = 0; i < n; i += len) {
      std::complex<double> w = 1.0;
      for (size_t j = 0; j < len / 2; j++, w *= wlen) {
        auto u = a[i + j], v = a[i + j + len / 2] * w;
        a[i + j] = u + v;
        a[i + j + len / 2] = u - v;
      }
    }
  }
  a.resize(n / 2 + 1);
  for (auto& v : a) v /= double(n);
  return a;
}

static double snr_db(double signal_energy, double error_energy) {
  if (error_energy == 0.0) return INFINITY;
  return 10.0 * std::log10(signal_energy / error_energy);
}
// ---------------------------------------------------------------------------

int main(int argc, char** argv) {
  std::string input_filename = argc > 1 ? argv[1] : "audio/mixture.wav";
  std::string output_filename = argc > 2 ? argv[2] : "output/reconstructed_q15.wav";
  long n_fft = argc > 3 ? std::stol(argv[3]) : 1024;
  long n_shift = argc > 4 ? std::stol(argv[4]) : n_fft / 2;

  wav16::Wav wav;
  try {
    wav = wav16::read(input_filename);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
  const long n_samples = wav.n_samples, n_chan = wav.n_chan;
  const std::vector<int16_t>& x = wav.samples;
  std::printf("Input %s: %ld channels, %ld samples at %d Hz (16-bit)\n", input_filename.c_str(),
              n_chan, n_samples, wav.samplerate);

  // ---- Q15 processing (integer arithmetic) ----
  std::vector<int16_t> analysis_window;
  q15::SynthesisWindow synthesis_window;
  q15::Stft X;
  std::vector<int16_t> y;
  long n_out = 0;
  try {
    analysis_window = q15::hann(size_t(n_fft));
    synthesis_window = q15::make_dual(analysis_window, size_t(n_shift));
    X = q15::stft(x, n_samples, n_chan, analysis_window, n_shift);
    // ... process X here ...
    y = q15::istft(X, synthesis_window, n_shift, n_out);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }

  std::printf("Analysis window: Q15, synthesis window: Q%d\n", synthesis_window.frac_bits);
  std::printf("STFT shape (n_freq, n_chan, n_frames): (%ld, %ld, %ld), stored as int16\n",
              X.n_freq, X.n_chan, X.n_frames);
  std::printf("iSTFT output: %ld samples x %ld channels (input: %ld samples)\n", n_out, n_chan,
              n_samples);

  y.resize(size_t(n_samples * n_chan));
  try {
    auto dir = std::filesystem::path(output_filename).parent_path();
    if (!dir.empty()) std::filesystem::create_directories(dir);
    wav16::write(output_filename, y, n_chan, wav.samplerate);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
  std::printf("Wrote %s\n\n", output_filename.c_str());

  // ---- accuracy (double precision, for reporting only) ----
  double sig = 0.0, err = 0.0;
  long max_err = 0, n_diff = 0;
  for (size_t i = 0; i < y.size(); i++) {
    long d = long(y[i]) - long(x[i]);
    sig += double(x[i]) * double(x[i]);
    err += double(d) * double(d);
    max_err = std::max(max_err, std::labs(d));
    n_diff += d != 0;
  }
  std::printf("Reconstruction vs input:\n");
  std::printf("  SNR                  %.1f dB\n", snr_db(sig, err));
  std::printf("  max error            %ld LSB (of 32768)\n", max_err);
  std::printf("  samples changed      %ld of %zu\n", n_diff, y.size());

  // Q15 spectrum vs exact spectrum of the same frames
  double spec_sig = 0.0, spec_err = 0.0;
  std::vector<double> frame(n_fft);
  for (long t = 0; t < X.n_frames; t++) {
    long s1 = t * n_shift + n_shift - n_fft;
    for (long c = 0; c < n_chan; c++) {
      for (long i = 0; i < n_fft; i++) {
        long s = s1 + i;
        double sample = (s >= 0 && s < n_samples) ? x[size_t(s * n_chan + c)] / 32768.0 : 0.0;
        frame[i] = sample * 0.5 * (1.0 - std::cos(2.0 * M_PI * double(i) / double(n_fft)));
      }
      auto ref = reference_rfft(frame);
      for (long f = 0; f < X.n_freq; f++) {
        std::complex<double> q(X.re[X.idx(f, c, t)] / 32768.0, X.im[X.idx(f, c, t)] / 32768.0);
        spec_sig += std::norm(ref[f]);
        spec_err += std::norm(q - ref[f]);
      }
    }
  }
  std::printf("Q15 spectrum vs exact double FFT:\n");
  std::printf("  SNR                  %.1f dB\n", snr_db(spec_sig, spec_err));
  std::printf("(the same STFT in double precision reconstructs the input exactly)\n");
}
