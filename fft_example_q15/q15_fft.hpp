// 16-bit fixed-point (Q15) FFT, STFT and iSTFT.
//
// Q15: an int16_t value v represents v / 32768, i.e. the range [-1, 1).
//
// Data that is stored (samples, windows, cos/sin tables, spectra) is int16_t.
// Products and sums are computed in int32_t / int64_t and then shifted back,
// as fixed-point DSP libraries (e.g. CMSIS-DSP q15) do. The processing is
// integer-only; floating point is used only once, to create the tables.
//
// Scaling (so that 16-bit spectra cannot overflow):
//   forward  X_q15[k] = (1/N) sum_n x[n] e^{-j 2 pi k n / N}   (halved at every stage)
//   inverse  x[n]     =       sum_k X[k] e^{+j 2 pi k n / N}   (no scaling)
// so forward followed by inverse gives back x, like numpy rfft / irfft.
//
// STFT framing and windows follow the original piva C++ code
// (include/piva/stft.hpp): periodic Hann analysis window, its dual as the
// synthesis window, frame t covers samples [t*shift + shift - N, t*shift + shift).
// N must be a power of two.
#pragma once

#define _USE_MATH_DEFINES  // M_PI on MSVC
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace q15 {

inline int16_t saturate16(int64_t v) {
  if (v > 32767) return 32767;
  if (v < -32768) return -32768;
  return int16_t(v);
}

// round(v / 2^shift) for signed integers
inline int64_t round_shift(int64_t v, int shift) {
  return (v + (int64_t(1) << (shift - 1))) >> shift;
}

// double in [-1, 1] -> Q15, used only when building tables
inline int16_t from_double(double v, int frac_bits = 15) {
  return saturate16(std::llround(v * double(int64_t(1) << frac_bits)));
}

class FFT {
 public:
  explicit FFT(size_t n) : n_(n), cos_(n / 2), sin_(n / 2), rev_(n) {
    if (n < 2 || (n & (n - 1)) != 0) throw std::invalid_argument("FFT size must be a power of two");
    for (size_t k = 0; k < n / 2; k++) {
      cos_[k] = from_double(std::cos(2.0 * M_PI * double(k) / double(n)));
      sin_[k] = from_double(std::sin(2.0 * M_PI * double(k) / double(n)));
    }
    size_t bits = 0;
    while ((size_t(1) << bits) < n) bits++;
    for (size_t i = 0; i < n; i++) {
      size_t r = 0;
      for (size_t b = 0; b < bits; b++)
        if (i & (size_t(1) << b)) r |= size_t(1) << (bits - 1 - b);
      rev_[i] = r;
    }
    wr32_.resize(n);
    wi32_.resize(n);
    wr64_.resize(n);
    wi64_.resize(n);
  }

  size_t size() const { return n_; }

  // x: N real Q15 samples -> re, im: N/2 + 1 Q15 bins, scaled by 1/N
  void forward(const int16_t* x, int16_t* re, int16_t* im) {
    for (size_t i = 0; i < n_; i++) {
      wr32_[rev_[i]] = x[i];
      wi32_[rev_[i]] = 0;
    }
    for (size_t len = 2; len <= n_; len <<= 1) {
      size_t half = len / 2, step = n_ / len;
      for (size_t i = 0; i < n_; i += len) {
        for (size_t j = 0; j < half; j++) {
          int32_t c = cos_[j * step], s = sin_[j * step];
          int32_t xr = wr32_[i + j + half], xi = wi32_[i + j + half];
          // v = x * e^{-j theta}
          int32_t vr = int32_t(round_shift(int64_t(xr) * c + int64_t(xi) * s, 15));
          int32_t vi = int32_t(round_shift(int64_t(xi) * c - int64_t(xr) * s, 15));
          int32_t ur = wr32_[i + j], ui = wi32_[i + j];
          // halve at every stage: magnitudes stay within 16 bits
          wr32_[i + j] = (ur + vr + 1) >> 1;
          wi32_[i + j] = (ui + vi + 1) >> 1;
          wr32_[i + j + half] = (ur - vr + 1) >> 1;
          wi32_[i + j + half] = (ui - vi + 1) >> 1;
        }
      }
    }
    for (size_t k = 0; k <= n_ / 2; k++) {
      re[k] = saturate16(wr32_[k]);
      im[k] = saturate16(wi32_[k]);
    }
  }

  // re, im: N/2 + 1 Q15 bins (as produced by forward) -> x: N samples (int32,
  // Q15 scale; can exceed 16 bits before windowing and overlap-add)
  void inverse(const int16_t* re, const int16_t* im, int32_t* x) {
    // full spectrum of a real signal: X[N-k] = conj(X[k])
    for (size_t k = 0; k < n_; k++) {
      size_t src = k <= n_ / 2 ? k : n_ - k;
      int64_t r = re[src], i = k <= n_ / 2 ? im[src] : -im[src];
      wr64_[rev_[k]] = r;
      wi64_[rev_[k]] = i;
    }
    for (size_t len = 2; len <= n_; len <<= 1) {
      size_t half = len / 2, step = n_ / len;
      for (size_t i = 0; i < n_; i += len) {
        for (size_t j = 0; j < half; j++) {
          int64_t c = cos_[j * step], s = sin_[j * step];
          int64_t xr = wr64_[i + j + half], xi = wi64_[i + j + half];
          // v = x * e^{+j theta}
          int64_t vr = round_shift(xr * c - xi * s, 15);
          int64_t vi = round_shift(xi * c + xr * s, 15);
          int64_t ur = wr64_[i + j], ui = wi64_[i + j];
          wr64_[i + j] = ur + vr;
          wi64_[i + j] = ui + vi;
          wr64_[i + j + half] = ur - vr;
          wi64_[i + j + half] = ui - vi;
        }
      }
    }
    for (size_t i = 0; i < n_; i++) x[i] = int32_t(wr64_[i]);
  }

 private:
  size_t n_;
  std::vector<int16_t> cos_, sin_;
  std::vector<size_t> rev_;
  std::vector<int32_t> wr32_, wi32_;
  std::vector<int64_t> wr64_, wi64_;
};

// Periodic Hann window in Q15
inline std::vector<int16_t> hann(size_t length) {
  std::vector<int16_t> w(length);
  for (size_t n = 0; n < length; n++)
    w[n] = from_double(0.5 * (1.0 - std::cos(2.0 * M_PI * double(n) / double(length))));
  return w;
}

// Synthesis window (same construction as piva::windows::make_dual), computed
// from the quantized analysis window. Its peak can exceed 1, so it is stored
// with `frac_bits` fractional bits (Q15, Q14, ...) chosen to fit in int16.
struct SynthesisWindow {
  std::vector<int16_t> w;
  int frac_bits = 15;
};

inline SynthesisWindow make_dual(const std::vector<int16_t>& awindow, size_t n_shift) {
  long length = long(awindow.size());
  long shift = long(n_shift);
  std::vector<double> a(length), norm(length, 0.0), dual(length);
  for (long m = 0; m < length; m++) a[m] = awindow[m] / 32768.0;

  long n = 0;
  while (n - shift > -length) n -= shift;
  for (; n < length; n += shift)
    for (long m = std::max(0L, n); m < std::min(length, n + length); m++)
      norm[m] += a[m - n] * a[m - n];

  double peak = 0.0;
  for (long m = 0; m < length; m++) {
    dual[m] = norm[m] > 0.0 ? a[m] / norm[m] : 0.0;
    peak = std::max(peak, std::abs(dual[m]));
  }

  SynthesisWindow sw;
  while (sw.frac_bits > 0 && peak * double(1 << sw.frac_bits) > 32767.0) sw.frac_bits--;
  sw.w.resize(length);
  for (long m = 0; m < length; m++) sw.w[m] = from_double(dual[m], sw.frac_bits);
  return sw;
}

// Q15 STFT, shape (n_freq, n_chan, n_frames)
struct Stft {
  long n_freq = 0, n_chan = 0, n_frames = 0;
  std::vector<int16_t> re, im;
  size_t idx(long f, long c, long t) const { return size_t((f * n_chan + c) * n_frames + t); }
};

// td_signals: (n_samples, n_chan) Q15, row-major. Same framing as piva::impl::stft.
inline Stft stft(const std::vector<int16_t>& td_signals, long n_samples, long n_chan,
                 const std::vector<int16_t>& window, long n_shift) {
  long n_fft = long(window.size());
  FFT fft{size_t(n_fft)};

  long n_frames = (n_samples + n_fft - n_shift) / n_shift;
  if (n_frames * n_shift - n_fft + n_shift < n_samples) n_frames++;

  Stft out;
  out.n_freq = n_fft / 2 + 1;
  out.n_chan = n_chan;
  out.n_frames = n_frames;
  out.re.assign(size_t(out.n_freq * n_chan * n_frames), 0);
  out.im.assign(out.re.size(), 0);

  std::vector<int16_t> frame(n_fft), re(out.n_freq), im(out.n_freq);

  for (long t = 0; t < n_frames; t++) {
    long s2 = t * n_shift + n_shift;  // end of the analysis window
    long s1 = s2 - n_fft;             // beginning (negative for the first frames)
    for (long c = 0; c < n_chan; c++) {
      for (long i = 0; i < n_fft; i++) {
        long s = s1 + i;
        int32_t sample = (s >= 0 && s < n_samples) ? td_signals[size_t(s * n_chan + c)] : 0;
        frame[i] = saturate16(round_shift(int64_t(sample) * window[i], 15));  // Q15 * Q15
      }
      fft.forward(frame.data(), re.data(), im.data());
      for (long f = 0; f < out.n_freq; f++) {
        out.re[out.idx(f, c, t)] = re[f];
        out.im[out.idx(f, c, t)] = im[f];
      }
    }
  }
  return out;
}

// Same weighted overlap-add as piva::impl::istft. Returns (n_samples, n_chan)
// Q15 with n_samples = n_frames * n_shift - (N - n_shift).
inline std::vector<int16_t> istft(const Stft& X, const SynthesisWindow& window, long n_shift,
                                  long& n_samples_out) {
  long n_fft = long(window.w.size());
  if (n_fft / 2 + 1 != X.n_freq) throw std::invalid_argument("window does not match the STFT");
  FFT fft{size_t(n_fft)};

  long n_samples = X.n_frames * n_shift + (n_shift - n_fft);
  n_samples_out = n_samples;
  std::vector<int32_t> acc(size_t(n_samples * X.n_chan), 0);  // overlap-add headroom
  std::vector<int16_t> re(X.n_freq), im(X.n_freq);
  std::vector<int32_t> frame(n_fft);

  for (long t = 0; t < X.n_frames; t++) {
    long s2 = t * n_shift + n_shift;
    long s1 = s2 - n_fft;
    for (long c = 0; c < X.n_chan; c++) {
      for (long f = 0; f < X.n_freq; f++) {
        re[f] = X.re[X.idx(f, c, t)];
        im[f] = X.im[X.idx(f, c, t)];
      }
      fft.inverse(re.data(), im.data(), frame.data());
      for (long i = 0; i < n_fft; i++) {
        long s = s1 + i;
        if (s < 0 || s >= n_samples) continue;
        acc[size_t(s * X.n_chan + c)] +=
            int32_t(round_shift(int64_t(frame[i]) * window.w[i], window.frac_bits));
      }
    }
  }

  std::vector<int16_t> out(acc.size());
  for (size_t i = 0; i < acc.size(); i++) out[i] = saturate16(acc[i]);
  return out;
}

}  // namespace q15
