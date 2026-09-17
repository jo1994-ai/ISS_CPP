// Streaming STFT / iSTFT for multichannel audio, built on FFTW.
//
// StftAnalysis takes `hop` new samples per channel and returns one STFT
// frame. StftSynthesis takes one frame and returns `hop` output samples.
// The analysis window is a periodic Hann window and the synthesis window is
// its dual (same construction as piva::windows::make_dual), so a
// pass-through reconstructs the input delayed by n_fft - hop samples.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
#pragma once

#include <fftw3.h>

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include "piva_online/common.hpp"

namespace piva_online {

inline std::vector<double> hann_window(size_t length) {
  std::vector<double> w(length);
  for (size_t n = 0; n < length; n++)
    w[n] = 0.5 * (1.0 - std::cos(2.0 * M_PI * double(n) / double(length)));
  return w;
}

inline std::vector<double> dual_window(const std::vector<double>& awin, size_t hop) {
  long length = long(awin.size());
  std::vector<double> norm(length, 0.0);
  for (long shift = -long((length - 1) / hop) * long(hop); shift < length; shift += long(hop))
    for (long m = std::max(0L, shift); m < std::min(length, length + shift); m++)
      norm[m] += awin[m - shift] * awin[m - shift];
  std::vector<double> swin(length);
  for (long m = 0; m < length; m++) swin[m] = awin[m] / norm[m];
  return swin;
}

class StftAnalysis {
 public:
  StftAnalysis(size_t n_fft, size_t hop, size_t n_chan)
      : n_fft_(n_fft),
        hop_(hop),
        n_chan_(n_chan),
        n_freq_(n_fft / 2 + 1),
        window_(hann_window(n_fft)),
        buffer_(n_chan * n_fft, 0.0) {
    if (hop == 0 || hop > n_fft) throw std::invalid_argument("invalid hop size");
    in_ = fftw_alloc_real(n_fft);
    out_ = fftw_alloc_complex(n_freq_);
    plan_ = fftw_plan_dft_r2c_1d(int(n_fft), in_, out_, FFTW_MEASURE);
  }
  ~StftAnalysis() {
    fftw_destroy_plan(plan_);
    fftw_free(in_);
    fftw_free(out_);
  }
  StftAnalysis(const StftAnalysis&) = delete;
  StftAnalysis& operator=(const StftAnalysis&) = delete;

  // samples: (hop, n_chan) interleaved; frame: (n_freq, n_chan)
  void push(const double* samples, cplx* frame) {
    for (size_t c = 0; c < n_chan_; c++) {
      double* buf = &buffer_[c * n_fft_];
      std::copy(buf + hop_, buf + n_fft_, buf);
      for (size_t i = 0; i < hop_; i++) buf[n_fft_ - hop_ + i] = samples[i * n_chan_ + c];

      for (size_t i = 0; i < n_fft_; i++) in_[i] = buf[i] * window_[i];
      fftw_execute(plan_);
      for (size_t f = 0; f < n_freq_; f++)
        frame[f * n_chan_ + c] = cplx(out_[f][0], out_[f][1]);
    }
  }

  size_t n_freq() const { return n_freq_; }

 private:
  size_t n_fft_, hop_, n_chan_, n_freq_;
  std::vector<double> window_, buffer_;
  double* in_;
  fftw_complex* out_;
  fftw_plan plan_;
};

class StftSynthesis {
 public:
  StftSynthesis(size_t n_fft, size_t hop, size_t n_chan)
      : n_fft_(n_fft),
        hop_(hop),
        n_chan_(n_chan),
        n_freq_(n_fft / 2 + 1),
        window_(dual_window(hann_window(n_fft), hop)),
        accum_(n_chan * n_fft, 0.0) {
    in_ = fftw_alloc_complex(n_freq_);
    out_ = fftw_alloc_real(n_fft);
    plan_ = fftw_plan_dft_c2r_1d(int(n_fft), in_, out_, FFTW_MEASURE);
  }
  ~StftSynthesis() {
    fftw_destroy_plan(plan_);
    fftw_free(in_);
    fftw_free(out_);
  }
  StftSynthesis(const StftSynthesis&) = delete;
  StftSynthesis& operator=(const StftSynthesis&) = delete;

  // frame: (n_freq, n_chan); samples: (hop, n_chan) interleaved
  void push(const cplx* frame, double* samples) {
    double inv_n = 1.0 / double(n_fft_);  // FFTW's inverse is unnormalized
    for (size_t c = 0; c < n_chan_; c++) {
      for (size_t f = 0; f < n_freq_; f++) {
        in_[f][0] = frame[f * n_chan_ + c].real();
        in_[f][1] = frame[f * n_chan_ + c].imag();
      }
      fftw_execute(plan_);

      double* acc = &accum_[c * n_fft_];
      for (size_t i = 0; i < n_fft_; i++) acc[i] += out_[i] * inv_n * window_[i];
      for (size_t i = 0; i < hop_; i++) samples[i * n_chan_ + c] = acc[i];
      std::copy(acc + hop_, acc + n_fft_, acc);
      std::fill(acc + n_fft_ - hop_, acc + n_fft_, 0.0);
    }
  }

 private:
  size_t n_fft_, hop_, n_chan_, n_freq_;
  std::vector<double> window_, accum_;
  fftw_complex* in_;
  double* out_;
  fftw_plan plan_;
};

}  // namespace piva_online
