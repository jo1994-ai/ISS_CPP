// Shared pieces of the online AuxIVA-IP / AuxIVA-ISS C++ implementations.
//
// These are C++ ports of piva/auxiva_ip_online.py and
// piva/auxiva_iss_online.py. Both algorithms share the recursive weighted
// covariance, the auxiliary variable and the per-frame iteration scheme;
// they only differ in the demixing update rule (see OnlineAuxIVABase and
// the two derived classes).
//
// Memory layout (row-major, same axis order as the Python code):
//   x, y : (n_freq, n_chan)                 one STFT frame
//   W    : (n_freq, n_chan, n_chan)         row k stores w_k^H, so y = W x
//   V    : (n_src, n_freq, n_chan, n_chan)  per-source weighted covariance
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace piva_online {

// Frequency bins are independent, so they can be processed in parallel.
// Compiled with OpenMP (-fopenmp) this uses threads, otherwise everything
// below falls back to plain serial loops.
inline int max_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}
inline int this_thread() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

using cplx = std::complex<double>;

enum class Model { Laplace, Gauss };

inline Model model_from_string(const std::string& name) {
  if (name == "laplace") return Model::Laplace;
  if (name == "gauss") return Model::Gauss;
  throw std::invalid_argument("No such model " + name);
}

// Auxiliary-variable weight, identical to _phi in auxiva_iss_online.py
inline double phi(double r, Model model, size_t n_freq, double eps) {
  if (model == Model::Laplace) return 1.0 / std::max(eps, 2.0 * r);
  return 1.0 / std::max(eps, (r * r) / double(n_freq));
}

// Solve A x = b in place for a small dense complex system with partial
// pivoting. A (n x n, row-major) is destroyed; the solution is left in b.
// Returns false if the matrix is numerically singular.
inline bool solve_inplace(cplx* A, cplx* b, size_t n) {
  for (size_t col = 0; col < n; col++) {
    size_t piv = col;
    double best = std::abs(A[col * n + col]);
    for (size_t r = col + 1; r < n; r++) {
      double v = std::abs(A[r * n + col]);
      if (v > best) {
        best = v;
        piv = r;
      }
    }
    if (best == 0.0) return false;
    if (piv != col) {
      for (size_t c = 0; c < n; c++) std::swap(A[col * n + c], A[piv * n + c]);
      std::swap(b[col], b[piv]);
    }
    for (size_t r = col + 1; r < n; r++) {
      cplx f = A[r * n + col] / A[col * n + col];
      if (f == cplx(0.0)) continue;
      for (size_t c = col; c < n; c++) A[r * n + c] -= f * A[col * n + c];
      b[r] -= f * b[col];
    }
  }
  for (size_t i = n; i-- > 0;) {
    cplx s = b[i];
    for (size_t c = i + 1; c < n; c++) s -= A[i * n + c] * b[c];
    b[i] = s / A[i * n + i];
  }
  return true;
}

class OnlineAuxIVABase {
 public:
  OnlineAuxIVABase(size_t n_freq, size_t n_chan, size_t n_iter = 3,
                   double alpha = 0.96, Model model = Model::Laplace,
                   double eps = 1e-10)
      : n_freq_(n_freq),
        n_chan_(n_chan),
        n_src_(n_chan),
        n_iter_(n_iter),
        alpha_(alpha),
        model_(model),
        eps_(eps),
        W_(n_freq * n_chan * n_chan),
        V_(n_chan * n_freq * n_chan * n_chan),
        V_prev_(V_.size()),
        xxH_(n_freq * n_chan * n_chan),
        a_(size_t(max_threads()), std::vector<cplx>(n_chan * n_chan)),
        e_(size_t(max_threads()), std::vector<cplx>(n_chan)) {
    // Threads only pay off once there is real work per frame: the cost per
    // bin grows with n_chan^3, and starting a parallel region costs tens of
    // microseconds. Measured on an 8-core M1: 2 channels is ~2x SLOWER with
    // threads, 4 channels at n_fft 1024 breaks even, 8 channels gains ~3x.
    // The threshold below keeps the small cases serial; set_parallel()
    // overrides it.
    parallel_ = n_freq_ * n_chan_ * n_chan_ * n_chan_ >= 100000 && max_threads() > 1;
    reset();
  }

  // Override the automatic choice (mainly for benchmarking)
  void set_parallel(bool on) { parallel_ = on; }
  bool parallel() const { return parallel_; }

  virtual ~OnlineAuxIVABase() = default;

  // W = identity per bin, V = 0.01 * identity per source/bin
  void reset() {
    size_t n = n_chan_;
    std::fill(W_.begin(), W_.end(), cplx(0.0));
    std::fill(V_.begin(), V_.end(), cplx(0.0));
    for (size_t f = 0; f < n_freq_; f++)
      for (size_t c = 0; c < n; c++) W_[w_idx(f, c, c)] = 1.0;
    for (size_t k = 0; k < n_src_; k++)
      for (size_t f = 0; f < n_freq_; f++)
        for (size_t c = 0; c < n; c++) V_[v_idx(k, f, c, c)] = 0.01;
  }

  // Process one STFT frame. x and y are (n_freq, n_chan), row-major.
  // y is the demixed frame, without any scale correction.
  void process_frame(const cplx* x, cplx* y) {
    size_t n = n_chan_;

    // x x^H per frequency bin
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (parallel_)
#endif
    for (long long ff = 0; ff < (long long)n_freq_; ff++) {
      size_t f = size_t(ff);
      for (size_t c = 0; c < n; c++)
        for (size_t d = 0; d < n; d++)
          xxH_[(f * n + c) * n + d] = x[f * n + c] * std::conj(x[f * n + d]);
    }

    // snapshot of V_{t-1}; the recursion is anchored to the previous frame
    std::copy(V_.begin(), V_.end(), V_prev_.begin());

    for (size_t it = 0; it < n_iter_; it++) {
      for (size_t k = 0; k < n_src_; k++) {
        update_covariance(x, k);
        update_demixing(k);
      }
    }

    demix(x, y);
  }

  // Scale-corrected output by the minimum distortion principle:
  // y_k <- A[ref, k] y_k with A = W^{-1}, using the current W only, so it
  // stays causal. Writes (n_freq, n_src) into y_scaled.
  void project_back_frame(const cplx* y, cplx* y_scaled, size_t ref = 0) {
    size_t n = n_chan_;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (parallel_)
#endif
    for (long long ff = 0; ff < (long long)n_freq_; ff++) {
      size_t f = size_t(ff);
      cplx* a = a_[size_t(this_thread())].data();  // scratch, one per thread
      cplx* e = e_[size_t(this_thread())].data();
      // column k of A solves W a = e_k; row `ref` of A solves W^T z = e_ref
      for (size_t r = 0; r < n; r++)
        for (size_t c = 0; c < n; c++) a[r * n + c] = W_[w_idx(f, c, r)];
      for (size_t c = 0; c < n; c++) e[c] = (c == ref) ? 1.0 : 0.0;
      bool ok = solve_inplace(a, e, n);
      for (size_t k = 0; k < n; k++)
        y_scaled[f * n + k] = ok ? e[k] * y[f * n + k] : y[f * n + k];
    }
  }

  size_t n_freq() const { return n_freq_; }
  size_t n_chan() const { return n_chan_; }
  size_t n_src() const { return n_src_; }
  const std::vector<cplx>& demixing() const { return W_; }

 protected:
  size_t w_idx(size_t f, size_t r, size_t c) const {
    return (f * n_chan_ + r) * n_chan_ + c;
  }
  size_t v_idx(size_t k, size_t f, size_t r, size_t c) const {
    return ((k * n_freq_ + f) * n_chan_ + r) * n_chan_ + c;
  }

  // auxiliary variable and recursive covariance for source k
  void update_covariance(const cplx* x, size_t k) {
    size_t n = n_chan_;
    double r2 = 0.0;
    // r_k sums over all frequencies, so this is a reduction
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+ : r2) if (parallel_)
#endif
    for (long long ff = 0; ff < (long long)n_freq_; ff++) {
      size_t f = size_t(ff);
      cplx yk = 0.0;
      for (size_t c = 0; c < n; c++) yk += W_[w_idx(f, k, c)] * x[f * n + c];
      r2 += std::norm(yk);
    }
    double weight = (1.0 - alpha_) * phi(std::sqrt(r2), model_, n_freq_, eps_);

    size_t block = n * n;
    size_t off = k * n_freq_ * block;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (parallel_)
#endif
    for (long long i = 0; i < (long long)(n_freq_ * block); i++)
      V_[off + size_t(i)] = alpha_ * V_prev_[off + size_t(i)] + weight * xxH_[size_t(i)];
  }

  void demix(const cplx* x, cplx* y) const {
    size_t n = n_chan_;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (parallel_)
#endif
    for (long long ff = 0; ff < (long long)n_freq_; ff++) {
      size_t f = size_t(ff);
      for (size_t r = 0; r < n; r++) {
        cplx s = 0.0;
        for (size_t c = 0; c < n; c++) s += W_[w_idx(f, r, c)] * x[f * n + c];
        y[f * n + r] = s;
      }
    }
  }

  virtual void update_demixing(size_t k) = 0;

  size_t n_freq_, n_chan_, n_src_, n_iter_;
  double alpha_;
  Model model_;
  double eps_;

  bool parallel_ = false;

  std::vector<cplx> W_, V_, V_prev_, xxH_;
  // scratch for projection back, one buffer per thread
  std::vector<std::vector<cplx>> a_, e_;
};

}  // namespace piva_online
