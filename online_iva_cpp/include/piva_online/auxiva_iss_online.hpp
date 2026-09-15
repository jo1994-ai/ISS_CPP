// Online (frame-recursive) AuxIVA with iterative source steering (ISS).
// C++ port of piva/auxiva_iss_online.py (arXiv:2209.00937, Algorithm 1).
//
// ISS update for pivot source k, per frequency bin, for every source m:
//   v_m = (w_m^H V_m w_k) / (w_k^H V_m w_k)        for m != k
//   v_k = 1 - 1 / sqrt(w_k^H V_k w_k)
//   w_m^H <- w_m^H - v_m w_k^H
// No matrix inversion or linear solve is needed.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
#pragma once

#include "piva_online/common.hpp"

namespace piva_online {

class AuxIVAISSOnline : public OnlineAuxIVABase {
 public:
  AuxIVAISSOnline(size_t n_freq, size_t n_chan, size_t n_iter = 3,
                  double alpha = 0.96, Model model = Model::Laplace,
                  double eps = 1e-10)
      : OnlineAuxIVABase(n_freq, n_chan, n_iter, alpha, model, eps),
        wk_(n_chan),
        v_(n_chan) {}

 protected:
  void update_demixing(size_t k) override {
    size_t n = n_chan_;

    for (size_t f = 0; f < n_freq_; f++) {
      // frozen pivot row (stores w_k^H)
      for (size_t c = 0; c < n; c++) wk_[c] = W_[w_idx(f, k, c)];

      // all v_m are computed from the W before this update
      for (size_t m = 0; m < n_src_; m++) {
        const cplx* Vm = &V_[v_idx(m, f, 0, 0)];
        double denom = 0.0;
        cplx numer = 0.0;
        for (size_t r = 0; r < n; r++) {
          cplx Vw = 0.0;  // (V_m conj(W[f,k,:]))_r
          for (size_t c = 0; c < n; c++) Vw += Vm[r * n + c] * std::conj(wk_[c]);
          denom += std::real(wk_[r] * Vw);
          numer += W_[w_idx(f, m, r)] * Vw;
        }
        denom = std::max(denom, eps_);
        v_[m] = (m == k) ? cplx(1.0 - 1.0 / std::sqrt(denom)) : numer / denom;
      }

      for (size_t m = 0; m < n_src_; m++)
        for (size_t c = 0; c < n; c++) W_[w_idx(f, m, c)] -= v_[m] * wk_[c];
    }
  }

 private:
  std::vector<cplx> wk_, v_;
};

}  // namespace piva_online
