// Online (frame-recursive) AuxIVA with iterative projection (IP).
// C++ port of piva/auxiva_ip_online.py.
//
// IP update for source k, per frequency bin:
//   solve (W V_k) w_k = e_k,  then  w_k <- w_k / sqrt(w_k^H V_k w_k)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
#pragma once

#include "piva_online/common.hpp"

namespace piva_online {

class AuxIVAIPOnline : public OnlineAuxIVABase {
 public:
  AuxIVAIPOnline(size_t n_freq, size_t n_chan, size_t n_iter = 3,
                 double alpha = 0.96, Model model = Model::Laplace,
                 double eps = 1e-10)
      : OnlineAuxIVABase(n_freq, n_chan, n_iter, alpha, model, eps),
        WV_(n_chan * n_chan),
        w_(n_chan) {}

 protected:
  void update_demixing(size_t k) override {
    size_t n = n_chan_;

    for (size_t f = 0; f < n_freq_; f++) {
      const cplx* Vk = &V_[v_idx(k, f, 0, 0)];

      // WV = W V_k, row m is w_m^H V_k
      for (size_t r = 0; r < n; r++)
        for (size_t c = 0; c < n; c++) {
          cplx s = 0.0;
          for (size_t d = 0; d < n; d++) s += W_[w_idx(f, r, d)] * Vk[d * n + c];
          WV_[r * n + c] = s;
        }

      // w = (W V_k)^{-1} e_k, the true (unconjugated) vector
      for (size_t c = 0; c < n; c++) w_[c] = (c == k) ? 1.0 : 0.0;
      if (!solve_inplace(WV_.data(), w_.data(), n)) continue;

      // normalize by sqrt(w^H V_k w)
      double denom = 0.0;
      for (size_t r = 0; r < n; r++) {
        cplx Vw = 0.0;
        for (size_t c = 0; c < n; c++) Vw += Vk[r * n + c] * w_[c];
        denom += std::real(std::conj(w_[r]) * Vw);
      }
      double scale = 1.0 / std::sqrt(std::max(denom, eps_));

      // store back as a row, i.e. as w_k^H
      for (size_t c = 0; c < n; c++) W_[w_idx(f, k, c)] = std::conj(w_[c] * scale);
    }
  }

 private:
  std::vector<cplx> WV_, w_;
};

}  // namespace piva_online
