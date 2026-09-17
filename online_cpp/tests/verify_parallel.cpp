// Checks that the OpenMP (parallel over frequency bins) path gives the same
// result as the serial path, and measures the speed-up.
//
// Build:  ./build.sh   (the script builds this too), or by hand:
//   c++ -std=c++17 -O3 -Iinclude tests/verify_parallel.cpp -fopenmp -o verify_parallel
#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <random>
#include <vector>

#include "piva_online/auxiva_ip_online.hpp"
#include "piva_online/auxiva_iss_online.hpp"

using namespace piva_online;

static std::unique_ptr<OnlineAuxIVABase> make(const std::string& algo, size_t n_freq,
                                              size_t n_chan, bool parallel) {
  std::unique_ptr<OnlineAuxIVABase> bss;
  if (algo == "ip")
    bss = std::make_unique<AuxIVAIPOnline>(n_freq, n_chan, 3, 0.96);
  else
    bss = std::make_unique<AuxIVAISSOnline>(n_freq, n_chan, 3, 0.96);
  bss->set_parallel(parallel);
  return bss;
}

// returns seconds per frame, and writes the last frame's output into y
static double run(const std::string& algo, size_t n_freq, size_t n_chan, bool parallel,
                  const std::vector<cplx>& X, size_t n_frames, std::vector<cplx>& y) {
  auto bss = make(algo, n_freq, n_chan, parallel);
  std::vector<cplx> out(n_freq * n_chan);
  size_t stride = n_freq * n_chan;

  auto tic = std::chrono::steady_clock::now();
  for (size_t t = 0; t < n_frames; t++) {
    bss->process_frame(&X[t * stride], out.data());
    bss->project_back_frame(out.data(), y.data(), 0);
  }
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - tic).count() /
         double(n_frames);
}

int main() {
#ifdef _OPENMP
  std::printf("OpenMP enabled, %d threads available\n\n", max_threads());
#else
  std::printf("Built WITHOUT OpenMP: everything runs serially\n\n");
#endif
  std::printf("%-5s %6s %8s %12s %12s %9s %12s\n", "algo", "n_fft", "channels", "serial ms",
              "parallel ms", "speed-up", "max rel diff");

  std::mt19937 rng(0);
  std::normal_distribution<double> nd;
  const size_t n_frames = 60;

  for (size_t n_fft : {1024u, 4096u}) {
    size_t n_freq = n_fft / 2 + 1;
    for (size_t n_chan : {2u, 4u, 8u}) {
      std::vector<cplx> X(n_frames * n_freq * n_chan);
      for (auto& v : X) v = cplx(nd(rng), nd(rng));

      std::vector<cplx> y_ser(n_freq * n_chan), y_par(n_freq * n_chan);
      double t_ser = run("iss", n_freq, n_chan, false, X, n_frames, y_ser);
      double t_par = run("iss", n_freq, n_chan, true, X, n_frames, y_par);

      double num = 0.0, den = 0.0;
      for (size_t i = 0; i < y_ser.size(); i++) {
        num = std::max(num, std::abs(y_par[i] - y_ser[i]));
        den = std::max(den, std::abs(y_ser[i]));
      }
      std::printf("%-5s %6zu %8zu %12.3f %12.3f %8.2fx %12.1e\n", "iss", n_fft, n_chan,
                  1e3 * t_ser, 1e3 * t_par, t_ser / t_par, den > 0 ? num / den : 0.0);

      t_ser = run("ip", n_freq, n_chan, false, X, n_frames, y_ser);
      t_par = run("ip", n_freq, n_chan, true, X, n_frames, y_par);
      num = den = 0.0;
      for (size_t i = 0; i < y_ser.size(); i++) {
        num = std::max(num, std::abs(y_par[i] - y_ser[i]));
        den = std::max(den, std::abs(y_ser[i]));
      }
      std::printf("%-5s %6zu %8zu %12.3f %12.3f %8.2fx %12.1e\n", "ip", n_fft, n_chan,
                  1e3 * t_ser, 1e3 * t_par, t_ser / t_par, den > 0 ? num / den : 0.0);
    }
  }
  std::printf(
      "\n(differences come only from the order of the r_k sum in the covariance\n"
      " update; every frequency bin itself is computed identically)\n");
}
