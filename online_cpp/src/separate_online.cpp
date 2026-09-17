// Streaming blind source separation of a multichannel wav file with online
// AuxIVA-IP or AuxIVA-ISS.
//
// Audio is fed hop by hop, as it would arrive from a sound card:
//   hop samples -> streaming STFT -> online IVA -> projection back
//               -> streaming iSTFT -> hop output samples
//
// Usage:
//   separate_online <input.wav> <output_prefix> <ip|iss> [n_fft=4096]
//                   [hop=n_fft/2] [n_iter=3] [alpha=0.96] [model=laplace]
//
// Writes <output_prefix>_src<k>.wav for every source, time-aligned with
// the input (the n_fft - hop samples of STFT latency are removed).
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
#include <sndfile.h>

#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "piva_online/auxiva_ip_online.hpp"
#include "piva_online/auxiva_iss_online.hpp"
#include "piva_online/stft_online.hpp"

using namespace piva_online;

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0]
              << " <input.wav> <output_prefix> <ip|iss> [n_fft=4096] [hop=n_fft/2]"
                 " [n_iter=3] [alpha=0.96] [model=laplace]\n";
    return 1;
  }
  std::string input_filename(argv[1]);
  std::string output_prefix(argv[2]);
  std::string algo(argv[3]);
  size_t n_fft = argc > 4 ? std::stoul(argv[4]) : 4096;
  size_t hop = argc > 5 ? std::stoul(argv[5]) : n_fft / 2;
  size_t n_iter = argc > 6 ? std::stoul(argv[6]) : 3;
  double alpha = argc > 7 ? std::stod(argv[7]) : 0.96;
  Model model = model_from_string(argc > 8 ? argv[8] : "laplace");

  // read the input file
  SF_INFO info{};
  SNDFILE* in_file = sf_open(input_filename.c_str(), SFM_READ, &info);
  if (!in_file) {
    std::cerr << "Could not open " << input_filename << ": " << sf_strerror(nullptr) << "\n";
    return 1;
  }
  size_t n_chan = size_t(info.channels);
  size_t n_samples = size_t(info.frames);
  std::vector<double> audio(n_samples * n_chan);
  sf_readf_double(in_file, audio.data(), info.frames);
  sf_close(in_file);

  std::cout << "Input: " << n_chan << " channels, " << n_samples << " samples at "
            << info.samplerate << " Hz\n";

  std::unique_ptr<OnlineAuxIVABase> bss;
  size_t n_freq = n_fft / 2 + 1;
  if (algo == "ip")
    bss = std::make_unique<AuxIVAIPOnline>(n_freq, n_chan, n_iter, alpha, model);
  else if (algo == "iss")
    bss = std::make_unique<AuxIVAISSOnline>(n_freq, n_chan, n_iter, alpha, model);
  else {
    std::cerr << "Unknown algorithm " << algo << " (use ip or iss)\n";
    return 1;
  }

  StftAnalysis analysis(n_fft, hop, n_chan);
  StftSynthesis synthesis(n_fft, hop, n_chan);

  size_t latency = n_fft - hop;
  size_t n_frames = (n_samples + latency + hop - 1) / hop;

  std::vector<double> in_hop(hop * n_chan), out_hop(hop * n_chan);
  std::vector<cplx> X(n_freq * n_chan), Y(n_freq * n_chan), Y_pb(n_freq * n_chan);
  std::vector<double> output(n_frames * hop * n_chan, 0.0);

  double bss_seconds = 0.0;

  for (size_t t = 0; t < n_frames; t++) {
    // next hop of input, zero-padded past the end of the file
    for (size_t i = 0; i < hop; i++) {
      size_t s = t * hop + i;
      for (size_t c = 0; c < n_chan; c++)
        in_hop[i * n_chan + c] = s < n_samples ? audio[s * n_chan + c] : 0.0;
    }

    analysis.push(in_hop.data(), X.data());

    auto tic = std::chrono::steady_clock::now();
    bss->process_frame(X.data(), Y.data());
    bss->project_back_frame(Y.data(), Y_pb.data(), 0);
    bss_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - tic).count();

    synthesis.push(Y_pb.data(), out_hop.data());
    std::copy(out_hop.begin(), out_hop.end(), output.begin() + long(t * hop * n_chan));
  }

  double frame_ms = 1000.0 * bss_seconds / double(n_frames);
  double hop_ms = 1000.0 * double(hop) / double(info.samplerate);
  std::cout << "Online AuxIVA-" << (algo == "ip" ? "IP" : "ISS") << ": " << n_frames
            << " frames, " << frame_ms << " ms/frame (hop is " << hop_ms
            << " ms, real-time factor " << frame_ms / hop_ms << ")\n";

  // write one file per source, peak-normalized to 0.9 like the Python examples
  for (size_t k = 0; k < n_chan; k++) {
    std::vector<double> src(n_samples);
    double peak = 0.0;
    for (size_t s = 0; s < n_samples; s++) {
      src[s] = output[(s + latency) * n_chan + k];
      peak = std::max(peak, std::abs(src[s]));
    }
    if (peak > 0.0)
      for (auto& v : src) v *= 0.9 / peak;

    std::string filename = output_prefix + "_src" + std::to_string(k + 1) + ".wav";
    SF_INFO out_info{};
    out_info.samplerate = info.samplerate;
    out_info.channels = 1;
    out_info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
    SNDFILE* out_file = sf_open(filename.c_str(), SFM_WRITE, &out_info);
    if (!out_file) {
      std::cerr << "Could not write " << filename << "\n";
      return 1;
    }
    sf_writef_double(out_file, src.data(), sf_count_t(n_samples));
    sf_close(out_file);
    std::cout << "Wrote " << filename << "\n";
  }
}
