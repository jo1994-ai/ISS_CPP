// Reading and writing 16-bit PCM WAV files directly as int16_t (Q15) samples.
#pragma once

#include <cstdint>
#include <cstring>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

namespace wav16 {

struct Wav {
  int samplerate = 0;
  long n_chan = 0;
  long n_samples = 0;
  std::vector<int16_t> samples;  // interleaved, (n_samples, n_chan)
};

inline uint32_t u32(const unsigned char* p) {
  return uint32_t(p[0]) | (uint32_t(p[1]) << 8) | (uint32_t(p[2]) << 16) | (uint32_t(p[3]) << 24);
}
inline uint16_t u16(const unsigned char* p) { return uint16_t(p[0] | (p[1] << 8)); }

inline Wav read(const std::string& filename) {
  std::ifstream in(filename, std::ios::binary);
  if (!in) throw std::runtime_error("Could not open " + filename);
  std::vector<unsigned char> b((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  if (b.size() < 12 || std::memcmp(b.data(), "RIFF", 4) || std::memcmp(b.data() + 8, "WAVE", 4))
    throw std::runtime_error(filename + " is not a WAV file");

  uint16_t format = 0, n_chan = 0, bits = 0;
  uint32_t rate = 0;
  const unsigned char* data = nullptr;
  size_t data_size = 0;
  for (size_t pos = 12; pos + 8 <= b.size();) {
    const unsigned char* chunk = b.data() + pos;
    uint32_t size = u32(chunk + 4);
    if (!std::memcmp(chunk, "fmt ", 4) && size >= 16) {
      format = u16(chunk + 8);
      n_chan = u16(chunk + 10);
      rate = u32(chunk + 12);
      bits = u16(chunk + 22);
      if (format == 0xFFFE && size >= 26) format = u16(chunk + 32);
    } else if (!std::memcmp(chunk, "data", 4)) {
      data = chunk + 8;
      data_size = std::min<size_t>(size, b.size() - pos - 8);
    }
    pos += 8 + size + (size & 1);
  }
  if (!data || n_chan == 0) throw std::runtime_error(filename + ": missing fmt or data chunk");
  if (format != 1 || bits != 16)
    throw std::runtime_error(filename + ": only 16-bit PCM WAV files are supported");

  Wav w;
  w.samplerate = int(rate);
  w.n_chan = n_chan;
  w.n_samples = long(data_size / (2 * n_chan));
  w.samples.resize(size_t(w.n_samples * n_chan));
  for (size_t i = 0; i < w.samples.size(); i++) w.samples[i] = int16_t(u16(data + 2 * i));
  return w;
}

inline void write(const std::string& filename, const std::vector<int16_t>& samples, long n_chan,
                  int samplerate) {
  std::ofstream out(filename, std::ios::binary);
  if (!out) throw std::runtime_error("Could not write " + filename);
  auto put32 = [&](uint32_t v) { for (int i = 0; i < 4; i++) out.put(char((v >> (8 * i)) & 0xff)); };
  auto put16 = [&](uint16_t v) { out.put(char(v & 0xff)); out.put(char(v >> 8)); };

  uint32_t data_size = uint32_t(samples.size() * 2);
  out.write("RIFF", 4);
  put32(36 + data_size);
  out.write("WAVEfmt ", 8);
  put32(16);
  put16(1);
  put16(uint16_t(n_chan));
  put32(uint32_t(samplerate));
  put32(uint32_t(samplerate * n_chan * 2));
  put16(uint16_t(n_chan * 2));
  put16(16);
  out.write("data", 4);
  put32(data_size);
  for (int16_t s : samples) put16(uint16_t(s));
}

}  // namespace wav16
