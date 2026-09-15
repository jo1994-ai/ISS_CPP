#!/bin/sh
# Build and run the Q15 (16-bit fixed-point) FFT example.
# Needs only a C++17 compiler (macOS: xcode-select --install,
# Linux: sudo apt install build-essential). No FFTW or other libraries.
set -e
cd "$(dirname "$0")"

CXX=${CXX:-c++}
$CXX -std=c++17 -O2 -Wall -Wextra fft_example_q15.cpp -o fft_example_q15
echo "Built ./fft_example_q15"
./fft_example_q15 "$@"
