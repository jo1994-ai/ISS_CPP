#!/bin/sh
# Build and run fft_example using the FFTW bundled in third_party/.
# Nothing needs to be downloaded or installed except a C/C++ compiler
# (macOS: xcode-select --install, Linux: sudo apt install build-essential).
#
# The first build compiles FFTW from source (1-3 minutes) into
# third_party/fftw/, later builds reuse it.
set -e
cd "$(dirname "$0")"

FFTW_VERSION=3.3.11
ROOT=$(pwd)
TP=$ROOT/third_party
FFTW=$TP/fftw
CC=${CC:-cc}
CXX=${CXX:-c++}

# compile and link against the bundled static FFTW; fails if it's missing or
# was built for a different machine
build_example() {
  [ -f "$FFTW/lib/libfftw3.a" ] || return 1
  $CXX -std=c++17 -O2 -Wall -Wextra fft_example.cpp -I"$FFTW/include" \
      "$FFTW/lib/libfftw3.a" -lm -o fft_example 2>"$TP/link.log"
}

build_fftw() {
  echo "Building FFTW $FFTW_VERSION from third_party/ (first time only, 1-3 minutes)..."
  rm -rf "$FFTW" "$TP/fftw-$FFTW_VERSION"
  tar -xzf "$TP/fftw-$FFTW_VERSION.tar.gz" -C "$TP"
  cd "$TP/fftw-$FFTW_VERSION"

  # SIMD speeds up the FFT; fall back to a plain build if unsupported
  case "$(uname -m)" in
    arm64|aarch64) SIMD="--enable-neon" ;;
    x86_64|amd64) SIMD="--enable-sse2 --enable-avx --enable-avx2" ;;
    *) SIMD="" ;;
  esac
  COMMON="--prefix=$FFTW --enable-static --disable-shared --disable-fortran --disable-doc CC=$CC"

  if ! ./configure $COMMON $SIMD >"$TP/configure.log" 2>&1; then
    echo "  SIMD configure failed, building without SIMD"
    ./configure $COMMON >"$TP/configure.log" 2>&1
  fi
  JOBS=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)
  make -j"$JOBS" >"$TP/make.log" 2>&1
  make install >>"$TP/make.log" 2>&1

  cd "$ROOT"
  rm -rf "$TP/fftw-$FFTW_VERSION"  # keep only the installed result
  echo "FFTW built into third_party/fftw/"
}

if ! build_example; then
  build_fftw
  if ! build_example; then
    cat "$TP/link.log"
    exit 1
  fi
fi
rm -f "$TP/link.log"

echo "Built ./fft_example"
mkdir -p output
./fft_example
