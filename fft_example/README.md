fft_example
===========

Standalone STFT / iSTFT in C++ with FFTW. The analysis and synthesis
match the original piva C++ code (`include/piva/stft.hpp`):

* analysis window: periodic Hann
* synthesis window: dual window (`make_dual`), for perfect reconstruction
* `stft()` returns `(n_freq, n_chan, n_frames)`
* `istft()` does weighted overlap-add, inverse FFT scaled by 1/n_fft

Checked against `piva::stft` / `piva::istft`: identical output.

The demo reads a speech recording (`audio/mixture.wav`, the same 2-mic
mixture used for the online IVA algorithms), transforms it to the STFT
domain and back, and writes the reconstructed audio.

Folder contents
---------------

    fft_example.cpp                    STFT / iSTFT and the demo
    wav_io.hpp                         WAV reading / writing
    audio/mixture.wav                  2-channel speech mixture, 16 kHz
    third_party/fftw-3.3.11.tar.gz     FFTW source, from https://fftw.org
                                       (sha256 5630c24c...f239a1)
    third_party/fftw/                  FFTW compiled by build.sh (static library)
    build.sh, Makefile

Dependencies
------------

FFTW is included, so nothing has to be downloaded. You only need a C/C++
compiler:

| OS            | Install                                |
|---------------|----------------------------------------|
| macOS         | `xcode-select --install`               |
| Ubuntu/Debian | `sudo apt install build-essential`     |
| Fedora        | `sudo dnf install gcc gcc-c++ make`    |

Build and run
-------------

    ./build.sh

or

    make run

The first time on a machine, `build.sh` compiles FFTW from
`third_party/fftw-3.3.11.tar.gz` (about 1 minute, with SIMD: NEON on ARM,
AVX2 on Intel/AMD) and installs it into `third_party/fftw/`. After that a
build takes under a second. If `third_party/fftw/` was built on a different
kind of machine, it is rebuilt automatically. FFTW is linked statically, so
the `fft_example` program does not need FFTW installed to run.

`make distclean` removes the compiled FFTW and the program, to start fresh.

Expected output:

    Input audio/mixture.wav: 2 channels, 72764 samples at 16000 Hz
    STFT shape (n_freq, n_chan, n_frames): (513, 2, 144)
    iSTFT output: 73216 samples x 2 channels (input: 72764 samples)
    Reconstruction max error: 3.33e-16
    Wrote output/reconstructed.wav

`output/reconstructed.wav` is identical to the input, sample for sample.

Other files and frame sizes:

    ./fft_example <input.wav> <output.wav> [n_analysis=1024] [n_shift=n_analysis/2]
    ./fft_example audio/mixture.wav output/recon_4096.wav 4096 2048

Using it
--------

    auto aw = hann(n_analysis);
    auto sw = make_dual(aw, n_shift);
    Stft X  = stft(x, n_samples, n_chan, aw, n_shift, 0, 0);
    // ... process X ...
    auto y  = istft(X, sw, n_shift, 0, 0, n_out);

`x` and `y` are `(n_samples, n_chan)` stored row-major:
`x[s * n_chan + c]`.

License
-------

FFTW is GPL-2.0-or-later (see `COPYING` inside the FFTW tarball). Fine for
research use; a closed-source product needs a commercial FFTW license.
