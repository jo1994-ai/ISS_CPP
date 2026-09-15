fft_example (Windows)
=====================

Windows version of `fft_example`: STFT / iSTFT in C++ with FFTW, using the
same analysis and synthesis as the original piva C++ code
(periodic Hann analysis window, dual synthesis window, overlap-add).

It reads a 2-channel speech mixture, goes to the STFT domain and back, and
writes the reconstructed audio, which is identical to the input.

FFTW is included as source code and compiled together with the example, so
nothing has to be downloaded apart from Visual Studio.

Folder contents
---------------

    fft_example.cpp                  STFT / iSTFT and the demo
    wav_io.hpp                       WAV reading / writing
    audio\mixture.wav                2-channel speech mixture, 16 kHz
    third_party\fftw-3.3.11.tar.gz   FFTW source, from https://fftw.org
                                     (sha256 5630c24c...f239a1)
    CMakeLists.txt                   builds FFTW (static) + the example
    build.bat                        build and run in one step

Step 1: install Visual Studio (once)
------------------------------------

1. Download **Visual Studio 2022 Community** (free):
   https://visualstudio.microsoft.com/downloads/
2. In the installer, tick **Desktop development with C++**.
   This includes the C++ compiler and CMake.
3. Install (restart if asked).

Step 2: copy the folder
-----------------------

Put `fft_example_windows` on the Desktop, e.g.
`C:\Users\<you>\Desktop\fft_example_windows`.
If it came as a zip, right-click > **Extract All** first.

Step 3: build and run
---------------------

**Option A, double-click:** open the folder and double-click `build.bat`.

**Option B, Command Prompt:**

    cd %USERPROFILE%\Desktop\fft_example_windows
    build.bat

The first build compiles FFTW (a few minutes), later builds take seconds.
Expected output:

    Input audio/mixture.wav: 2 channels, 72764 samples at 16000 Hz
    STFT shape (n_freq, n_chan, n_frames): (513, 2, 144)
    iSTFT output: 73216 samples x 2 channels (input: 72764 samples)
    Reconstruction max error: 3.33e-16
    Wrote output/reconstructed.wav

Listen to the result:

    start output\reconstructed.wav

Other files and frame sizes (after building once):

    fft_example.exe <input.wav> <output.wav> [n_analysis=1024] [n_shift=n_analysis/2]
    fft_example.exe audio\mixture.wav output\recon_4096.wav 4096 2048

Option C: Visual Studio IDE
---------------------------

**File > Open > Folder...**, choose `fft_example_windows`, wait for CMake to
finish, select `fft_example.exe` as the startup item and press **F5**.

Troubleshooting
---------------

| Message                                | Fix                                                          |
|----------------------------------------|--------------------------------------------------------------|
| `CMake / Visual Studio not found`      | Install Visual Studio with "Desktop development with C++"    |
| `No CMAKE_CXX_COMPILER could be found` | Same as above: the C++ workload is missing                   |
| Build fails after moving the folder    | Delete the `build` folder and run `build.bat` again          |
| Windows SmartScreen blocks `build.bat` | Click **More info > Run anyway** (or run it from cmd)        |

Using it
--------

    auto aw = hann(n_analysis);
    auto sw = make_dual(aw, n_shift);
    Stft X  = stft(x, n_samples, n_chan, aw, n_shift, 0, 0);
    // ... process X ...
    auto y  = istft(X, sw, n_shift, 0, 0, n_out);

`x` and `y` are `(n_samples, n_chan)` stored row-major: `x[s * n_chan + c]`.

License
-------

FFTW is GPL-2.0-or-later (see `COPYING` inside the FFTW tarball). Fine for
research use; a closed-source product needs a commercial FFTW license.
