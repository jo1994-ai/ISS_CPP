Online AuxIVA-IP and AuxIVA-ISS in C++
======================================

C++ ports of `piva/auxiva_ip_online.py` and `piva/auxiva_iss_online.py`,
with a streaming STFT so audio can be processed hop by hop.

Layout
------

    include/piva_online/
        common.hpp             shared state, recursive covariance, projection back
        auxiva_ip_online.hpp   AuxIVAIPOnline  (IP update: per-bin linear solve)
        auxiva_iss_online.hpp  AuxIVAISSOnline (ISS update: inverse-free)
        stft_online.hpp        StftAnalysis / StftSynthesis (FFTW)
    src/separate_online.cpp    streaming separation of a wav file
    python/bindings.cpp        piva_online_cpp Python module
    tests/validate_vs_python.py

Header-only; the only dependency of the algorithms is the C++17 standard
library (FFTW for the STFT, libsndfile for the example program).

Usage in C++
------------

```cpp
#include "piva_online/auxiva_iss_online.hpp"

piva_online::AuxIVAISSOnline bss(n_freq, n_chan, /*n_iter=*/3, /*alpha=*/0.96);

// for every incoming STFT frame, x and y are (n_freq, n_chan) row-major
bss.process_frame(x, y);             // demixed frame
bss.project_back_frame(y, y_out, 0); // optional scale fix to mic 0
```

Differences from the Python versions
------------------------------------

* `process_frame` gives exactly the same output as the Python code with
  `proj_back=False` (checked by `tests/validate_vs_python.py`).
* Python applies the batch `project_back` on the whole output at the end.
  That is not causal, so `project_back_frame` instead uses the minimum
  distortion principle with the current demixing matrix
  (`y_k <- (W^-1)[ref, k] y_k`). This needs a small solve per bin, also for ISS.
* Only the determined case (`n_src == n_chan`) is supported, as in Python.

Build
-----

Using the `piva` conda environment:

    conda activate piva
    cmake -S . -B build
    cmake --build build -j

Run
---

    ./build/separate_online <input.wav> <output_prefix> <ip|iss> \
        [n_fft=4096] [hop=n_fft/2] [n_iter=3] [alpha=0.96] [model=laplace]

Writes `<output_prefix>_src<k>.wav` for each source.

Check against Python
--------------------

    python tests/validate_vs_python.py
