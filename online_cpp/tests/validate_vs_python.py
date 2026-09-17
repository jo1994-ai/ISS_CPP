"""
Checks the C++ online AuxIVA-IP / AuxIVA-ISS against the Python reference
implementations (piva/auxiva_ip_online.py, piva/auxiva_iss_online.py).

Both are run on the same STFT of a simulated 2- and 3-source room mixture,
without projection back, and the outputs must agree to round-off. Also
reports the speed-up of C++ over Python.

Run from online_cpp/ after building:
    python tests/validate_vs_python.py
"""
import os
import sys
import time

import numpy as np
import pyroomacoustics as pra

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(HERE, "..", "build"))
sys.path.insert(0, ROOT)

import piva_online_cpp  # noqa: E402
from piva.auxiva_ip_online import auxiva_ip_online  # noqa: E402
from piva.auxiva_iss_online import auxiva_iss_online  # noqa: E402
from utils import make_room, samples  # noqa: E402

fs = 16000
framesize = 4096
hop = framesize // 2

all_ok = True

for n_src in [2, 3]:
    np.random.seed(0)
    room = make_room(
        [10.0, 7.5, 3.2],
        n_src,
        [samples[i][1] for i in range(n_src)],
        fs=fs,
        max_order=17,
        absorption=0.35,
        array_radius=0.05,
    )
    premix = room.simulate(return_premix=True)
    premix /= np.std(premix[:, 0, :], axis=1)[:, None, None]
    mix = np.sum(premix, axis=0) + 10 ** (-15 / 10) * n_src * np.random.randn(
        *premix.shape[1:]
    )
    X = pra.transform.analysis(mix.T, framesize, hop, win=np.hamming(framesize))
    X = X.astype(np.complex128)

    for name, py_fn, cpp_fn in [
        ("IP ", auxiva_ip_online, piva_online_cpp.auxiva_ip_online),
        ("ISS", auxiva_iss_online, piva_online_cpp.auxiva_iss_online),
    ]:
        for model in ["laplace", "gauss"]:
            tic = time.perf_counter()
            Y_py = py_fn(X, n_iter=3, alpha=0.96, model=model, proj_back=False)
            t_py = time.perf_counter() - tic

            tic = time.perf_counter()
            Y_cpp = cpp_fn(X, n_iter=3, alpha=0.96, model=model)
            t_cpp = time.perf_counter() - tic

            err = np.max(np.abs(Y_py - Y_cpp)) / np.max(np.abs(Y_py))
            ok = err < 1e-8
            all_ok &= ok
            print(
                f"n_src={n_src} {name} {model:>7}: rel. max error {err:.1e} "
                f"[{'OK' if ok else 'MISMATCH'}]  python {t_py:6.3f} s  "
                f"C++ {t_cpp:6.3f} s  ({t_py / t_cpp:5.1f}x)"
            )

print("\nAll outputs match." if all_ok else "\nSome outputs do NOT match.")
sys.exit(0 if all_ok else 1)
