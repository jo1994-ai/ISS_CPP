// Python bindings for the online AuxIVA C++ implementations, with the same
// call signature as the Python versions so they can be compared directly.
//
//   import piva_online_cpp
//   Y = piva_online_cpp.auxiva_ip_online(X, n_iter=3, alpha=0.96, model="laplace")
//   Y = piva_online_cpp.auxiva_iss_online(X, n_iter=3, alpha=0.96, model="laplace")
//
// X is (n_frames, n_freq, n_chan) complex128; Y has the same shape and is
// the raw demixed output (no projection back), frame t computed causally.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "piva_online/auxiva_ip_online.hpp"
#include "piva_online/auxiva_iss_online.hpp"

namespace py = pybind11;
using namespace piva_online;

using carray = py::array_t<cplx, py::array::c_style | py::array::forcecast>;

template <class Algo>
carray run(carray X, size_t n_iter, double alpha, const std::string& model) {
  if (X.ndim() != 3) throw std::invalid_argument("X must be (n_frames, n_freq, n_chan)");
  size_t n_frames = X.shape(0), n_freq = X.shape(1), n_chan = X.shape(2);

  carray Y({n_frames, n_freq, n_chan});
  const cplx* x = X.data();
  cplx* y = Y.mutable_data();

  {
    py::gil_scoped_release release;
    Algo bss(n_freq, n_chan, n_iter, alpha, model_from_string(model));
    size_t stride = n_freq * n_chan;
    for (size_t t = 0; t < n_frames; t++) bss.process_frame(x + t * stride, y + t * stride);
  }
  return Y;
}

PYBIND11_MODULE(piva_online_cpp, m) {
  m.doc() = "Online AuxIVA-IP and AuxIVA-ISS (C++)";
  m.def("auxiva_ip_online", &run<AuxIVAIPOnline>, py::arg("X"), py::arg("n_iter") = 3,
        py::arg("alpha") = 0.96, py::arg("model") = "laplace");
  m.def("auxiva_iss_online", &run<AuxIVAISSOnline>, py::arg("X"), py::arg("n_iter") = 3,
        py::arg("alpha") = 0.96, py::arg("model") = "laplace");
}
