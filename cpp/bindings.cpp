#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "backend.hpp"

namespace py = pybind11;
using namespace mdn;

PYBIND11_MODULE(moldyn_core, m) {
    py::class_<FeatureSet>(m, "FeatureSet")
        .def_readonly("distances", &FeatureSet::distances)
        .def_readonly("hbond_scores", &FeatureSet::hbond_scores)
        .def_readonly("pocket_center", &FeatureSet::pocket_center);

    m.def("compute_features", &compute_features,
          py::arg("coords_xyz"), py::arg("n_atoms"), py::arg("hbond_cutoff_angstrom")=3.5);

    m.def("run_mini_md", &run_mini_md, py::arg("pdb_path"), py::arg("steps")=500, py::arg("tempK")=300.0);
}
