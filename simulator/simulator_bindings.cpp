
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "FIELD_3D.h"
#include "PARTICLE_SYSTEM.h"

namespace py = pybind11;

class Simulator {
public:
    FIELD_3D field;
    PARTICLE_SYSTEM particles;
    int grid_size;

    Simulator(int grid_size) : grid_size(grid_size) {
        field = FIELD_3D(grid_size);
        particles = PARTICLE_SYSTEM(grid_size);
    }

    void step(int action) {
        particles.update(action);
        field.apply_suction(particles);
    }

    py::array_t<float> get_state() {
        std::vector<float> state(grid_size * grid_size);
        for (int i = 0; i < grid_size; ++i) {
            for (int j = 0; j < grid_size; ++j) {
                state[i * grid_size + j] = field.get_value(i, j);
            }
        }
        return py::array_t<float>(state.size(), state.data());
    }

    void reset() {
        particles.reset();
        field.reset();
    }
};

PYBIND11_MODULE(simulator, m) {
    py::class_<Simulator>(m, "Simulator")
        .def(py::init<int>())
        .def("step", &Simulator::step)
        .def("get_state", &Simulator::get_state)
        .def("reset", &Simulator::reset);
}
