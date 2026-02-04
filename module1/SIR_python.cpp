#include <vector>
#include <cmath>
#include <stdexcept>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

// Function to take a step in the SIR model
// state: vector of S, I, R
// beta: infection rate
// gamma: recovery rate
// dt: time step
std::vector<double> take_step(const std::vector<double>& state, double beta, double gamma, double dt) {
    double S = state[0];
    double I = state[1];
    double R = state[2];

    double N = S + I + R;

    double dSdt = -beta * S * I / N;
    double dIdt =  beta * S * I / N - gamma * I;
    double dRdt =  gamma * I;

    double S_new = S + dSdt * dt;
    double I_new = I + dIdt * dt;
    double R_new = R + dRdt * dt;

    return {S_new, I_new, R_new};
}

// Function simulating num_steps of the SIR model, saving the state every return_every steps and returning the results
// Returns an array with columns: [t, S, I, R]
pybind11::array_t<double> integrate_system(double S0, double I0, double R0,
                                           double beta, double gamma, double dt,
                                           int num_steps, int return_every) {
    if (num_steps < 0) {
        throw pybind11::value_error("num_steps must be >= 0");
    }
    if (return_every <= 0) {
        throw pybind11::value_error("return_every must be >= 1");
    }
    if (dt <= 0.0) {
        throw pybind11::value_error("dt must be > 0");
    }

    std::vector<double> state = {S0, I0, R0};

    // how many rows we will store (include initial state at step=0)
    int num_records = (num_steps / return_every) + 1;

    // numpy array: shape (num_records, 4) -> [t, S, I, R]
    pybind11::array_t<double> out({num_records, 4});
    auto out_m = out.mutable_unchecked<2>();

    int rec = 0;

    // store initial
    out_m(rec, 0) = 0.0;
    out_m(rec, 1) = state[0];
    out_m(rec, 2) = state[1];
    out_m(rec, 3) = state[2];
    rec += 1;

    // integrate
    for (int step = 1; step <= num_steps; step++) {
        state = take_step(state, beta, gamma, dt);

        if (step % return_every == 0) {
            double t = step * dt;
            out_m(rec, 0) = t;
            out_m(rec, 1) = state[0];
            out_m(rec, 2) = state[1];
            out_m(rec, 3) = state[2];
            rec += 1;
        }
    }

    return out;
}

PYBIND11_MODULE(SIR_python, m) {
    m.doc() = "Python bindings for a simple SIR model (Forward Euler)";

    m.def("integrate_system", &integrate_system,
          pybind11::arg("S0"), pybind11::arg("I0"), pybind11::arg("R0"),
          pybind11::arg("beta"), pybind11::arg("gamma"), pybind11::arg("dt"),
          pybind11::arg("num_steps"), pybind11::arg("return_every"));
}
