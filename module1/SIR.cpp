#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <string>
#include <filesystem>

// state: vector of S, I, R
std::vector<double> take_step(std::vector<double> state, double beta, double gamma, double dt) {
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

int main(int argc, char* argv[]) {
    // parameters
    double beta = 0.2;
    double gamma = pow(10, -1);   // 0.1
    int N = 1000;

    int population = N;

    double infected = 1.0;
    double susceptible = population - infected;
    double immune = 0.0;

    // state is S, I, R
    std::vector<double> state = {susceptible, infected, immune};

    // timestep
    double dt = 0.1;
    if (argc >= 2) {
        dt = std::stod(argv[1]);
    }

    // fixed simulation end time for fair comparison
    double t_end = 200.0;

    // write output every 1 day (must be divisible by dt)
    double write_dt = 1.0;

    double ratio = write_dt / dt;
    long long write_every = llround(ratio);

    double tol = 1e-9 * std::max(1.0, std::fabs(ratio));
    if (std::fabs(ratio - static_cast<double>(write_every)) > tol) {
        std::cerr << "ERROR: write_dt/dt is not an integer (within tolerance).\n";
        std::cerr << "write_dt = " << write_dt << ", dt = " << dt << ", ratio = " << ratio << "\n";
        return 1;
    }

    // filename based on dt
    std::ostringstream base;
    base << "filename_dt" << std::fixed << std::setprecision(6) << dt;
    std::string base_name = base.str();
    for (char& c : base_name) if (c == '.') c = 'p';

    std::filesystem::create_directories("data");
    std::string filename = "data/" + base_name + ".csv";

    std::ofstream outFile(filename);
    outFile << "Time,Infected,Susceptibile,Recovered\n";

    long long step = 0;

    while (true) {
        double current_time = step * dt;
        if (current_time > t_end + 1e-12) break;

        if (step % write_every == 0) {
            outFile << std::fixed << std::setprecision(12) << current_time << ","
                    << state[1] << ","   // Infected
                    << state[0] << ","   // Susceptible
                    << state[2] << "\n"; // Recovered
        }

        state = take_step(state, beta, gamma, dt);
        step += 1;
    }

    outFile.close();
    return 0;
}
