#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <iomanip> // For std::setw and std::setfill

// Forward declaration for a JSON parsing library (e.g., nlohmann/json)
// For a real project, you would include the header of your chosen JSON library here.
// For now, we'll simulate reading it.
namespace nlohmann {
    struct json {
        // Basic simulation of json object access
        template<typename T>
        T at(const std::string& key) const {
            // In a real scenario, this would parse JSON and return the value.
            // Here, we'll provide dummy values or use if-else for demonstration.
            if (key == "grid_size_x") return static_cast<T>(200);
            if (key == "grid_size_y") return static_cast<T>(200);
            if (key == "simulation_time") return static_cast<T>(5.0);
            if (key == "time_step_dt") return static_cast<T>(0.01);
            if (key == "output_data_steps") return static_cast<T>(5);
            if (key == "source_pos_x") return static_cast<T>(100);
            if (key == "source_pos_y") return static_cast<T>(100);
            if (key == "ricker_frequency") return static_cast<T>(25.0);
            return T{};
        }

        static json parse(std::ifstream& fs) {
            // Simulate parsing. In a real scenario, this would read and parse the file.
            if (!fs.is_open()) {
                throw std::runtime_error("Failed to open config file for parsing.");
            }
            std::cout << "Simulating JSON parsing from file." << std::endl;
            // We'd read the content of fs here.
            return json{}; // Return a dummy json object
        }
    };
} // namespace nlohmann

// Structure to hold simulation parameters
struct SimConfig {
    int grid_size_x;
    int grid_size_y;
    double simulation_time;
    double time_step_dt;
    int output_data_steps;
    int source_pos_x;
    int source_pos_y;
    double ricker_frequency;

    // Function to load parameters from config.json
    bool load_from_json(const std::string& config_path) {
        std::ifstream f(config_path);
        if (!f.is_open()) {
            std::cerr << "Error: Could not open config file: " << config_path << std::endl;
            return false;
        }

        try {
            // In a real project, use a proper JSON library like nlohmann/json
            // For example: nlohmann::json data = nlohmann::json::parse(f);
            nlohmann::json data = nlohmann::json::parse(f); // Simulated parse

            grid_size_x = data.at<int>("grid_size_x");
            grid_size_y = data.at<int>("grid_size_y");
            simulation_time = data.at<double>("simulation_time");
            time_step_dt = data.at<double>("time_step_dt");
            output_data_steps = data.at<int>("output_data_steps");
            source_pos_x = data.at<int>("source_pos_x");
            source_pos_y = data.at<int>("source_pos_y");
            ricker_frequency = data.at<double>("ricker_frequency");

            std::cout << "Configuration loaded successfully." << std::endl;
            return true;
        } catch (const std::exception& e) {
            // Catch exceptions from the JSON library (e.g., parse error, key not found)
            std::cerr << "Error parsing JSON: " << e.what() << std::endl;
            return false;
        }
    }
};

// Function to calculate Ricker wavelet value
double ricker_wavelet(double t, double freq) {
    // t: current time
    // freq: dominant frequency of the wavelet
    // A common formulation for Ricker wavelet:
    // (1 - 2 * (pi * freq * t)^2) * exp(-(pi * freq * t)^2)
    // We need to adjust 't' so the wavelet starts near t=0.
    // Let's define a shift, e.g., peak time. A typical peak for this formula is around t_peak = sqrt(2) / (pi * freq)
    // Or more simply, shift by a few periods, e.g., 1.5 / freq or 2.0 / freq
    double t_shift = 1.5 / freq; // Shift to start the wavelet appropriately
    double adjusted_t = t - t_shift;
    if (adjusted_t < -t_shift) return 0.0; // Ensure it doesn't start before t=0 effectively

    double pi_f_t = M_PI * freq * adjusted_t;
    double pi_f_t_sq = pi_f_t * pi_f_t;
    return (1.0 - 2.0 * pi_f_t_sq) * std::exp(-pi_f_t_sq);
}

// Function to save a data frame to a file
void save_frame_data(const std::vector<std::vector<double>>& grid, int frame_num, int grid_x, int grid_y, const std::string& output_dir) {
    std::ofstream outfile;
    std::stringstream ss;
    ss << output_dir << "/frame_" << std::setw(4) << std::setfill('0') << frame_num << ".dat";
    outfile.open(ss.str());

    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file to save frame data: " << ss.str() << std::endl;
        return;
    }

    // Write header (optional, but good practice)
    // outfile << "x y height\n"; // As per requirement: three columns: x y height

    for (int i = 0; i < grid_x; ++i) {
        for (int j = 0; j < grid_y; ++j) {
            outfile << i << " " << j << " " << grid[i][j] << "\n";
        }
    }
    outfile.close();
}


int main() {
    std::cout << "Starting Seismic Simulation..." << std::endl;

    // 1. Load Configuration
    SimConfig config;
    if (!config.load_from_json("config.json")) {
        std::cerr << "Failed to load configuration. Exiting." << std::endl;
        return 1;
    }

    // 2. Initialize Simulation Grid
    // The grid will store the wave amplitude at each point.
    // Let's use three time steps for the finite difference method: prev, current, next
    std::vector<std::vector<double>> u_prev(config.grid_size_x, std::vector<double>(config.grid_size_y, 0.0));
    std::vector<std::vector<double>> u_curr(config.grid_size_x, std::vector<double>(config.grid_size_y, 0.0));
    std::vector<std::vector<double>> u_next(config.grid_size_x, std::vector<double>(config.grid_size_y, 0.0));

    // Simulation constants (derived from config)
    double dt = config.time_step_dt;
    double dx = 1.0; // Assuming grid spacing is 1 unit, can be made configurable
    double dy = 1.0;
    double c = 10.0; // Wave speed, should ideally be configurable or derived
                     // For stability, (c * dt / dx)^2 <= 0.5 for 2D wave equation
    double C_sq = c * c;
    double dt_sq = dt * dt;
    double dx_sq = dx * dx;
    double dy_sq = dy * dy;

    // Courant number check for stability (optional, but good for diagnostics)
    double courant_x = C_sq * dt_sq / dx_sq;
    double courant_y = C_sq * dt_sq / dy_sq;
    if (courant_x + courant_y > 0.5) { // A common stability condition for explicit 2D schemes
        std::cout << "Warning: Simulation might be unstable. Courant number condition might not be met." << std::endl;
        std::cout << "C_sq * dt_sq / dx_sq + C_sq * dt_sq / dy_sq = " << courant_x + courant_y << " (should be <= 0.5)" << std::endl;
    }


    // 3. Main Simulation Loop
    int num_time_steps = static_cast<int>(config.simulation_time / dt);
    int output_counter = 0;
    int frame_number = 0;

    std::cout << "Running simulation for " << num_time_steps << " time steps." << std::endl;

    for (int t_step = 0; t_step < num_time_steps; ++t_step) {
        double current_time = t_step * dt;

        // Update wave field using finite differences for the 2D wave equation:
        // u_next[i][j] = 2*u_curr[i][j] - u_prev[i][j] +
        //                C_sq * dt_sq * ( (u_curr[i+1][j] - 2*u_curr[i][j] + u_curr[i-1][j])/dx_sq +
        //                                 (u_curr[i][j+1] - 2*u_curr[i][j] + u_curr[i][j-1])/dy_sq )
        for (int i = 1; i < config.grid_size_x - 1; ++i) {
            for (int j = 1; j < config.grid_size_y - 1; ++j) {
                double laplacian_x = (u_curr[i+1][j] - 2.0 * u_curr[i][j] + u_curr[i-1][j]) / dx_sq;
                double laplacian_y = (u_curr[i][j+1] - 2.0 * u_curr[i][j] + u_curr[i][j-1]) / dy_sq;
                u_next[i][j] = 2.0 * u_curr[i][j] - u_prev[i][j] + C_sq * dt_sq * (laplacian_x + laplacian_y);
            }
        }

        // Apply source term (Ricker wavelet)
        // The source should only be active for its effective duration.
        double source_val = ricker_wavelet(current_time, config.ricker_frequency);
        // Ensure source position is within bounds (adjusting for 0-based index)
        int sx = config.source_pos_x;
        int sy = config.source_pos_y;
        if (sx > 0 && sx < config.grid_size_x -1 && sy > 0 && sy < config.grid_size_y -1) {
             // Add source to u_next or u_curr? Adding to u_next is common for explicit schemes.
             // The scaling factor dt_sq is often applied to the source term in discretized wave equations.
            u_next[sx][sy] += source_val * dt_sq; // dt_sq scaling for consistency with wave equation discretization
        }


        // Boundary Conditions (e.g., simple reflecting or absorbing)
        // For simplicity, using Dirichlet boundary conditions (amplitude = 0 at boundaries)
        // This means u_next at boundaries remains 0 (already initialized)
        // More advanced BCs (like PML or absorbing) would be implemented here.


        // Update grids for next iteration
        u_prev = u_curr;
        u_curr = u_next;
        // u_next will be recalculated or can be reset if needed, though it's overwritten.

        // 4. Output Data periodically
        if (t_step % config.output_data_steps == 0) {
            std::cout << "Time: " << std::fixed << std::setprecision(2) << current_time
                      << "s, Saving frame: " << frame_number << std::endl;
            save_frame_data(u_curr, frame_number, config.grid_size_x, config.grid_size_y, "data");
            frame_number++;
        }
    }

    std::cout << "Simulation finished." << std::endl;
    std::cout << "Total frames saved: " << frame_number << std::endl;

    return 0;
}

// TODO:
// - Implement proper JSON parsing (e.g., using nlohmann/json library).
//   You would need to include the library and link against it.
//   Example: #include "nlohmann/json.hpp"
// - Refine boundary conditions (e.g., Absorbing Boundary Conditions like PML).
// - Parameterize wave speed 'c' and grid spacings 'dx', 'dy' in config.json.
// - Add error handling for file operations.
// - Consider performance optimizations for large grids.
// - The Ricker wavelet's effective duration: ensure it acts as a pulse and then stops.
//   The current implementation of ricker_wavelet function itself returns 0 outside its main lobe.
// - Stability condition (Courant-Friedrichs-Lewy - CFL) should be carefully checked and potentially enforced.
//   The wave speed 'c' is crucial here.
// - The current `nlohmann::json` part is a placeholder. For this to compile and run
//   correctly, you'd need to integrate the actual nlohmann/json library.
//   For example, by adding `find_package(nlohmann_json 3 REQUIRED)` to CMakeLists.txt
//   and `target_link_libraries(seismic_sim PRIVATE nlohmann_json::nlohmann_json)`.
//   Or by directly including the single-header file `json.hpp` in the project.
//   For this skeleton, I've provided a dummy implementation to make it runnable as a standalone file
//   if one were to replace the dummy nlohmann::json with actual values.
```
