#include "json.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <iomanip>

using json = nlohmann::json;

struct SimConfig {
    int grid_size_x, grid_size_y;
    double simulation_time, time_step_dt;
    int output_data_steps;
    int source_pos_x, source_pos_y;
    double ricker_frequency, source_amplitude;
    double wave_speed_c;
    double grid_spacing_dx;
    double grid_spacing_dy;

    bool load_from_json(const std::string& config_path) {
        std::ifstream f(config_path);
        if (!f.is_open()) {
            std::cerr << "Error: Could not open config file: " << config_path << std::endl;
            return false;
        }
        try {
            json data = json::parse(f);
            grid_size_x = data.at("grid_size_x");
            grid_size_y = data.at("grid_size_y");
            simulation_time = data.at("simulation_time");
            time_step_dt = data.at("time_step_dt");
            output_data_steps = data.at("output_data_steps");
            source_pos_x = data.at("source_pos_x");
            source_pos_y = data.at("source_pos_y");
            ricker_frequency = data.at("ricker_frequency");
            source_amplitude = data.at("source_amplitude");
            wave_speed_c = data.at("wave_speed_c");
            grid_spacing_dx = data.at("grid_spacing_dx");
            grid_spacing_dy = data.at("grid_spacing_dy");
            std::cout << "Configuration loaded successfully from " << config_path << std::endl;
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Error processing JSON: " << e.what() << std::endl;
            return false;
        }
    }
};

void save_frame_data(const std::vector<std::vector<double>>& grid, int frame_num, int grid_x, int grid_y, const std::string& output_dir) {
    std::ofstream outfile;
    std::stringstream ss;
    ss << output_dir << "/frame_" << std::setw(4) << std::setfill('0') << frame_num << ".dat";
    outfile.open(ss.str());
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file to save frame data: " << ss.str() << std::endl;
        return;
    }
    for (int i = 0; i < grid_x; ++i) {
        for (int j = 0; j < grid_y; ++j) {
            outfile << i << " " << j << " " << grid[i][j] << "\n";
        }
    }
    outfile.close();
}


double ricker_wavelet(double t, double freq) {
    double t_shift = 1.5 / freq;
    double adjusted_t = t - t_shift;
    if (adjusted_t < -t_shift) return 0.0;
    double pi_f_t = M_PI * freq * adjusted_t;
    double pi_f_t_sq = pi_f_t * pi_f_t;
    return (1.0 - 2.0 * pi_f_t_sq) * std::exp(-pi_f_t_sq);
}


int main() {
    std::cout << "Starting Seismic Simulation..." << std::endl;
    SimConfig config;
    if (!config.load_from_json("config.json")) {
        std::cerr << "Failed to load configuration. Exiting." << std::endl;
        return 1;
    }

    std::vector<std::vector<double>> u_prev(config.grid_size_x, std::vector<double>(config.grid_size_y, 0.0));
    std::vector<std::vector<double>> u_curr(config.grid_size_x, std::vector<double>(config.grid_size_y, 0.0));
    std::vector<std::vector<double>> u_next(config.grid_size_x, std::vector<double>(config.grid_size_y, 0.0));

    double dt = config.time_step_dt;
    double dx = config.grid_spacing_dx;
    double dy = config.grid_spacing_dy;
    double c = config.wave_speed_c;
    
    double C_sq = c * c;
    double dt_sq = dt * dt;
    double dx_sq = dx * dx;
    double dy_sq = dy * dy;

    double courant_sum = C_sq * dt_sq * (1/dx_sq + 1/dy_sq);
    if (courant_sum > 1.0) {
        std::cout << "Warning: Simulation might be unstable. Courant condition not met." << std::endl;
        std::cout << "Value: " << courant_sum << " (should be <= 1.0)" << std::endl;
    } else {
        std::cout << "Courant condition OK. Value: " << courant_sum << std::endl;
    }

    int num_time_steps = static_cast<int>(config.simulation_time / dt);
    int frame_number = 0;
    std::cout << "Running simulation for " << num_time_steps << " time steps." << std::endl;

    for (int t_step = 0; t_step < num_time_steps; ++t_step) {
        double current_time = t_step * dt;
        for (int i = 1; i < config.grid_size_x - 1; ++i) {
            for (int j = 1; j < config.grid_size_y - 1; ++j) {
                double laplacian_x = (u_curr[i+1][j] - 2.0 * u_curr[i][j] + u_curr[i-1][j]) / dx_sq;
                double laplacian_y = (u_curr[i][j+1] - 2.0 * u_curr[i][j] + u_curr[i][j-1]) / dy_sq;
                u_next[i][j] = 2.0 * u_curr[i][j] - u_prev[i][j] + C_sq * dt_sq * (laplacian_x + laplacian_y);
            }
        }

        double source_val = ricker_wavelet(current_time, config.ricker_frequency);
        int sx = config.source_pos_x;
        int sy = config.source_pos_y;
        if (sx > 0 && sx < config.grid_size_x -1 && sy > 0 && sy < config.grid_size_y -1) {
             u_next[sx][sy] += config.source_amplitude * source_val;
        }

        u_prev = u_curr;
        u_curr = u_next;

        if (t_step % config.output_data_steps == 0) {
            // ----- ESPÍA 1: MENSAJE DENTRO DEL IF -----
            std::cout << "Saving frame " << frame_number << " at time step " << t_step << "..." << std::endl;
            save_frame_data(u_curr, frame_number++, config.grid_size_x, config.grid_size_y, "data");
        }
    }

    std::cout << "Simulation finished." << std::endl;
    // ----- ESPÍA 2: MENSAJE FINAL -----
    std::cout << "Total frames saved: " << frame_number << std::endl;
    return 0;
}
