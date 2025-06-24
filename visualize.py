import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import imageio
import os
import glob
import re

# --- Configuration ---
DATA_DIR = "data"
OUTPUT_GIF = "simulation_output.gif"
FPS = 10 # Frames per second for the GIF

def get_simulation_time_from_filename(filename, dt, output_data_steps):
    """
    Extracts frame number from filename and calculates simulation time.
    Assumes filename format like 'frame_XXXX.dat'.
    """
    match = re.search(r"frame_(\d+)\.dat", os.path.basename(filename))
    if match:
        frame_num = int(match.group(1))
        # Time = frame_index * (dt * output_data_steps)
        # This assumes frame 0 is t=0, frame 1 is t=output_data_steps*dt, etc.
        return frame_num * dt * output_data_steps
    return None

def find_global_z_min_max(data_files):
    """
    Finds the global minimum and maximum z-values across all data files
    to ensure consistent Z-axis scaling in the animation.
    """
    global_min_z = float('inf')
    global_max_z = float('-inf')

    if not data_files:
        return 0, 0 # Default if no files

    for file_path in data_files:
        try:
            # Read data: x, y, height
            df = pd.read_csv(file_path, delim_whitespace=True, header=None, names=['x', 'y', 'height'])
            if not df.empty:
                global_min_z = min(global_min_z, df['height'].min())
                global_max_z = max(global_max_z, df['height'].max())
        except pd.errors.EmptyDataError:
            print(f"Warning: Data file {file_path} is empty. Skipping for z-range calculation.")
        except Exception as e:
            print(f"Error reading or processing {file_path} for z-range: {e}")
            continue

    if global_min_z == float('inf') or global_max_z == float('-inf'): # handles case where all files were empty
        return -1, 1 # Default range if no valid data found
    return global_min_z, global_max_z


def create_frame(data_file_path, frame_num, total_frames, sim_time, z_min, z_max):
    """
    Generates a single frame for the animation.
    - data_file_path: Path to the .dat file for this frame.
    - frame_num: Current frame number (for progress).
    - total_frames: Total number of frames (for progress).
    - sim_time: Simulation time for this frame.
    - z_min, z_max: Global min/max for Z-axis scaling.
    """
    print(f"Processing frame {frame_num + 1}/{total_frames} (Time: {sim_time:.2f}s) from {data_file_path}")

    try:
        # Read data: x, y, height
        df = pd.read_csv(data_file_path, delim_whitespace=True, header=None, names=['x', 'y', 'height'])
        if df.empty:
            print(f"Warning: Data file {data_file_path} is empty. Skipping frame.")
            return None # Skip if no data

        # Create a 2D grid for plotting
        # Determine grid dimensions (assuming they are consistent)
        grid_x_max = df['x'].max()
        grid_y_max = df['y'].max()

        # Create pivot table for Z values (height)
        # This might be slow for very large grids directly.
        # Consider if the C++ output can be made more grid-like if performance is an issue.
        # For now, pivot_table is convenient.
        grid_z = df.pivot_table(index='y', columns='x', values='height').values

        # Create X, Y coordinate matrices
        x_coords = np.arange(grid_z.shape[1]) # df['x'].unique()
        y_coords = np.arange(grid_z.shape[0]) # df['y'].unique()
        X, Y = np.meshgrid(x_coords, y_coords)

    except pd.errors.EmptyDataError:
        print(f"Error: The data file {data_file_path} is empty.")
        return None
    except Exception as e:
        print(f"Error processing data from {data_file_path}: {e}")
        return None

    # --- Plotting ---
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Surface plot
    surf = ax.plot_surface(X, Y, grid_z, cmap='viridis', edgecolor='none')

    # Set labels and title
    ax.set_xlabel("X coordinate")
    ax.set_ylabel("Y coordinate")
    ax.set_zlabel("Amplitude (Height)")
    ax.set_title(f"Seismic Wave Propagation - Time: {sim_time:.2f} s")

    # Set consistent Z-axis limits
    ax.set_zlim(z_min, z_max)

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5, label="Amplitude")

    # Save plot to a buffer
    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    plt.close(fig) # Close the figure to free memory
    return image


def main():
    """
    Main function to generate the GIF.
    """
    print("Starting visualization script...")

    # 1. Check if data directory exists and has files
    if not os.path.isdir(DATA_DIR):
        print(f"Error: Data directory '{DATA_DIR}' not found.")
        print("Please ensure the C++ simulation has been run and generated data.")
        return

    data_files = sorted(glob.glob(os.path.join(DATA_DIR, "frame_*.dat")))

    if not data_files:
        print(f"Error: No data files (frame_*.dat) found in '{DATA_DIR}'.")
        print("Execute first the C++ simulation to generate data.")
        return

    print(f"Found {len(data_files)} data files.")

    # --- Read config.json to get dt and output_data_steps for accurate time display ---
    # This is a simplified way to get these values.
    # A more robust way would be to use a proper JSON parser.
    dt_sim = 0.01 # Default, ideally read from config.json
    output_steps_sim = 5 # Default, ideally read from config.json
    try:
        import json
        with open("config.json", 'r') as f:
            config_data = json.load(f)
            dt_sim = config_data.get("time_step_dt", dt_sim)
            output_steps_sim = config_data.get("output_data_steps", output_steps_sim)
        print(f"Loaded dt={dt_sim}, output_data_steps={output_steps_sim} from config.json for time calculation.")
    except FileNotFoundError:
        print("Warning: config.json not found. Using default dt and output_data_steps for time calculation.")
    except Exception as e:
        print(f"Warning: Could not read config.json ({e}). Using default dt and output_data_steps.")


    # 2. Determine global Z min/max for consistent scaling
    print("Calculating global Z-axis range...")
    z_min, z_max = find_global_z_min_max(data_files)
    if z_min == z_max : # If all values are flat or only one file with one value
        z_min -= 0.5
        z_max += 0.5 # Add some padding if flat
    print(f"Global Z-axis range: [{z_min:.2f}, {z_max:.2f}]")


    # 3. Generate frames for GIF
    images = []
    total_files = len(data_files)
    for i, file_path in enumerate(data_files):
        sim_time = get_simulation_time_from_filename(file_path, dt_sim, output_steps_sim)
        if sim_time is None: # Fallback if regex fails
            sim_time = i * dt_sim * output_steps_sim # Approximate time
            print(f"Warning: Could not parse frame number from {file_path}. Approximating time to {sim_time:.2f}s.")

        frame_image = create_frame(file_path, i, total_files, sim_time, z_min, z_max)
        if frame_image is not None:
            images.append(frame_image)

    # 4. Save GIF
    if images:
        print(f"Saving GIF to {OUTPUT_GIF} with {len(images)} frames at {FPS} FPS...")
        try:
            imageio.mimsave(OUTPUT_GIF, images, fps=FPS, loop=0) # loop=0 means infinite loop
            print(f"GIF saved successfully: {OUTPUT_GIF}")
        except Exception as e:
            print(f"Error saving GIF: {e}")
            print("You might need to install the imageio ffmpeg plugin: pip install imageio[ffmpeg]")
    else:
        print("No frames were generated. GIF not created.")


if __name__ == "__main__":
    # // Entry point for the script
    main()

# TODO:
# - Consider more robust parsing of data files, especially if they can be very large.
#   (e.g., iterative processing if memory becomes an issue)
# - Add more customization options (e.g., colormap, viewing angle) via command-line arguments.
# - The current method of reading dt and output_data_steps from config.json is basic.
#   A shared configuration library or passing these values during C++ execution could be more robust.
# - Ensure grid dimensions are correctly inferred or read if they can vary.
#   The current script assumes x and y in the .dat files are 0-indexed and continuous.
```
