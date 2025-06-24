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
FPS = 20 # Puedes ajustar esto como prefieras

def get_simulation_time_from_filename(filename, dt, output_data_steps):
    match = re.search(r"frame_(\d+)\.dat", os.path.basename(filename))
    if match:
        frame_num = int(match.group(1))
        return frame_num * dt * output_data_steps
    return None

def find_global_z_range_with_percentiles(data_files, low_perc=1, high_perc=99):
    all_heights = []
    print("Reading all data files to calculate value ranges...")
    for file_path in data_files:
        try:
            df = pd.read_csv(file_path, sep=r'\s+', header=None, names=['x', 'y', 'height'])
            if not df.empty:
                all_heights.extend(df['height'].tolist())
        except Exception as e:
            print(f"Warning: Could not process {file_path}: {e}")
            continue

    if not all_heights:
        return -1, 1, -1, 1

    all_heights_series = pd.Series(all_heights)
    
    # Rango para la escala de colores Y para el eje Z (basado en percentiles)
    vmin = all_heights_series.quantile(low_perc / 100.0)
    vmax = all_heights_series.quantile(high_perc / 100.0)
    
    # Asegurarse de que el rango no sea cero
    if vmin == vmax:
        vmin -= 0.5
        vmax += 0.5

    print(f"Using percentile range for Z-Axis and Color Scale ({low_perc}th-{high_perc}th): [{vmin:.4f}, {vmax:.4f}]")
    
    return vmin, vmax

def create_frame(data_file_path, frame_num, total_frames, sim_time, z_lim_min, z_lim_max, v_min, v_max):
    print(f"Processing frame {frame_num + 1}/{total_frames} (Time: {sim_time:.2f}s) from {data_file_path}")
    try:
        df = pd.read_csv(data_file_path, sep=r'\s+', header=None, names=['x', 'y', 'height'])
        if df.empty:
            return None
        grid_z = df.pivot_table(index='y', columns='x', values='height').values
        x_coords = np.arange(grid_z.shape[1])
        y_coords = np.arange(grid_z.shape[0])
        X, Y = np.meshgrid(x_coords, y_coords)
    except Exception as e:
        print(f"Error processing data from {data_file_path}: {e}")
        return None

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(X, Y, grid_z, cmap='viridis', edgecolor='none', vmin=v_min, vmax=v_max)

    ax.set_xlabel("X coordinate")
    ax.set_ylabel("Y coordinate")
    ax.set_zlabel("Amplitude (Height)")
    ax.set_title(f"Seismic Wave Propagation - Time: {sim_time:.2f} s")
    
    # Usamos el mismo rango para los límites del eje Z
    ax.set_zlim(z_lim_min, z_lim_max)

    fig.colorbar(surf, shrink=0.5, aspect=5, label="Amplitude")
    fig.canvas.draw()
    rgba_buffer = fig.canvas.buffer_rgba()
    image = np.asarray(rgba_buffer)[:, :, :3]
    plt.close(fig)
    return image

def main():
    print("Starting visualization script...")
    if not os.path.isdir(DATA_DIR):
        print(f"Error: Data directory '{DATA_DIR}' not found.")
        return
    data_files = sorted(glob.glob(os.path.join(DATA_DIR, "frame_*.dat")))
    if not data_files:
        print(f"Error: No data files found in '{DATA_DIR}'.")
        return

    print(f"Found {len(data_files)} data files.")
    dt_sim, output_steps_sim = 0.01, 5
    try:
        import json
        with open("config.json", 'r') as f:
            config_data = json.load(f)
            dt_sim = config_data.get("time_step_dt", dt_sim)
            output_steps_sim = config_data.get("output_data_steps", output_steps_sim)
        print(f"Loaded dt={dt_sim}, output_data_steps={output_steps_sim} from config.json.")
    except Exception as e:
        print(f"Warning: Could not read config.json ({e}). Using default values.")

    # Obtenemos el rango de percentiles
    vmin, vmax = find_global_z_range_with_percentiles(data_files)

    images = []
    total_files = len(data_files)
    for i, file_path in enumerate(data_files):
        sim_time = get_simulation_time_from_filename(file_path, dt_sim, output_steps_sim)
        if sim_time is None:
            sim_time = i * dt_sim * output_steps_sim
        
        # --- CAMBIO IMPORTANTE AQUÍ ---
        # Pasamos el mismo rango (vmin, vmax) a los límites del eje Y a los de los colores.
        frame_image = create_frame(file_path, i, total_files, sim_time, vmin, vmax, vmin, vmax)
        
        if frame_image is not None:
            images.append(frame_image)

    if images:
        print(f"Saving GIF to {OUTPUT_GIF} with {len(images)} frames at {FPS} FPS...")
        imageio.mimsave(OUTPUT_GIF, images, fps=FPS, loop=0)
        print(f"GIF saved successfully: {OUTPUT_GIF}")
    else:
        print("No frames were generated. GIF not created.")

if __name__ == "__main__":
    main()
