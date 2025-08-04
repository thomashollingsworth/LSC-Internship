import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import re
import matplotlib.animation as animation

"""Data is arranged as  x0 y0
                        x1 y0
                        x2 y0
                        x3 y0
                        etc.

                        x1 y0
                        ...
                        """

source_dir = "/Users/tomhollingsworth/Documents/LSC Internship/KelvinHelmholtz/"
save_name = (
    "/Users/tomhollingsworth/Documents/LSC Internship/KelvinHelmholtz/BxAnimation"
)

files = [f for f in os.listdir(source_dir) if f.startswith("Plot_")]


def extract_time(filename):
    match = re.search(r"Plot_([0-9.]+)", filename)
    return float(match.group(1)) if match else float("inf")


# Sort files by extracted time
sorted_files = sorted(files, key=extract_time)

column_names = [
    "x",
    "y",
    "density",
    "v_x",
    "v_y",
    "v_z",
    "pressure",
    "B_x",
    "B_y",
    "B_z",
    "psi",
]

Bx_arrays = []


def convert_to_2D(row, num_xvalues, num_yvalues):
    array = row[:].to_numpy().reshape((num_yvalues, num_xvalues))
    return array


for file in sorted_files:
    raw_data = pd.read_csv(
        os.path.join(
            "/Users/tomhollingsworth/Documents/LSC Internship/KelvinHelmholtz/", file
        ),
        sep=r"\s+",
        header=None,
        skip_blank_lines=True,
        names=column_names,
    )

    num_x = len(np.unique(raw_data["x"]))
    num_y = len(np.unique(raw_data["y"]))

    Bx_data = convert_to_2D(raw_data["B_x"], num_x, num_y)
    Bx_arrays.append(Bx_data)

fig, ax = plt.subplots(figsize=(8, 6))


def init():
    plt.figure(figsize=(12, 7))
    sns.heatmap(Bx_data, cmap="magma", cbar=False, ax=ax)
    plt.title(f"B_x")
    ax.set_xticks([])
    ax.set_yticks([])
    return []


def update(frame):
    ax.clear()
    sns.heatmap(Bx_arrays[frame], cmap="magma", cbar=False, ax=ax)
    ax.set_title(f"Time step: {frame}")
    ax.set_xticks([])
    ax.set_yticks([])
    return []


ani = animation.FuncAnimation(
    fig, update, frames=len(Bx_arrays), init_func=init, blit=False, repeat=False
)

plt.tight_layout()
ani.save(f"{save_name}.mp4", writer="ffmpeg", fps=5)
