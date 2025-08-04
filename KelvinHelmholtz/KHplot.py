import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns

"""Data is arranged as  x0 y0
                        x1 y0
                        x2 y0
                        x3 y0
                        etc.

                        x1 y0
                        ...
                        """
source_data = (
    "/Users/tomhollingsworth/Documents/LSC Internship/KelvinHelmholtz/Plot_12.025506"
)
save_dir = "/Users/tomhollingsworth/Documents/LSC Internship/KelvinHelmholtz/T_12/"
os.makedirs(save_dir, exist_ok=True)


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
raw_data = pd.read_csv(
    source_data, sep=r"\s+", header=None, skip_blank_lines=True, names=column_names
)

num_x = len(np.unique(raw_data["x"]))
num_y = len(np.unique(raw_data["y"]))


def convert_to_2D(row, num_xvalues, num_yvalues):
    array = row[:].to_numpy().reshape((num_yvalues, num_xvalues))
    return array


for column_name in raw_data.columns[2:]:
    heat_data = convert_to_2D(raw_data[column_name], num_x, num_y)

    plt.figure(figsize=(12, 7))
    ax = sns.heatmap(
        heat_data,
        cmap="magma",
    )
    plt.title(f"{column_name}, t=12")
    ax.set_xticks([])
    ax.set_yticks([])
    plt.savefig(os.path.join(save_dir, f"{column_name}"))
    print(f"Saved file {column_name}")
