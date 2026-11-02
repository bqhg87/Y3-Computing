import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from constants import *

def create_orbital_plot(system):
  # lock roll (by enforcing azimuth / elevation only rotation style)
  mpl.rcParams["axes3d.mouserotationstyle"] = "azel"

  fig = plt.figure(figsize=(9, 7), facecolor="black")
  ax = fig.add_subplot(projection="3d")

  ax.set_facecolor("black")

  if SHOW_GRID:
    ax.xaxis.set_pane_color((0.1, 0.1, 0.1, 0.5))
    ax.yaxis.set_pane_color((0.1, 0.1, 0.1, 0.5))
    ax.zaxis.set_pane_color((0.1, 0.1, 0.1, 0.5))
    ax.xaxis._axinfo["grid"]["color"] = (0.3, 0.3, 0.3, 0.5)
    ax.yaxis._axinfo["grid"]["color"] = (0.3, 0.3, 0.3, 0.5)
    ax.zaxis._axinfo["grid"]["color"] = (0.3, 0.3, 0.3, 0.5)
  else:
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.grid(False)

  # plot all bodies from the dictionary
  scatters = {}
  for name, data in system.items():
    pos = data["pos"]
    scatters[name] = ax.scatter(
        pos[0], pos[1], pos[2],
        color=data.get("color", "blue"),
        s=data.get("size", 50),
        marker=data.get("marker", "o"),
        label=name,
    )

  # labels and legend
  if SHOW_AXIS_LABELS:
    ax.set_xlabel("X (m)", color="white")
    ax.set_ylabel("Y (m)", color="white")
    ax.set_zlabel("Z (m)", color="white")
  if SHOW_TICK_LABELS:
    ax.tick_params(colors="white", labelcolor="white")
  ax.view_init(elev=30, azim=-60, roll=0)
  if SHOW_LEGEND:
    legend = ax.legend(loc="upper right", facecolor=(0.1, 0.1, 0.1, 0.5), edgecolor="none", labelcolor="white", prop={"family": "monospace", "size": 9})


  plt.tight_layout()
  plt.ion()
  plt.show(block=False)

  return fig, ax, scatters