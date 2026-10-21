import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def create_orbital_plot(system):
  # lock roll (by enforcing azimuth / elevation only rotation style)
  mpl.rcParams["axes3d.mouserotationstyle"] = "azel"

  fig = plt.figure(figsize=(9, 7))
  ax = fig.add_subplot(projection="3d")

  # plot all bodies from the dictionary
  for name, data in system.items():
    pos = data["pos"]
    ax.scatter(
        pos[0],
        pos[1],
        pos[2],
        color=data.get("color", "blue"),
        s=data.get("size", 50),
        marker=data.get("marker", "o"),
        label=name,
    )

  # labels and legend
  ax.set_xlabel("X (m)")
  ax.set_ylabel("Y (m)")
  ax.set_zlabel("Z (m)")
  ax.legend(loc="upper right")
  ax.view_init(elev=30, azim=-60, roll=0)

  # origin-centered scroll-to-zoom handler
  def on_scroll(event):
    if event.inaxes == ax:
      scale = 0.9 if event.button == "up" else 1.1
      for get_lim, set_lim in [
          (ax.get_xlim, ax.set_xlim),
          (ax.get_ylim, ax.set_ylim),
          (ax.get_zlim, ax.set_zlim),
      ]:
        span = (get_lim()[1] - get_lim()[0]) * scale / 2.0
        set_lim([-span, span])
      fig.canvas.draw_idle()

  fig.canvas.mpl_connect("scroll_event", on_scroll)

  plt.tight_layout()
  plt.show()