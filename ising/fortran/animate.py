import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

# Load all frames
data = [np.loadtxt(f"snapshot{i:03d}.txt") for i in range(100)]  # Assumes files are named 000.npy to 099.npy

# Setup plot
fig, ax = plt.subplots()
im = ax.imshow(data[0], cmap='gray', vmin=-1, vmax=1)

def update(frame):
    im.set_array(data[frame])
    return [im]

ani = animation.FuncAnimation(fig, update, frames=len(data), interval=50, blit=True)

plt.show()

# Optional: save as mp4
# ani.save('animation.mp4', fps=20)
