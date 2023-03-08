import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Load data from file
x, y, sol = np.loadtxt('sol.dat', delimiter=' ', skiprows=1, unpack=True)

# Determine the dimensions of the grid
nx = len(np.unique(x))
ny = len(np.unique(y))

# Reshape the data into 2D arrays
X = x.reshape((nx,ny))
Y = y.reshape((nx,ny))
Z = sol.reshape((nx,ny))

# Create a figure and axis object
fig, ax = plt.subplots(figsize=(6, 5))

# Plot the 2D contour
levels = np.linspace(np.min(Z), np.max(Z), 10)
im = ax.contourf(X, Y, Z, levels=levels, cmap='inferno')

# Add a colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cb = plt.colorbar(im, cax=cax)
cb.ax.tick_params(labelsize=10)

# Set axis labels and title
ax.set_xlabel('x', fontsize=14)
ax.set_ylabel('y', fontsize=14)
ax.set_title('Contour Plot', fontsize=14)

# Adjust plot layout
fig.tight_layout()

# Save the figure in png format at 300 DPI
plt.savefig('contour_plot.png', dpi=300)
