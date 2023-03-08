import os
import numpy as np
import matplotlib.pyplot as plt

# Set the directory containing the data files
data_dir = 'data'

# Get a list of the data files
data_files = sorted(os.listdir(data_dir))

# Set plot settings
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['lines.linewidth'] = 2.0

# Loop over the data files and plot each one
for i, data_file in enumerate(data_files):
    # Load the data from the file
    data = np.loadtxt(os.path.join(data_dir, data_file))

    # Extract the time, x, u, and uexact data
    time = data[0, 0]
    x = data[:, 1]
    u = data[:, 2]
    uexact = data[:, 3]

    # Create a figure with a single subplot
    fig, ax = plt.subplots()

    # Plot the numerical and exact solutions
    ax.plot(x, u, label='Numerical')
    ax.plot(x, uexact, label='Exact')

    # Set the plot title and labels
    ax.set_title('Time = {:.2f}'.format(time))
    ax.set_xlabel('x')
    ax.set_ylabel('u')
    ax.legend()
    ax.set_ylim(-1, 1)

    # Save the figure as a high-quality PNG image
    fig.savefig('plot/plot_{:05d}.png'.format(i), dpi=600, bbox_inches='tight')

    # Close the figure to free up memory
    plt.close(fig)
