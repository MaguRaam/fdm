import numpy as np
import matplotlib.pyplot as plt
from math import log

def plot_convergence(n, error, filename):
    fig, ax = plt.subplots()
    ax.loglog(n, error, '-bo', label='average error')
    ax.set_xlabel(r'$N$')
    ax.set_ylabel(r'$E$')
    ax.set_title(r'average error at t = 0.4 and rd = 1/6')
    ax.legend()
    fig.savefig(filename, dpi=600, bbox_inches='tight')


def convergence_rate(error):
    return [(log(error[i-1]) - log(error[i])) / log(2) for i in range(1, np.size(error))]

# load n and error
n, error = np.loadtxt("error.dat", unpack=True)

# plot convergence
plot_convergence(n, error, "error.png")

# calculate and print convergence rates
rates = convergence_rate(error)
print("Convergence rates:", rates)
