import os
import numpy as np
import matplotlib.pyplot as plt

source_values = [100, 500, 1000, 1500]
colors = ['-r', '-b', '-m', '-g']

fig, ax = plt.subplots()
for source, color in zip(source_values, colors):
    r, T = np.loadtxt(f'{source}.000000.dat', unpack=True)
    ax.plot(r, T, color, label=r'$S = {}$'.format(source))

ax.set_xlabel(r'$r$')
ax.set_ylabel(r'$T$')
ax.legend()
ax.set_title(r'Plot $T(r)$ for different sources')

plt.savefig('plot.png', dpi=500)
plt.clf()

