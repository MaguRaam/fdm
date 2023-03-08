#plot x(t) for different damping coefficients:

import os
import numpy as np
import matplotlib.pyplot as plt

# load files for different damping coefficient c:
t, x5 = np.loadtxt('5.000000.dat', unpack=True)
t, x40 = np.loadtxt('40.000000.dat', unpack=True)
t, x200 = np.loadtxt('200.000000.dat', unpack=True)

# plot for c = 5, 40, 200 in the same graph
plt.plot(t, x5, '-r', label='$c = 5$, underdamped')
plt.plot(t, x40, '-b', label='$c = 40$, critically damped')
plt.plot(t, x200, '-m', label='$c = 200$, overdamped')
plt.xlabel(r'$t$')
plt.ylabel(r'$x$')
plt.legend()
plt.title(r'Plot $x(t)$ for different damping coefficients')
plt.savefig('damping_coeff.png', dpi=500)
plt.clf()
