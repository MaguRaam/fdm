import numpy as np
import matplotlib.pyplot as plt

x, y = np.loadtxt('sol.dat', unpack=True)

plt.plot(x, y, color='red')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.savefig('plot.png', dpi=500)
plt.clf()
