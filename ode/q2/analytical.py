import math
import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(0, 2,1000)
y = [math.exp((x**3/3) - 1.1*x) for x in t]

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(t, y, '-r', label=r'$y = e^{(\frac{t^3}{3} - 1.1t)}$')
ax.set_xlabel('$t$', fontsize=14)
ax.set_ylabel('$y(t)$', fontsize=14)
ax.legend(fontsize=12)
ax.set_title(r'Solution to the ODE $dy/dt = yt^2 - 1.1y$', fontsize=16)
plt.tight_layout()
plt.savefig('ode_solution.png', dpi=300)
