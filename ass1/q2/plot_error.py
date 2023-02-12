# reads t vs y for euler, rk4 and ab2 schemes for different time steps and plots dt vs error in log-log plot:

import os
import numpy as np
import matplotlib.pyplot as plt
from math import *
import string

#--------------------functions-------------------------------------#

# analytical solution:


def f(t): return exp(-t)


# get scheme, dt, mean and max error from file:
def postprocess(filename, f):

    t, y = np.loadtxt(filename, unpack=True)
    F = np.vectorize(f)
    y_exact = F(t)
    error = np.absolute(y - y_exact)
    mean_error = np.mean(error)
    max_error = np.max(error)

    # get scheme and dt:
    filename = filename.replace('.dat', '')
    split = filename.split('_',)
    scheme = split[0]
    dt = float(split[1])

    return scheme, dt, mean_error, max_error

#print convergence rate:
def convergence_rate(error):
  error.sort(reverse=True)  
  for i in range(1,np.size(error)):
    covergence_rate = (log(error[i-1]) - log(error[i]))/log(2)
    print("order of convergence ~", covergence_rate)



#------------------------main program--------------------------------------#

# get scheme, dt, mean and max error from .dat files from current directory and append it to list:

# create list to store dt and error for each scheme:
euler_max, euler_mean, euler_dt = ([] for i in range(3))
rk4_max, rk4_mean, rk4_dt = ([] for i in range(3))
ab2_max, ab2_mean, ab2_dt = ([] for i in range(3))


# loop over .dat files in current directory
for file in os.listdir():
    if file.endswith(".dat"):

        # get file data
        scheme, dt, mean_error, max_error = postprocess(file, f)

        # append to lists
        if (scheme == 'euler'):
            euler_dt.append(dt)
            euler_max.append(max_error)
            euler_mean.append(mean_error)

        elif (scheme == 'rk4'):
            rk4_dt.append(dt)
            rk4_max.append(max_error)
            rk4_mean.append(mean_error)

        elif (scheme == 'ab2'):
            ab2_dt.append(dt)
            ab2_max.append(max_error)
            ab2_mean.append(mean_error)



# dt vs mean error plot for different schemes:

# create a first order accurate test error function C1dt
C1 = euler_mean[0]/euler_dt[0] + 0.01  # compute constant C1
test_error1 = [C1*dt_ for dt_ in euler_dt]  # compute C1dt

# create a fourth order accurate test error function C4dt^4
C4 = rk4_mean[0]/rk4_dt[0] + 0.001  # compute constant C4
test_error4 = [C4*dt_**4 for dt_ in rk4_dt]  # compute C4dt^4

# create a second order accurate test error function C2dt^2
C2 = ab2_mean[0]/ab2_dt[0] + 0.001  # compute constant C2
test_error2 = [C2*dt_**2 for dt_ in ab2_dt]  # compute C2dt^2



plt.loglog(euler_dt, euler_mean, '-r', label='Euler')
plt.loglog(euler_dt, test_error1, linestyle='--', label='first-order slope')
plt.loglog(rk4_dt, rk4_mean, '-b', label='RK4')
plt.loglog(rk4_dt, test_error4, linestyle='--', label='fourth-order slope')
plt.loglog(ab2_dt, ab2_mean, '-m', label='AB2')
plt.loglog(ab2_dt, test_error2, linestyle='--', label='second-order slope')
plt.xlabel(r'$dt$')
plt.ylabel(r'$|E(t)|_{avg}$')
plt.title(r'dt vs average error')
plt.legend()
plt.savefig("average_error.png", dpi=600)
plt.clf()


# dt vs max error plot for different schemes:

# create a first order accurate test error function C1dt
C1 = euler_max[0]/euler_dt[0] + 0.05  # compute constant C1
test_error1 = [C1*dt_ for dt_ in euler_dt]  # compute C1dt

# create a fourth order accurate test error function C4dt^4
C4 = rk4_max[0]/rk4_dt[0] + 0.005  # compute constant C4
test_error4 = [C4*dt_**4 for dt_ in rk4_dt]  # compute C4dt^4


# create a second order accurate test error function C2dt^2
C2 = ab2_max[0]/ab2_dt[0] + 0.001  # compute constant C2
test_error2 = [C2*dt_**2 for dt_ in ab2_dt]  # compute C2dt^2


plt.loglog(euler_dt, euler_max, '-r', label='Euler')
plt.loglog(euler_dt, test_error1, linestyle='--', label='first-order slope')
plt.loglog(rk4_dt, rk4_max, '-b', label='RK4')
plt.loglog(rk4_dt, test_error4, linestyle='--', label='fourth-order slope')
plt.loglog(ab2_dt, ab2_max, '-m', label='AB2')
plt.loglog(ab2_dt, test_error2, linestyle='--', label='second-order slope')
plt.xlabel(r'$dt$')
plt.ylabel(r'$E(t)_{max}$')
plt.title(r'dt vs maximum error')
plt.legend()
plt.savefig("max_error.png", dpi=500)
plt.clf()

#print order of convergence:
convergence_rate(ab2_mean)
