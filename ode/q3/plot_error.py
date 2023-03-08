#plot error vs dt for RK4 scheme and show fourth-order convergence:

import os
import numpy as np
import matplotlib.pyplot as plt
from math import *

#----------------------function definitions--------------------------------------#

# analytical solution:
f = lambda t : exp(-t/8.0)*( cos( (3.0*sqrt(7.0)*t)/8.0 ) + (1.0/ (3.0*sqrt(7.0)))*sin( (3.0*sqrt(7.0)*t)/8.0 ) )

#get dt, mean error from file:
def postprocess(filename, f):
    
    t, x = np.loadtxt(filename, unpack=True)
    F = np.vectorize(f)
    x_exact = F(t)
    error = np.absolute(x - x_exact)
    mean_error = np.mean(error)
   
    #get scheme and dt:
    filename = filename.replace('.dat', '')
    dt = float(filename)

    return dt, mean_error

#print convergence rate:
def convergence_rate(error):
  error.sort(reverse=True)  
  for i in range(1,np.size(error)):
    covergence_rate = (log(error[i-1]) - log(error[i]))/log(2)
    print("order of convergence ~", covergence_rate)



#-----------------------------main program------------------------------------#


#create list to store dt and error:
dt = []
error = []
test_error = []

# loop over .dat files in current directory and get dt and average error
for file in os.listdir():
    if file.endswith(".dat"):
        
        #get file data:
        dt_, error_ = postprocess(file, f)
        dt.append(dt_)
        error.append(error_)

#create a test error function Cdt^4
C = error[0]/dt[0]**4 + 0.01                    #compute constant C
test_error = [C*dt_**4 for dt_ in dt]           #compute Cdt4

# dt vs mean error plot
plt.loglog(dt, error, '-r', marker = 'o', label = 'rk4')
plt.loglog(dt, test_error, linestyle='--', label = 'fourth-order slope')
plt.xlabel(r'$dt$')
plt.ylabel(r'$|E(t)|_{avg}$')
plt.legend()
plt.title(r'$dt$ vs average error for $c = 5$')
plt.savefig("average_error.png", dpi=500)
plt.clf()

#print order of convergence:
convergence_rate(error)

