import scipy.io as sio
import numpy as np
from rk45 import RK45
import matplotlib.pyplot as plt
import os

### Inputs ###
h = 0.03356 # Time step
nt = 2244 # Number of snapshots for reconstruction

parameters = sio.loadmat('./data/parameters.mat') # Load parameters
X0 = np.genfromtxt('./data/V_matrix.dat', dtype = "float32") # Initial condition a(t=0)
n_var = X0.shape[0] # Number of variables
t = np.linspace(0,h*(nt-1),nt) # Time array (ROM)

### Reconstruction ###
X_s = RK45(X0,parameters,h,n_var,nt)
X_s = X_s.T

os.system('rm -r ./figs/') 
os.system('mkdir -p ./figs/')
   
### Plot figures ### 
for i in range(n_var):    
    fig, ax = plt.subplots()
    ax.plot(t,X_s[:,i],'c--',label='Prediction',linewidth=0.8)
    plt.xlabel('Time (s)')
    plt.ylabel('x' + str(i+1))
    ax.legend()
    ax.set_ylim([1.1*np.min(X_s[:,i]), 1.1*np.max(X_s[:,i])])
    fig.tight_layout()
    fig.subplots_adjust(right=0.9) 
    plt.savefig('./figs/a' + str(i+1) +'.pdf',dpi=600)
    plt.close()

### Write temporal modes ###
file = open('./data/temporal_modes.dat',"w")
for x in X_s[:nt]:
    for y in x[:n_var]:
        file.write('%.15E   '  % y )
    
    file.write(' \n') 

file.close()
#
### Write the number of snapshots (ROM)
file = open('./data/info.dat',"w")
file.write(str(h)  + ' \n')
file.write(str(nt)  + ' \n')



