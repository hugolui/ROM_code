import scipy.io as sio
import numpy as np
import os
import shutil
from rk45 import RK45
import matplotlib.pyplot as plt

#### folder path ###
path = os.path.dirname(os.getcwd())

### Delete and create folder ###
os.system('rm -r ./data')
os.system('mkdir -p ./data')
os.system('rm -r ./figs')
os.system('mkdir -p ./figs')

### Copy best_model.dat ###
f = path + '/regression/deep_learning/results/best_models.dat'
shutil.copy(f, './data/')
best_models = np.genfromtxt('./data/best_models.dat') # Best models

### Copy inputs.inp ###
f = path + '/regression/deep_learning/data/inputs.inp'
shutil.copy(f, './data/')

### Copy SVD files ###
files = []
files.append(path + '/pod/SVD_files/initial_time.dat')
files.append(path + '/pod/SVD_files/sigma.dat')
files.append(path + '/pod/SVD_files/V_matrix.dat')
for f in files:
  shutil.copy(f, './data/')

### POD modes ###
f = path + '/pod/POD_modes.cgns'
shutil.copy(f, './data/')

### Inputs ###
# Load inputs file
with open('./data/inputs.inp') as f:
    inputs = f.read()

# Split it into lines
inputs = inputs.splitlines()

# Find where POD INPUTS key word is
for i,lines in enumerate(inputs):
  if lines=='*DNN INPUTS*': break

# Remove all unwanted content
inputs = [line.replace(' ','').split('!')[0] for line in inputs[i+1:] if line!= '']

nt_s = int(inputs[9]) # Number of snapshots (ROM)

for i in range(len(best_models)):
    model = best_models[i]
    parameters = sio.loadmat(path + '/regression/deep_learning/results/parameters/params_' + str(int(model)) +'.mat') # Load parameters
    architecture = np.genfromtxt(path + '/regression/deep_learning/results/hyperparameters/architecture_' + str(int(model)) +'.dat') # Load DNN architecture
    with open(path + '/pod/SVD_files/V_matrix.dat') as V:
        V_matrix = V.read()
    V_matrix = V_matrix.splitlines() # Split it into lines
    h = float(V_matrix[3]) # Time step
    X = np.genfromtxt(path + '/pod/SVD_files/V_matrix.dat', dtype = "float32", skip_header = 4) # Load data
    X0 = X[0,:] # Initial solution
    n_var = X.shape[1] # Number of variables

    nt = X.shape[0] # Number of snapshots (FOM)
    t = np.linspace(0,h*(nt-1),nt) # Time array (FOM)

    h_s = h # Time step (ROM)
    t_s = np.linspace(0,h_s*(nt_s-1),nt_s) # Time array (ROM)

    ### Reconstruction (ROM)
    X_s = RK45(X0,parameters,h_s,n_var,nt_s)
    X_s = X_s.T

    os.system('mkdir -p ./figs/model_' + str(int(model)) + '/')

    ### Plot figures ### 
    TW_line_x = np.full((100, 1), t[X.shape[0]-1]) # End of the training window 
    for i in range(n_var):    
        fig, ax = plt.subplots()
        ax.plot(t,X[:,i],'k',label='Original',linewidth=0.8)
        ax.plot(t_s,X_s[:,i],'c--',label='Prediction',linewidth=0.8)
        ax.plot(TW_line_x, np.linspace(1.1*np.min(X[:,i]),1.1*np.max(X[:,i]),100),'b--', label = 'End of training window',linewidth=0.8)
        plt.xlabel('Time (s)')
        plt.ylabel('x' + str(i+1))
        ax.legend()
        ax.set_ylim([1.1*np.min(X[:,i]), 1.1*np.max(X[:,i])])
        fig.tight_layout()
        fig.subplots_adjust(right=0.9) 
        plt.savefig('./figs/model_' + str(int(model)) + '/x' + str(i+1) +'.pdf',dpi=600)
        plt.close()

    ### Write temporal modes ###
    file = open('./data/temporal_modes_'+ str(int(model)) + '.dat',"w")
    for x in X_s[:nt_s]:
        for y in x[:n_var]:
            file.write('%.15E   '  % y )
        
        file.write(' \n') 
    
    file.close()

### Write the number of snapshots (ROM)
file = open('./data/nsnap_reconst.dat',"w")
file.write('%d' % nt_s)
 
### Write inputs for reconstruction.f90 ###
# Load inputs file
with open('./data/inputs.inp') as f:
    inputs_reconst = f.read()

# Split it into lines
inputs_reconst = inputs_reconst.splitlines()
inputs_reconst1 = [line.replace(' ','').split('!')[0] for line in inputs_reconst if line!= ''] # Remove all unwanted content
inputs_reconst1 = [line.replace('','').split('!')[0] for line in inputs_reconst if line!= ''] # Remove all unwanted content

# Write inputs_recons #
best_model = '/data/temporal_modes_' + str(int(best_models[0])) + '.dat' 
inputs_reconst = []
inputs_reconst.append(inputs_reconst1[0])
inputs_reconst.append(inputs_reconst1[3])
inputs_reconst.append(inputs_reconst1[1])
inputs_reconst.append(best_model)
inputs_reconst.append(int(best_models[0]))
inputs_reconst.append(inputs_reconst1[8])
inputs_reconst.append(inputs_reconst1[9])
inputs_reconst.append(inputs_reconst1[11])
inputs_reconst.append(inputs_reconst1[12])
inputs_reconst.append(nt_s)
inputs_reconst.append(h_s)
inputs_reconst.append('F')

file = open('./inputs_reconst.inp',"w")
for x in inputs_reconst:
    file.write('%s   '  % x )
    file.write(' \n') 
file.close()

# Run reconst.f90
os.system('sh compile.sh')
os.system('./reconst.out')

for i in range(1,len(best_models)):
  
  # Load inputs file
  with open('inputs_reconst.inp') as f:
      inputs = f.read()
  
  # Split it into lines
  inputs = inputs.splitlines()
  
  # Replace lines
  inputs[3] = '/data/temporal_modes_'+ str(int(best_models[i])) +'.dat'
  inputs[4] =  str(int(best_models[i]))
  
  # Write new inputs 
  file = open('inputs_reconst.inp',"w")
  for x in inputs:
    file.write('%s'  % x + ' \n') 
  file.close()
  
  # Run reconst.f90
  os.system('sh compile.sh')
  os.system('./reconst.out')


### Find best model ###
L1_error = np.genfromtxt('./data/L1_error.dat') # L1_error of each model
min_pos_L1 = np.argmin(L1_error)
best_model_pos = best_models[min_pos_L1] # Best candidate model

# Load inputs file
with open('inputs_reconst.inp') as f:
    inputs = f.read()

# Split it into lines
inputs = inputs.splitlines()

# Replace lines
inputs[3] = '/data/temporal_modes_'+ str(int(best_model_pos)) +'.dat'
inputs[4] =  str(int(best_model_pos))
inputs[11] = 'T'

# Write inputs 
file = open('./inputs_reconst.inp',"w")
for x in inputs:
    file.write('%s   '  % x )
    file.write(' \n') 
file.close()

# Run reconst.f90
os.system('sh compile.sh')
os.system('./reconst.out')


