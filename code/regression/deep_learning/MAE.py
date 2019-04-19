import numpy as np
import scipy.io as sio
import heapq
import os

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

n_models = int(inputs[0]) # Number of models
n_models_eval = int(inputs[10]) # Number of best models to be evaluated 

X_train = np.genfromtxt('./data/V_matrix.dat', skip_header = 4)
m = X_train.shape[0] # Training size
n_var = X_train.shape[1] # Number of variables

MAE_score = np.zeros(n_models) # Mean Absolute Error (MAE)
### Load candidate models results ###
for i in range(n_models):
  X_dummy = np.genfromtxt('./results/temporal_modes/temporal_modes_' + str(i) +'.dat') 
  MAE = 0
  for j in range(n_var):
    X_s = X_dummy[:m,j]
    ### Mean Absolute Error (MAE)
    MAE += (1/m)*np.sum(abs(X_train[:,j] - X_s))
  MAE_score[i] = MAE
  
### Best models ###
best_models = heapq.nsmallest(int(n_models_eval), range(len(MAE_score)), MAE_score.take) 
file = open('./results/best_models.dat',"w")
for x in best_models:
  file.write('%d'  % x + ' \n') 
file.close()

