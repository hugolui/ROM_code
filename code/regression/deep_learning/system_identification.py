import numpy as np
from model import model
from rk45 import RK45
import os
import scipy
import scipy.io as sio
from func import func

def system_identification(X_train,Y_train,layers_dims,lambd,learning_rate,num_iter,h,num_model,nt_s):

    ### Time ### 
    h_FOM = h # Time step - FOM
    h_ROM = h_FOM # Time step - ROM 
    nt = X_train.shape[0] #  Number of snapshots - FOM

    ### Training data  ###
    X_train = X_train.T # X - Train
    Y_train = Y_train.T # Y - Train
    n_var = X_train.shape[0] # Number of variables
    
    ### Train algorithm using Deep Neural Networks ###
    parameters, costs = model(X_train,Y_train,layers_dims,learning_rate,num_iter,lambd,print_cost = True)

    ### System Identification ###
    X_s = np.zeros((n_var,nt), dtype = "float32") # Solution
    X_s = RK45(X_train[:,0],parameters,h_ROM,n_var,nt_s) # Runge Kutta 45
    
    ### Save cost values in .dat ###
    scipy.savetxt('./results/cost_function/cost_function_' +  str(num_model) + '.dat', costs, fmt='%25.18e' )
    ### Save DNN architecture ###
    scipy.savetxt('./results/hyperparameters/architecture_' + str(num_model) + '.dat', layers_dims, fmt='%25.18e')
    ### Save temporal modes in .dat ###
    scipy.savetxt('./results/temporal_modes/temporal_modes_' + str(num_model) + '.dat', X_s.T, fmt='%25.18e' )
    ### Save parameters in .mat ###
    sio.savemat('./results/parameters/params_' + str(num_model) + '.mat', parameters)
    ### Save the regularization parameter ###
    scipy.savetxt('./results/hyperparameters/lambda_' + str(num_model) + '.dat', [lambd] ,fmt='%25.18e')

    return 
