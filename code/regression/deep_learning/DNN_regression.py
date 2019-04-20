import numpy as np
import time
import random
import scipy
from system_identification import system_identification
import os

#start_time = time.time() # Start timer

print ('Computing DNN Regression')

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
n_layers_min = int(inputs[1]) # Minimum number of layers
n_layers_max = int(inputs[2]) # Maximum number of layers
n_hidden_min = int(inputs[3]) # Minimum number of hidden units per layers
n_hidden_max = int(inputs[4]) # Maximum number of hidden units per layers
lambd_min = int(inputs[5]) # Minimum magnitude order of the regularization parameter
lambd_max = int(inputs[6]) # Maximum magnitude order of the regularization parameter
learning_rate = float(inputs[7]) # Learning rate
num_iter = int(inputs[8]) # Number of iterations
nt_s = int(inputs[9]) # Number of snapshots for reconstruction
opt = int(inputs[11]) # Hyperparameter optimization: Random search (1); Bayesian Optimization (2)

### Load data ###
with open('./data/V_matrix.dat') as V:
    V_matrix = V.read()
norm = V_matrix[0] # Get POD Correlation matrix norm 
V_matrix = V_matrix.splitlines() # Split it into lines
h = float(V_matrix[3]) # Time step
X_train = np.genfromtxt('./data/V_matrix.dat', dtype = "float32", skip_header = 4) # Training data X
Y_train = np.genfromtxt('./data/da_' + norm + '.dat', dtype = "float32") # # Training data Y
n_var = X_train.shape[1] # Number of variables

# Delete and create folders
os.system('rm -r ./results/figs/')
os.system('rm -r ./results/cost_function/')
os.system('rm -r ./results/hyperparameters/')
os.system('rm -r ./results/parameters/')
os.system('rm -r ./results/temporal_modes/')
os.system('mkdir -p ./results/figs')
os.system('mkdir -p ./results/cost_function')
os.system('mkdir -p ./results/hyperparameters')
os.system('mkdir -p ./results/parameters')
os.system('mkdir -p ./results/temporal_modes')

### RANDOM SEARCH ###
if (opt == 1):

  for i in range(n_models): 
    lambd = 10**(-random.uniform(lambd_min,lambd_max)) # Regularization parameter
    layers_dims = np.zeros(int(random.uniform(n_layers_min,n_layers_max))) # DNN architecture
    layers_dims[0] = n_var # Input Layer size
    layers_dims[layers_dims.shape[0]-1] = n_var # Output layer size
    for j in range(1,layers_dims.shape[0]-1):
      layers_dims[j] = int(random.uniform(n_hidden_min,n_hidden_max)) # Hidden layers
      
    system_identification(X_train,Y_train,layers_dims,lambd,learning_rate,num_iter,h,i,nt_s) # Regression step via DNN

  print ('DNN Regression complete!')

### Bayesian Optimization ###
if (opt == 2):

  import numpy as np
  import skopt
  from skopt import gp_minimize, forest_minimize
  from skopt.space import Real, Categorical, Integer
  from skopt.plots import plot_convergence
  from skopt.plots import plot_objective, plot_evaluations
  from skopt.utils import use_named_args
  import scipy
  import tensorflow as tf
  import matplotlib.pyplot as plt
  from create_placeholders import create_placeholders 
  from init_params import init_params
  from forward_propagation import forward_propagation
  from compute_cost import compute_cost
  from rk45 import RK45
  import scipy.io as sio
  import os
  import time

  # Lambda space #
  dim_lambd = Real(low=10**(-lambd_max), high=10**(-lambd_min), prior='log-uniform',name='lambd')

  # Number of layers space #
  dim_num_layers = Integer(low=n_layers_min, high=n_layers_max, name='num_layers')

  # Number of units space #
  dim_num_nodes = Integer(low=n_hidden_min, high=n_hidden_max, name='num_nodes')

  dimensions = [dim_num_layers,dim_num_nodes,dim_lambd]

  default_parameters = [n_layers_min, n_hidden_min, 10**(-random.uniform(lambd_min,lambd_max))] # Initial guesss

  def model(X_train,Y_train,layers_dims,learning_rate,num_iter,lambd,print_cost):
      
      with tf.device('/device:GPU:0'):
        tf.reset_default_graph() # to be able to rerun the model without overwriting tf variables
        (n_x, m) = X_train.shape # Number of features and number of training examples
        n_y = Y_train.shape[0] # Number of classes
        n_hidden_layers = len(layers_dims) # Number of hidden layers
        costs = [] # Keep track of the cost
        
        ### Create Placheholders ###
        X, Y = create_placeholders(n_x,n_y)
        
        ### Initialize Parameters ###
        parameters = init_params(layers_dims)
        
        ### Foward propagation - Build the forward propagation in the tensorflow graph ###
        ZL = forward_propagation(X,parameters)
        
        ### Cost - Add cost function to tensorflow graph ###
        cost_function = compute_cost(ZL,Y,parameters,n_hidden_layers,lambd,m)
        
        ### Backpropagation - Define the tensorflow optimizer ###
        optimizer = tf.train.AdamOptimizer(learning_rate = learning_rate).minimize(cost_function)
        #optimizer = tf.train.GradientDescentOptimizer(learning_rate = learning_rate).minimize(cost_function)
        
        ### Initializer all the variables ###
        init = tf.global_variables_initializer()
        
        ### Start the session to compute the tensorflow graph ###
        with tf.Session() as sess:
            
            # Run the initialization
            sess.run(init)
            
            # Do Training lopp #
            for i in range(num_iter):

                # Run the session to execute the optimizer and the cost
                _, cost_value = sess.run([optimizer,cost_function], feed_dict = {X: X_train, Y: Y_train})
        
                # Print the cost every 1000 iterations
                #if print_cost == True and i % 1000 == 0:
                #    print ("Cost after iteration %i: %f" % (i, cost_value))
                if print_cost == True and i % 1000 == 0:
                    costs.append(cost_value)
            # Save the parameters in a variable
            parameters = sess.run(parameters)

      return parameters, costs
  
  num_model = 0  
  error_best = 1

  @use_named_args(dimensions=dimensions)
  def dl_sindy_bayesian(num_layers,num_nodes,lambd):

      global num_model 
      
      ### Load data ###
      with open('./data/V_matrix.dat') as V:
          V_matrix = V.read()
      V_matrix = V_matrix.splitlines() # Split it into lines
      h = float(V_matrix[3]) # Time step
      X = np.genfromtxt('./data/V_matrix.dat', dtype = "float32", skip_header = 4) # X
      Y = np.genfromtxt('./data/da_P.dat', dtype = "float32") # Y
      X_train = X # Training data - X
      Y_train = Y # Training data - Y
      n_var = X.shape[1] # Number of variables
          
      ### Time ### 
      h_FOM = h # Time step - FOM
      h_ROM = h_FOM # Time step - ROM 
      nt = X_train.shape[0] #  Number of snapshots - FOM
      t = np.linspace(0,h_FOM*(nt-1),nt) # Time array FOM 
      t_s =  np.linspace(0,h_ROM*(nt_s-1),nt_s) # Time array ROM

      ### Training data  ###
      X_train = X_train.T # X - Train
      Y_train = Y_train.T # Y - Train

      # Network architecture
      layers_dims = np.zeros(num_layers)
      layers_dims[1:-1] = num_nodes
      layers_dims[0] = n_var
      layers_dims[-1] = n_var    
      
      ### Train algorithm using Deep Neural Networks ###
      parameters, costs = model(X_train,Y_train,layers_dims,learning_rate,num_iter,lambd,print_cost = True)

      ### System Identification ###
      X_s = np.zeros((n_var,nt), dtype = "float32") # Solution
      X_s = RK45(X_train[:,0],parameters,h_ROM,n_var,nt_s) # Runge Kutta 45
      
      ## Plot Graphs ###
      os.system('mkdir -p ./results/figs/model_' + str(num_model) + '/')
      
      # # Plot figures #  
      # for i in range(n_var):    
      #     fig, ax = plt.subplots()
      #     ax.plot(t.T,np.expand_dims(X_train[i,:], axis = 0).T,'r--',label='Original',linewidth=1.0)
      #     ax.plot(t_s.T,np.expand_dims(X_s[i,:], axis = 0).T,'k',label='Prediction',linewidth=0.5)
      #     plt.xlabel('Time (s)')
      #     plt.ylabel('a' + str(i+1))
      #     plt.title('Modo temporal ' + str(i+1))
      #     ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
      #     fig.tight_layout()
      #     fig.subplots_adjust(right=0.8) 
      #     #plt.show()
      #     plt.savefig('./results/figs/model_' + str(num_model) + '/a' + str(i+1) +'.png',dpi=600)
      #     plt.close()
      
      # # plot the cost
      # plt.figure()
      # plt.semilogy(np.squeeze(costs))
      # plt.ylabel('Cost Value')
      # plt.xlabel('iterations (x1000)')
      # plt.title("Learning rate =" + str(learning_rate))
      # plt.savefig('./results/cost_function/cost_function_' + str(num_model) + '.png',dpi=600)
      # plt.close()

      ## Save cost values in .dat ###s
      scipy.savetxt('./results/cost_function/cost_function_' +  str(num_model) + '.dat', costs, fmt='%25.18e' )
      ### Save DNN architecture ###
      scipy.savetxt('./results/hyperparameters/architecture_' + str(num_model) + '.dat', layers_dims, fmt='%25.18e')
      ### Save temporal modes in .dat ###
      scipy.savetxt('./results/temporal_modes/temporal_modes_' + str(num_model) + '.dat', X_s.T, fmt='%25.18e' )
      ## Save parameters in .mat ###
      sio.savemat('./results/parameters/params_' + str(num_model) + '.mat', parameters)
      ### Save the regularization parameter ###
      scipy.savetxt('./results/hyperparameters/lambda_' + str(num_model) + '.dat', [lambd] ,fmt='%25.18e')

      erro_train = 0 # training error
      for i in range(n_var):
        erro_train += ((1/nt)*np.sum(np.square(X_s[i,:nt] - X_train[i,:])))

      error_train = erro_train/n_var
      
      global error_best
      # Save the best model in .dat
      if error_train < error_best:
        # Save the best model in .dat
        file = open('./results/best_model.dat',"w+")
        file.write('%.0d'  % num_model + ' \n') 
        file.close()
        error_best = error_train
      
      num_model = num_model + 1

      return error_train

  ### Bayesian Optimization ###
  search_result = gp_minimize(func=dl_sindy_bayesian,dimensions=dimensions,acq_func='EI', n_calls=n_models,x0=default_parameters)

  model = np.zeros(n_models)
  for i in range(n_models):
    model[i] = i+1
  models_sort = sorted(zip(search_result.func_vals, search_result.x_iters,model))

  with open('./results/best_models.dat','w') as file:
      for i in range(n_models):
          file.write('%i \n' % models_sort[i][2])

  print ('DNN Regression complete!')

### LATIN HYPERCUBE SAMPLING ###
if (opt == 3):

  from smt.sampling_methods import LHS

  lambd = np.array([10**-(lambd_min),10**-(lambd_max)])
  num_layers = np.array([n_layers_min,n_layers_max])
  num_hidden_units = np.array([n_hidden_min,n_hidden_max])
     
  xlimits = np.array([lambd,num_layers,num_hidden_units])
  sampling = LHS(xlimits=xlimits)
    
  XXX = sampling(n_models)
    
  for i in range(n_models):
    XXX[i,1] = int(XXX[i,1])
    XXX[i,2] = int(XXX[i,2])
    layers_dims = np.zeros([int(XXX[0,1])])
    layers_dims[0] = n_var # Input Layer size
    layers_dims[layers_dims.shape[0]-1] = n_var # Output layer size
    for j in range(1,layers_dims.shape[0]-1):
      layers_dims[j] = int(XXX[i,2]) # Hidden layers
    
    system_identification(X_train,Y_train,layers_dims,lambd,learning_rate,num_iter,h,i,nt_s) # Regression step via DNN

  print ('DNN Regression complete!')
