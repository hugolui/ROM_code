import tensorflow as tf
import numpy as np
# import matplotlib.pyplot as plt
from create_placeholders import create_placeholders 
from init_params import init_params
from forward_propagation import forward_propagation
from compute_cost import compute_cost

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
      
    #            # plot the cost
    #            plt.semilogy(np.squeeze(costs))
    #            plt.ylabel('Cost Value')
    #            plt.xlabel('iterations (x1000)')
    #            plt.title("Learning rate =" + str(learning_rate))
    #            #plt.show()
    #            plt.savefig('figs/cost_function.png',dpi=600)
    #            plt.close()
          
          # Save the parameters in a variable
          parameters = sess.run(parameters)

    return parameters, costs
