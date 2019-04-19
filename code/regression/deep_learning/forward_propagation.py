import tensorflow as tf

def forward_propagation(X, parameters):
    
    A = X # Inicializing the method A[0] = X
    L = len(parameters)//2 # Number of hidden layers + output layer
    
    for l in range(1,L+1):
        A_prev = A # A[l-1] 
        W = parameters["W" + str(l)] # W parameters
        b = parameters["b" + str(l)] # b parameters
        Z = tf.matmul(W,A_prev) + b # Compute Z
        #A = tf.nn.tanh(Z) # Compute activation function (TANH)
        A = tf.nn.elu(Z) # Compute activation function (ELU)
    
    return Z
