import numpy as np

def func(X, parameters):
    
    A = X # Inicializing the method A[0] = X
    L = len(parameters)//2 # Number of hidden layers + output layer

    for l in range(1,L+1):
        A_prev = A # A[l-1] 
        W = parameters["W" + str(l)] # W parameters
        b = parameters["b" + str(l)] # b parameters
        Z = np.matmul(W,A_prev) + b # Compute Z
        #A = np.tanh(Z)
        A = np.zeros(Z.shape)
        for j in range(Z.shape[0]):
            if Z[j] <=0:
                A[j] = np.exp(Z[j]) - 1.0
            else:
                A[j] = Z[j]

    return Z