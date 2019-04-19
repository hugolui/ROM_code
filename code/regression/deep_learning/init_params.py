import tensorflow as tf

def init_params(layers_dims):

    parameters = {} # Initialize parameters
    L = len(layers_dims) # Number of layer (Including the input layer)
    #with tf.variable_scope("init_params", reuse=tf.AUTO_REUSE):
    for l in range(1,L):
        parameters["W" + str(l)] = tf.get_variable("W" + str(l), [layers_dims[l], layers_dims[l-1]], initializer = tf.contrib.layers.xavier_initializer( uniform=False))
        parameters["b" + str(l)] = tf.get_variable("b" + str(l), [layers_dims[l], 1], initializer = tf.zeros_initializer())

    return parameters