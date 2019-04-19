import tensorflow as tf

def compute_cost(ZL,Y,parameters,n_hidden_layers,lambd,m):

    reg = 0
    for l in range(1,n_hidden_layers):
        reg += tf.nn.l2_loss(parameters["W" + str(l)])

    cost = tf.reduce_mean(tf.square(ZL - Y)) + (lambd*reg)/m

    return cost