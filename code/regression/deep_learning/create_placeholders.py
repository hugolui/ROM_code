import tensorflow as tf

def create_placeholders(n_x, n_y):

    # n_x - Number of features
    # n_y - number of classes

    X = tf.placeholder(tf.float32, shape = [n_x, None], name = "X")
    Y = tf.placeholder(tf.float32, shape = [n_y, None], name = "Y")

    return X, Y