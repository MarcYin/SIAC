import numpy as np
from numba import jit
from numba import jit
# the forward and backpropogation 
# are from https://medium.com/unit8-machine-learning-publication/computing-the-jacobian-matrix-of-a-neural-network-in-python-4f162e5db180
# but added jit for faster speed in the calculation

@jit(nopython=True)
def affine_forward(x, w, b):
    """
    Forward pass of an affine layer
    :param x: input of dimension (D, )
    :param w: weights matrix of dimension (D, M)
    :param b: biais vector of dimension (M, )
    :return output of dimension (M, ), and cache needed for backprop
    """
    out = np.dot(x, w) + b
    cache = (x, w)
    return out, cache

@jit(nopython=True)
def affine_backward(dout, cache):
    """
    Backward pass for an affine layer.
    :param dout: Upstream Jacobian, of shape (O, M)
    :param cache: Tuple of:
      - x: Input data, of shape (D, )
      - w: Weights, of shape (D, M)
    :return the jacobian matrix containing derivatives of the O neural network outputs with respect to
            this layer's inputs, evaluated at x, of shape (O, D)
    """
    x, w = cache
    dx = np.dot(dout, w.T)
    return dx

@jit(nopython=True)
def relu_forward(x):
    """ Forward ReLU
    """
    out = np.maximum(np.zeros(x.shape).astype(np.float32), x)
    cache = x
    return out, cache

@jit(nopython=True)
def relu_backward(dout, cache):
    """
    Backward pass of ReLU
    :param dout: Upstream Jacobian
    :param cache: the cached input for this layer
    :return: the jacobian matrix containing derivatives of the O neural network outputs with respect to
             this layer's inputs, evaluated at x.
    """
    x = cache
    dx = dout * np.where(x > 0, np.ones(x.shape).astype(np.float32), np.zeros(x.shape).astype(np.float32))
    return dx

@jit(nopython=True)
def sigmod_forward(x):
    """ Forward sigmod
    """
    out = 1 / (1 + np.exp(-1*x))
    return out

def forward_backward(x, Hidden_Layers, Output_Layers, cal_jac=False): 
    layer_to_cache = dict()  # for each layer, we store the cache needed for backward pass 
    [[w1, b1], [w2, b2]] = Hidden_Layers
    a1, cache_a1 = affine_forward(x, w1, b1) 
    r1, cache_r1 = relu_forward(a1) 
    a2, cache_a2 = affine_forward(r1, w2, b2) 
    rets = []
    for output_layer in Output_Layers:
        w3, b3 = output_layer
        r3, cache_r3 = relu_forward(a2)
        out, cache_out = affine_forward(r3, w3, b3)
        if cal_jac:
            dout = affine_backward(np.ones_like(out), cache_out)
            dout = relu_backward(dout, cache_r3) 
            dout = affine_backward(dout, cache_a2) 
            dout = relu_backward(dout, cache_r1) 
            dx = affine_backward(dout, cache_a1)
            ret = [out, dx]
        else:
            ret = out
        
        rets.append(ret)
    return rets

def training(X, targs, nodes = 64, epochs = 2000, batch_size=60):
    try:
        import tensorflow as tf              
        from tensorflow import keras         
        from tensorflow.keras import layers 
    except ImportError:
        raise ImportError('To train a nueral network, you need to install tensorflow and keras.')
    inputs = layers.Input(shape=(X.shape[1],))    
    x = layers.Dense(nodes, activation='relu')(inputs)
    x = layers.Dense(nodes, activation='relu')(x)
    outputs = []
    for i in range(targs.shape[1]):
        outputs.append(layers.Dense(1)(x))
    model = tf.keras.Model(inputs=inputs, outputs=outputs)
    optimizer = tf.keras.optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0, amsgrad=False)
    model.compile(optimizer=optimizer,   
                  loss='mean_squared_error',
                  metrics=['mean_absolute_error'])
    history = model.fit(X, [targs[:,i] for i in range(targs.shape[1])], epochs=epochs, batch_size=batch_size)
    return model, history

def relearn(X, targs, model, epochs = 2000, batch_size=60):
    try:
        import tensorflow as tf              
        from tensorflow import keras         
        from tensorflow.keras import layers 
    except ImportError:
        raise ImportError('To train a nueral network, you need to install tensorflow and keras.')
    history = model.fit(X, [targs[:,i] for i in range(targs.shape[1])], epochs=epochs, batch_size=batch_size)
    return model, history

def get_layers(model):
    l1 = model.get_layer(index=1).get_weights()
    l2 = model.get_layer(index=2).get_weights()
    Hidden_Layers = [l1, l2]
    Output_Layers = []  
    for layer in model.layers[3:]:
        l3 = layer.get_weights()
        Output_Layers.append(l3)
    return Hidden_Layers, Output_Layers

def save_tf_model(model, fname):
    try:
        import tensorflow as tf              
        from tensorflow import keras         
        from tensorflow.keras import layers 
    except ImportError:
        raise ImportError('To use a tensorflow model, you need to install tensorflow and keras.')
    model.save(fname)

def load_tf_Model(fname):
    try:
        import tensorflow as tf              
        from tensorflow import keras         
        from tensorflow.keras import layers 
    except ImportError:
        raise ImportError('To use a tensorflow model, you need to install tensorflow and keras.')
    model  = tf.keras.models.load_model(fname)
    return model

def save_np_model(fname, Hidden_Layers, Output_Layers):
    Hidden_Layers[0][1] = np.atleast_2d(Hidden_Layers[0][1]).astype(np.float32)
    Hidden_Layers[1][1] = np.atleast_2d(Hidden_Layers[1][1]).astype(np.float32)
    Output_Layers[0][1] = np.atleast_2d(Output_Layers[0][1]).astype(np.float32)
    np.savez(fname, Hidden_Layers = Hidden_Layers, Output_Layers = Output_Layers)

def load_np_model(fname):
    f = np.load(fname,  allow_pickle=True)
    Hidden_Layers=f.f.Hidden_Layers.tolist() 
    Output_Layers=f.f.Output_Layers.tolist()
    Hidden_Layers[0][1] = (Hidden_Layers[0][1]).ravel().astype(np.float32)
    Hidden_Layers[1][1] = (Hidden_Layers[1][1]).ravel().astype(np.float32)
    Output_Layers[0][1] = (Output_Layers[0][1]).ravel().astype(np.float32)
    return Hidden_Layers, Output_Layers


class Two_NN(object):
    def __init__(self,
                tf_model      = None,
                tf_model_file = None,
                np_model_file = None,
                Hidden_Layers = None,
                Output_Layers = None,
                ):

        if tf_model_file is not None:
            self.tf_model_file = tf_model_file
            self.tf_model = load_tf_Model(self.tf_model_file)
            self.Hidden_Layers, self.Output_Layers = get_layers(self.tf_model)

        if tf_model      is not None:
            self.tf_model      = tf_model
            self.Hidden_Layers, self.Output_Layers = get_layers(self.tf_model)
       
        if np_model_file is not None:
            self.np_model_file = np_model_file
            self.Hidden_Layers, self.Output_Layers = load_np_model(np_model_file)

        if (Hidden_Layers is not None) & (Output_Layers is not None):
            self.Hidden_Layers = Hidden_Layers 
            self.Output_Layers = Output_Layers 
    
    def train(self, X, targs, nodes = 64, iterations = 2000, batch_size=60, tf_fname = "model.h5", save_tf_model = False):
        if (X is not None) & (targs is not None):
            self.tf_model, self.history = training(X, targs, nodes = nodes, epochs = iterations, batch_size=batch_size)
            self.Hidden_Layers, self.Output_Layers = get_layers(self.tf_model)
            if save_tf_model:
                save_tf_model(model, tf_fname)
        else:
            raise IOError('X and targs need to have values')

    def relearn(self, X, targs, iterations = 2000, batch_size=60):
        if hasattr(self, 'tf_model'):
            self.tf_model, self.history = relearn(X, targs, self.tf_model, epochs = iterations, batch_size=batch_size)
            self.Hidden_Layers, self.Output_Layers = get_layers(self.tf_model)
        else:
            raise NameError('No tf model to relearn.')

    def predict(self, x, cal_jac = False):
        if hasattr(self, 'Hidden_Layers') and hasattr(self, 'Output_Layers'):
            x = x.astype(np.float32)
            rets = forward_backward(x, self.Hidden_Layers, self.Output_Layers, cal_jac=cal_jac)
        else:
            raise NameError('Hidden_Layers and Output_Layers have not yet been defined, and please try to train or load a model first.')
        return rets

    def save_tf_model(self, fname):
        if hasattr(self, 'tf_model'):
            save_tf_model(self.tf_model, fname)
            self.tf_model_file = fname
        else:
            raise NameError('No tf model to save.')
    def save_np_model(self, fname):
        if hasattr(self, 'Hidden_Layers') and hasattr(self, 'Output_Layers'):
            save_np_model(fname, self.Hidden_Layers, self.Output_Layers)
            self.np_model_file = fname
        else:
            raise NameError('Hidden_Layers and Output_Layers have not yet been defined, and please try to train or load a model first.')

if __name__ == '__main__' :
    f = np.load('/home/ucfafyi/DATA/Prosail/prosail_2NN.npz')
    v = np.load('/home/ucfafyi/DATA/Prosail/vals.npz')
    tnn = Two_NN(Hidden_Layers=f.f.Hidden_Layers, Output_Layers=f.f.Output_Layers)
    
    refs = tnn.predict(v.f.vals_x)
    from scipy.stats import linregress 
    import pylab as plt
    fig, axs = plt.subplots(ncols = 3, nrows = 3, figsize = (16,16)) 
    axs = axs.ravel() 
    #refs = model.predict(vals_x) 
    vals = v.f.vals 
    for i in range(9): 
        axs[i].plot(refs[i], vals[:,i], 'o', alpha=0.1) 
        #axs[i].set_title(s2a.iloc[100:2100, b_ind[i]+1].name) 
        lin = linregress(refs[i].ravel(), vals[:,i]) 
        print(lin.slope, lin.intercept, lin.rvalue, lin.stderr)

    plt.show()
