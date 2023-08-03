##################################################################################
# DNN gene expression autoencoder pathways network 6
##################################################################################
##################################################################################

from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Input, Dropout
from tensorflow.keras.optimizers import Adam
import tensorflow_model_optimization as tfmot
import numpy as np

def model_generator_train(x_train, y_train, gene_mask, params):

 model = Sequential()
 model.add(tfmot.sparsity.keras.prune_low_magnitude(Dense(gene_mask.shape[1]*params.prepath_neurons,
   input_dim = x_train.shape[1], activation = params.activation),
   tfmot.sparsity.keras.ConstantSparsity(0.5, 1000000, end_step =  1000000, frequency = 100)))
 model.add(tfmot.sparsity.keras.prune_low_magnitude(Dense(gene_mask.shape[1],
   activation = params.activation),
   tfmot.sparsity.keras.ConstantSparsity(0.5, 1000000, end_step =  1000000, frequency = 100)))
 if params.dense1_dropout:
    model.add(Dropout(params.dense1_dropout))
 model.add(Dense(x_train.shape[1]))

## Create gene mask for pre-pathways layer
 gene_mask2 = np.zeros((gene_mask.shape[0], gene_mask.shape[1]*params.prepath_neurons), float)
 for col in range(gene_mask.shape[1]):
      for i in range(params.prepath_neurons):
          gene_mask2[:, col*params.prepath_neurons + i] = gene_mask[:, col]

  ## Create mask for pre-pathways layer
 gene_mask3  = np.zeros((gene_mask.shape[1]*params.prepath_neurons, gene_mask.shape[1]), float)
 for col in range(gene_mask.shape[1]):
     gene_mask3[col*params.prepath_neurons:col*params.prepath_neurons + params.prepath_neurons, col] = 1

 ## Add gene mask
 w = model.get_weights()
 w[2] = gene_mask2
 w[7] = gene_mask3
 model.set_weights(w)
 opt = Adam(learning_rate = params.alpha)
 model.compile(loss='mse',
   optimizer = opt,
   metrics = ['mse'])
 return model
