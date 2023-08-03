##################################################################################
# DNN gene expression autoencoder pathways network 3
##################################################################################
##################################################################################

from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Input, Dropout
from tensorflow.keras.optimizers import Adam
import tensorflow_model_optimization as tfmot

def model_generator_train(x_train, y_train, gene_mask, params):

  model = Sequential()
  model.add(tfmot.sparsity.keras.prune_low_magnitude(Dense(gene_mask.shape[1],
    input_dim = x_train.shape[1], activation = params.activation),
    tfmot.sparsity.keras.ConstantSparsity(0.5, 1000000, end_step =  1000000, frequency = 100)))
  if params.dense1_dropout:
        model.add(Dropout(params.dense1_dropout))
  model.add(Dense(x_train.shape[1]))

  ## Add gene mask
  w = model.get_weights()
  w[2] = gene_mask
  model.set_weights(w)


  opt = Adam(learning_rate = params.alpha)

  model.compile(loss='mse',
    optimizer = opt,
    metrics = ['mse'])
  return model
