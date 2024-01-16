#! /usr/GenNet/env_GenNet/bin/python

#'#################################################################################
#'#################################################################################
#'  Train a neural network defined in Keras
#'#################################################################################
#'#################################################################################

import pickle
import csv
import os
import numpy as np
import sys
import scipy
import functools
import operator
import pandas as pd
import h5py
import tensorflow as tf
import time

from numpy import array, argmax
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten, Input, Reshape
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
from tensorflow.keras.regularizers import l2, l1
from scipy.sparse import coo_matrix


gpus = tf.config.list_physical_devices('GPU')
if gpus:
  try:
    # Currently, memory growth needs to be the same across GPUs
    for gpu in gpus:
      tf.config.experimental.set_memory_growth(gpu, True)
    logical_gpus = tf.config.list_logical_devices('GPU')
    print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
  except RuntimeError as e:
    # Memory growth must be set before GPUs have been initialized
    print(e)

sys.path.append("/usr/GenNet/") 
from GenNet_utils.LocallyDirectedConnected_tf2 import LocallyDirected1D
  
sys.path.append('./')
import params
import model

name = sys.argv[1]
jobs = int(sys.argv[2])


f = h5py.File('assay_reshaped.h5', 'r')
proj_labels = f['label'].attrs['labels']
f.close()

A = open('train.pb', 'rb')
[x_train, y_train] = pickle.load(A)
A.close()

A = open('test.pb', 'rb')
[x_test, y_test] = pickle.load(A)
A.close()

A = open('gene_mask.pb', 'rb')
gene_mask = pickle.load(A)
A.close()

num_classes = len(y_train[0])

# Train model ####
mod = model.model_generator_train(num_classes, gene_mask, params)
opt = Adam(learning_rate = params.alpha)

# Create a callback that saves the model's weights
checkpoint_path = "./checkpoints/"
os.mkdir(checkpoint_path)
class SaveWeights(tf.keras.callbacks.Callback):                                                 
 def on_epoch_end(self, epoch, logs=None):
   w = self.model.get_weights()
   pickle.dump(w, open(checkpoint_path + str(epoch) + '.pb', 'wb'))

mod.compile(loss='categorical_crossentropy',
  optimizer = opt,
  metrics = ['categorical_accuracy'])

mod.summary()

history = mod.fit(x_train, y_train,
  batch_size = params.batch_size,
  epochs = params.epochs,
  verbose = 1, callbacks = [SaveWeights()], validation_data = (x_test, y_test), workers = jobs)

history_dict = history.history
pickle.dump( history_dict, open( name + "_history_model.pb", "wb" ), protocol = 4 )
pickle.dump( proj_labels, open( name + "_labels.pb", "wb" ), protocol = 4 )
