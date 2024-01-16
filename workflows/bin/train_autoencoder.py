#! /usr/local/bin/python

#'#################################################################################
#'#################################################################################
#'  Train autoencoder from CNN network
#'#################################################################################
#'#################################################################################

import pickle
import csv
import numpy as np
import sys
import scipy
import functools
import operator
import pandas as pd
import h5py
import tensorflow as tf
import os
import re

from numpy import array, argmax
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Dense, Dropout, Activation, Flatten, Input
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint

name = sys.argv[1]
jobs = int(sys.argv[2])

sys.path.append('./')
import model
import params

A = open('train.pb', 'rb')
[x_train, y_train] = pickle.load(A)
A.close()

A = open('test.pb', 'rb')
[x_test, y_test] = pickle.load(A)
A.close()

## Create input for autoencoder
# Train model ####
model = model.model_generator_train(x_train, y_train, params)
opt = Adam(learning_rate = params.alpha)
model.compile(loss='mean_squared_error',
  optimizer = opt,
  metrics = ['mse'])

model.summary()
history = model.fit(x_train, x_train,
 batch_size = params.batch_size,
 epochs = params.epochs,
 verbose = 1, validation_data = (x_test, x_test), workers = jobs)

history_dict = history.history
pickle.dump( history_dict, open( name + "_history_model.pb", "wb" ), protocol = 4 )
model.save( name + '_autoencoder')
