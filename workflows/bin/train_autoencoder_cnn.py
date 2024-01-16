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

model = tf.keras.models.load_model('model/' + os.listdir('model/')[0]) 

A = open('train.pb', 'rb')
[x_train, y_train] = pickle.load(A)
A.close()

A = open('test.pb', 'rb')
[x_test, y_test] = pickle.load(A)
A.close()

## Create input for autoencoder
cnn_part = Model(inputs=model.input, outputs=model.layers[3].output)
x_train_auto = cnn_part.predict(x_train)
x_test_auto = cnn_part.predict(x_test)

model2 = Sequential()

## Add layers from original network
for i in range(len(model.layers) - 1):
  if re.search("dense", model.layers[i].name):
    if len(model2.layers) == 0:
      model2.add(Dense(model.layers[i].units, activation='relu', input_dim = x_train_auto.shape[1]))
    else:
      model2.add(Dense(model.layers[i].units, activation='relu'))
    model2.layers[-1].set_weights(model.layers[i].get_weights())
    
## Add layers for autoencoder
for i in range(len(model2.layers) - 2, -1, -1):
  model2.add(Dense(model2.layers[i].units, activation='relu'))

model2.add(Dense(x_train_auto.shape[1], activation='relu'))
opt = Adam(learning_rate = 0.001)
model2.compile(loss = 'mse',
  optimizer = opt,
  metrics = ['MeanSquaredError'])


model2.summary()

callbacks = [EarlyStopping(monitor = 'val_loss', patience = 5, verbose = 1,  min_delta = 1e-5)]

history_auto = model2.fit(x_train_auto, x_train_auto,
  batch_size = 128,
  epochs = 1000,
  verbose = 1, callbacks = callbacks, validation_data = (x_test_auto, x_test_auto), workers = jobs)

j = 0
for i in range(len(model.layers) - 1):
  if re.search("dense", model.layers[i].name):
    model.layers[i].set_weights(model2.layers[j].get_weights())
    j += 1

model.compile(loss='categorical_crossentropy',
  optimizer = opt,
  metrics = ['categorical_accuracy'])
callbacks = [EarlyStopping(monitor = 'val_loss', patience = 5, verbose = 1)]

history_ori = model.fit(x_train, y_train,
  batch_size = 128,
  epochs = 1,
  verbose = 1, callbacks = callbacks, validation_data = (x_test, y_test), workers = jobs)


history_dict_auto = history_auto.history
history_dict_ori = history_ori.history

pickle.dump( history_dict_ori, open( name + "_with_autoencode_history_model.pb", "wb" ), protocol = 4 )
pickle.dump( history_dict_auto, open( name + "_autoencoder_history_model.pb", "wb" ), protocol = 4 )
model.save( name + '_with_autoencode')
model2.save( name + '_autoencoder')

 
