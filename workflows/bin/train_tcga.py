#! /usr/local/bin/python

#'#################################################################################
#'#################################################################################
#'  Train TCGA network
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

from numpy import array, argmax
from sklearn.model_selection import RandomizedSearchCV
from keras.models import Sequential, Model
from keras.optimizers import Adam
from keras.layers import Conv1D, MaxPooling1D, Dense, Dropout, Activation, Flatten, Input
from keras.wrappers.scikit_learn import KerasClassifier
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.regularizers import l2, l1

sys.path.append('./')
import network_config

f = h5py.File('assay_reshaped.h5', 'r')
proj_labels = f['label'].attrs['labels']
f.close()

A = open('train.pb', 'rb')
[x_train, y_train] = pickle.load(A)
A.close()

A = open('test.pb', 'rb')
[x_test, y_test] = pickle.load(A)
A.close()

input_shape = (x_train.shape[1], 1)
num_classes = len(y_train[0])

# Train model ####
model = Sequential()
## *********** First layer Conv
model.add(Conv1D(network_config.filters, kernel_size = network_config.kernel_size,
    strides = network_config.stride,
    input_shape = input_shape))
model.add(Activation('relu'))
model.add(MaxPooling1D(network_config.pool))
## ********* Classification layer
model.add(Flatten())
model.add(Dense(network_config.dense_layer_sizes, activation='relu', kernel_regularizer=l2(0.001)))
model.add(Dense(num_classes, activation='softmax'))
opt = Adam(learning_rate = network_config.alpha)
model.compile(loss='categorical_crossentropy',
  optimizer = opt,
  metrics = ['categorical_accuracy'])
callbacks = [EarlyStopping(monitor = 'val_loss', patience = 5, verbose = 1)]
history = model.fit(x_train, y_train,
  batch_size = 128,
  epochs = 1000,
  verbose = 1, callbacks = callbacks, validation_data = (x_test, y_test))

history_dict = history.history
pickle.dump( history_dict, open( "history_model.pb", "wb" ), protocol = 4 )
pickle.dump( [ model, proj_labels], open( "model.pb", "wb" ), protocol = 4 )
