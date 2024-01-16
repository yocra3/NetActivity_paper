#! /usr/local/bin/python

#'#################################################################################
#'#################################################################################
#'  Transfer learning using features from TCGA model
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

from numpy import array, argmax
from sklearn.model_selection import RandomizedSearchCV
from keras.models import Sequential, Model
from keras.layers import Conv1D, MaxPooling1D, Dense, Dropout, Activation, Flatten, Input
from keras.wrappers.scikit_learn import KerasClassifier
from keras.callbacks import EarlyStopping, ModelCheckpoint

name = sys.argv[1]
last_layer = sys.argv[2]

 
model = tf.keras.models.load_model('model/' + os.listdir('model/')[0]) 

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

for i in range(len(model.layers)):
  model.layers[i].trainable = False
  if model.layers[i].name == last_layer:
    break


model.pop()
model.add(Dense(num_classes, activation="softmax", name="dense_output"))
model.summary()
model.compile(loss='categorical_crossentropy',
  optimizer = 'adam',
  metrics = ['categorical_accuracy'])
callbacks = [EarlyStopping(monitor = 'val_loss', patience = 20, verbose = 1)]

history = model.fit(x_train, y_train,
  batch_size = 64,
  epochs = 1000,
  verbose = 1, callbacks = callbacks, validation_data = (x_test, y_test))


history_dict = history.history
pickle.dump( history_dict, open( name + "_transfer_model_history.pb", "wb" ), protocol = 4 )
pickle.dump( proj_labels, open( name + "_labels.pb", "wb" ), protocol = 4 )
model.save( name )
