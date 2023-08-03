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

from numpy import array, argmax
from sklearn.model_selection import RandomizedSearchCV
from keras.models import Sequential, Model
from keras.optimizers import Adam
from keras.layers import Conv1D, MaxPooling1D, Dense, Dropout, Activation, Flatten, Input
from keras.wrappers.scikit_learn import KerasClassifier
from keras.callbacks import EarlyStopping, ModelCheckpoint

f = h5py.File('assay_reshaped.h5', 'r')
proj_labels = f['label'].attrs['labels']
f.close()

A = open('train.pb', 'rb')
[x_train, y_train] = pickle.load(A)
A.close()

A = open('test.pb', 'rb')
[x_test, y_test] = pickle.load(A)
A.close()

A = open('model.pb', 'rb')
[model, labels] = pickle.load(A)
A.close()

input_shape = (x_train.shape[1], 1)
num_classes = len(y_train[0])

for i in range(6):
    model.layers[i].trainable = False

ll = model.layers[-3].output
ll = Dense(256, activation="relu", name="dense_layer")(ll)
ll = Dense(num_classes, activation="softmax", name="dense_output")(ll)

new_model = Model(inputs=model.input, outputs=ll)
new_model.summary()
opt = Adam(learning_rate = 0.0001)
new_model.compile(loss='categorical_crossentropy',
  optimizer = opt,
  metrics = ['categorical_accuracy'])
callbacks = [EarlyStopping(monitor = 'val_loss', patience = 10, verbose = 1)]

history = new_model.fit(x_train, y_train,
  batch_size = 64,
  epochs = 1000,
  verbose = 1, callbacks = callbacks, validation_data = (x_test, y_test))


history_dict = history.history
pickle.dump( history_dict, open( "history_model.pb", "wb" ), protocol = 4 )
pickle.dump( [ new_model, proj_labels], open( "new_model.pb", "wb" ), protocol = 4 )
