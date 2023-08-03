#! /usr/local/bin/python

#'#################################################################################
#'#################################################################################
#' Run random search on tcga data
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

from numpy import array, argmax
from sklearn.model_selection import RandomizedSearchCV
from keras.models import Sequential, Model
from keras.optimizers import Adam
from keras.layers import Conv1D, MaxPooling1D, Dense, Dropout, Activation, Flatten, Input
from keras.wrappers.scikit_learn import KerasClassifier

sys.path.append('./')
import randomconfig

A = open('input.pb', 'rb')
[x_train, y_train] = pickle.load(A)
A.close()

name = sys.argv[2]



## Model
def make_model(dense_layer_sizes, filters, kernel_size, stride, alpha, pool, dense_layer_sizes2):
    input_shape = (x_train.shape[1], 1)
    num_classes = len(y_train[0])

    model = Sequential()
    ## *********** First layer Conv
    model.add(Conv1D(filters, kernel_size = kernel_size, strides = min(stride, kernel_size),
      input_shape = input_shape))
    model.add(Activation('relu'))
    model.add(MaxPooling1D(pool))
    model.output_shape

    ## ********* Classification layer
    model.add(Flatten())
    model.add(Dense(dense_layer_sizes, activation='relu'))
    model.add(Dense(dense_layer_sizes2, activation='relu'))
    model.add(Dense(num_classes, activation='softmax'))
    model.output_shape
    opt = Adam(learning_rate = alpha)
    model.compile(loss = 'categorical_crossentropy',
                  optimizer = opt,
                  metrics = ['accuracy'])
    model.summary()
    return model

my_classifier = KerasClassifier(make_model, batch_size = 128)
random = RandomizedSearchCV(my_classifier,
          randomconfig.param_distributions,
          scoring = 'neg_log_loss',
          pre_dispatch = 1,
          n_jobs = 1,
          n_iter = 60,
          cv = 5,
          random_state = 42)
random.fit(x_train, y_train)

with open('results_50iter_' + name + '.tsv', 'w') as csv_file:
    writer = csv.writer(csv_file, delimiter = '\t', escapechar=' ', quoting = csv.QUOTE_NONE)
    for key, value in random.cv_results_.items():
       writer.writerow([key + '\t' + '\t'.join([str(item) for item in value ])])

pickle.dump( random, open( "model.pb", "wb" ) )
