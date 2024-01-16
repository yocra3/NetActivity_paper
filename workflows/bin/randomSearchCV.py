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
import functools
import operator

from numpy import array, argmax
from sklearn.model_selection import RandomizedSearchCV
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Dense, Dropout, Activation, Flatten, Input
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier

sys.path.append('./')
import randomconfig
from model import model_function

A = open('input.pb', 'rb')
[x_train, y_train] = pickle.load(A)
A.close()

name = sys.argv[1]
jobs = int(sys.argv[2])

randomconfig.param_distributions['input_shape'] = [ (x_train.shape[1], 1) ]
randomconfig.param_distributions['num_classes'] = [ len(y_train[0]) ]

my_classifier = KerasClassifier(model_function, batch_size = 128)
random = RandomizedSearchCV(my_classifier,
          randomconfig.param_distributions,
          scoring = 'neg_log_loss',
          pre_dispatch = 1,
          n_jobs = 1,
          n_iter = 60,
          cv = 5,
          random_state = 42)
random.fit(x_train, y_train, epochs = 2, workers = jobs)

del random.best_estimator_
with open( name + '_results_iter.tsv', 'w') as csv_file:
    writer = csv.writer(csv_file, delimiter = '\t', escapechar=' ', quoting = csv.QUOTE_NONE)
    for key, value in random.cv_results_.items():
       writer.writerow([key + '\t' + '\t'.join([str(item) for item in value ])])

pickle.dump( random, open( name + "_model.pb", "wb" ) )
