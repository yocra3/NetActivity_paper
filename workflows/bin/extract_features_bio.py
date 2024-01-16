#! /usr/GenNet/env_GenNet/bin/python

#'#################################################################################
#'#################################################################################
#'  Extract features from bio DNN network
#'#################################################################################
#'#################################################################################

import pickle
import h5py
import csv
import numpy as np
import sys
import scipy
import functools
import operator
import pandas as pd
import tensorflow as tf
import os

from numpy import array, argmax
from tensorflow.keras.models import Model
from tensorflow.keras.optimizers import Adam

sys.path.append('./')
import model
import params

f = h5py.File('assays.h5', 'r')
meth_matrix = f['assay001']
input_matrix = meth_matrix[...]
f.close()

A = open('gene_mask.pb', 'rb')
gene_mask = pickle.load(A)
A.close()


A = open('checkpoints/0.pb', 'rb')
weights = pickle.load(A)
A.close()

num_classes = weights[len(weights) - 1].shape[0]

# Train model ####
mod = model.model_generator_train(num_classes, gene_mask, params)
opt = Adam(learning_rate = params.alpha)

if os.path.isdir('checkpoints/'):
  checks = os.listdir('checkpoints/')
  epochs = [file.split(".")[0] for file in checks]
  for epoch in set(epochs):
    print("Loading checkpoints")
    A = open('checkpoints/' + epoch + '.pb', 'rb')
    weights = pickle.load(A)
    A.close()
    print("Add checkpoints")
    mod.set_weights(weights)

    mod_names =  [i.name for i in mod.layers]
    indices = [i for i, s in enumerate(mod_names) if 'dense' in s or 'flatten' in s]
    indices.pop()
    for i in indices:
      new_model = Model(inputs=mod.input, outputs=mod.layers[i].output)
      print("Prediction")
      Y_pred = new_model.predict(input_matrix)
      df = pd.DataFrame(Y_pred)
      df.to_csv(mod.layers[i].name + '.' + epoch + '.tsv',  sep = "\t", index = False)

else:
  for i in range(len(model.layers) - 4, len(model.layers) - 1):
    new_model = Model(inputs=model.input, outputs=model.layers[i].output)
    Y_pred = new_model.predict(input_matrix)
    df = pd.DataFrame(Y_pred)
    df.to_csv(model.layers[i].name + '.tsv',  sep = "\t", index = False)
