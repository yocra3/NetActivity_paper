#! /usr/local/bin/python

#'#################################################################################
#'#################################################################################
#'  Extract features from network
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
from tensorflow.keras.models import Sequential, Model

f = h5py.File('./assay_reshaped.h5', 'r')
methy = f['methy'][...]
f.close()

model = tf.keras.models.load_model('model/' + os.listdir('model/')[0])

mod_names =  [i.name for i in model.layers]
indices = [i for i, s in enumerate(mod_names) if 'dense' in s or 'flatten' in s]
indices.pop()
for i in indices:
  new_model = Model(inputs=model.input, outputs=model.layers[i].output)
  Y_pred = new_model.predict(methy)
  df = pd.DataFrame(Y_pred)
  df.to_csv(model.layers[i].name + '.tsv',  sep = "\t", index = False)
