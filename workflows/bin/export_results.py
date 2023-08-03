#! /usr/local/bin/python

#'#################################################################################
#'#################################################################################
#' Export tcga model results
#'#################################################################################
#'#################################################################################


import pickle
import os
import csv
import numpy as np
import sys
import pandas as pd
import h5py
import tensorflow as tf


from tensorflow.keras.models import Sequential
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.models import Sequential, Model
from keras.layers import Conv1D, MaxPooling1D, Dense, Dropout, Activation, Flatten, Input
from keras.wrappers.scikit_learn import KerasClassifier
from keras.callbacks import EarlyStopping, ModelCheckpoint
from sklearn.metrics import confusion_matrix, classification_report

name = sys.argv[1]
auto = sys.argv[2]

model = tf.keras.models.load_model('model/' + os.listdir('model/')[0])

f = h5py.File('assay_reshaped.h5', 'r')
proj_labels = f['label'].attrs['labels']
f.close()


A = open('history_model.pb', 'rb')
history = pickle.load(A)
A.close()

A = open('labels.pb', 'rb')
labels = pickle.load(A)
A.close()

A = open('test.pb', 'rb')
[x_test, y_test] = pickle.load(A)
A.close()

with open(name + '_training_evaluation.txt', 'w') as csv_file:
    writer = csv.writer(csv_file, delimiter = '\t', escapechar=' ', quoting = csv.QUOTE_NONE)
    for key, value in history.items():
       writer.writerow([key + '\t' + '\t'.join([str(item) for item in value ])])

if auto != "autoencoder":
    Y_pred = model.predict(x_test)
    y_pred = np.argmax(Y_pred, axis=1)
    y_class = np.argmax(y_test, axis=1)

    cm = confusion_matrix(y_class, y_pred)
    df = pd.DataFrame(cm, columns = proj_labels)
    df.to_csv(name + '_confussionMatrix.txt',  sep = "\t", index = False)

    cr = classification_report(y_class, y_pred, target_names=proj_labels)
    text_file = open(name + "_classificationReport.txt", "w")
    text_file.write(cr)
    text_file.close()

    d = {'Real': labels[y_class], 'Predicted': labels[y_pred]}
    df = pd.DataFrame(data = d )
    df.to_csv(name + '_prediction.txt',  sep = "\t", index = False)
