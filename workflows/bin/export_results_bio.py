#! /usr/GenNet/env_GenNet/bin/python

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


from numpy import array, argmax
from tensorflow.keras.models import Model
from tensorflow.keras.optimizers import Adam
from sklearn.metrics import confusion_matrix, classification_report

sys.path.append('./')
import model
import params

name = sys.argv[1]


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

A = open('gene_mask.pb', 'rb')
gene_mask = pickle.load(A)
A.close()

A = open('checkpoints/15.pb', 'rb')
weights = pickle.load(A)
A.close()

num_classes = weights[len(weights) - 1].shape[0]

# Train model ####
mod = model.model_generator_train(num_classes, gene_mask, params)
opt = Adam(learning_rate = params.alpha)

mod.set_weights(weights)


with open(name + '_training_evaluation.tsv', 'w') as csv_file:
    writer = csv.writer(csv_file, delimiter = '\t', escapechar=' ', quoting = csv.QUOTE_NONE)
    for key, value in history.items():
       writer.writerow([key + '\t' + '\t'.join([str(item) for item in value ])])


Y_pred = mod.predict(x_test)
y_pred = np.argmax(Y_pred, axis=1)
y_class = np.argmax(y_test, axis=1)

cm = confusion_matrix(y_class, y_pred)
df = pd.DataFrame(cm, columns = proj_labels)
df.to_csv(name + '_confussionMatrix.tsv',  sep = "\t", index = False)

cr = classification_report(y_class, y_pred, target_names=proj_labels)
text_file = open(name + "_classificationReport.txt", "w")
text_file.write(cr)
text_file.close()

d = {'Real': labels[y_class], 'Predicted': labels[y_pred]}
df = pd.DataFrame(data = d )
df.to_csv(name + '_prediction.tsv',  sep = "\t", index = False)
