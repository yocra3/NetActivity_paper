#! /usr/local/bin/python

#'#################################################################################
#'#################################################################################
#'  Train a neural network based on kegg pathways
#'#################################################################################
#'#################################################################################

import pickle
import csv
import os
import numpy as np
import sys
import scipy
import functools
import operator
import pandas as pd
import h5py
import tensorflow as tf
import time

from numpy import array, argmax
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten, Input, Reshape
import tensorflow_model_optimization as tfmot


gpus = tf.config.list_physical_devices('GPU')
if gpus:
  try:
    # Currently, memory growth needs to be the same across GPUs
    for gpu in gpus:
      tf.config.experimental.set_memory_growth(gpu, True)
    logical_gpus = tf.config.list_logical_devices('GPU')
    print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
  except RuntimeError as e:
    # Memory growth must be set before GPUs have been initialized
    print(e)

sys.path.append('./')
import params
import model

name = sys.argv[1]
jobs = int(sys.argv[2])


f = h5py.File('assay_reshaped.h5', 'r')
proj_labels = f['label'].attrs['labels']
f.close()

A = open('train.pb', 'rb')
[x_train, y_train] = pickle.load(A)
A.close()

A = open('test.pb', 'rb')
[x_test, y_test] = pickle.load(A)
A.close()

## Load genes
with open('input_genes.txt','r') as file:
    genes = file.read()
genes = genes.split('\n')[0:-1]

## Select map rows with cpgs in input
with open('pathway_map.tsv','r') as file:
    gene_map = pd.read_csv(file, delimiter = '\t', index_col = False)


gene_map_filt = gene_map[gene_map['Symbol'].isin(genes)]
pathway_count = gene_map_filt['PathwayID'].value_counts()

gene_dict = dict(zip(genes, range(len(genes))))
gene_map_filt['idx'] = [gene_dict[gene] for gene in gene_map_filt['Symbol']]

## Create list of inputs - pathway
pathway_dict = {}
for g, idx in zip(gene_map_filt['PathwayID'].values, gene_map_filt['idx'].values):
  if g not in pathway_dict:
    pathway_dict[g] = [idx]
  else:
    pathway_dict[g].append(idx)

## Define gene mask ####
ngenes = len(genes)
npaths = len(pathway_dict.keys())

mask_d  = np.zeros((ngenes, npaths), float)

pcol = 0

pathways = pathway_dict.keys()
for path in pathways:
  mask_d[pathway_dict[path], pcol] = 1
  pcol += 1

df = pd.DataFrame(data = list(pathways) )
df.to_csv('pathways_names.txt',  sep = "\t", index = False)


# Define model ####
mod = model.model_generator_train(x_train, y_train, mask_d, params)
mod.summary()

# Pre-train wieghts in small model ####
callbacks = [tfmot.sparsity.keras.UpdatePruningStep()]

if params.epochs_prime:
    mask_log = mask_d == 1
    mini_weights = list()
    opt = Adam(learning_rate = params.alpha_prime)
    for i in range(npaths):
        print(i)
        minimodel = Sequential()
        if params.prepath_neurons:
            minimodel.add(Dense(params.prepath_neurons, input_dim = np.sum(mask_log[:, i])))
            minimodel.add(Dense(1))
        else:
            minimodel.add(Dense(1, input_dim = np.sum(mask_log[:, i])))
        minimodel.add(Dense(np.sum(mask_log[:, i])))
        minimodel.compile(loss='mse',
            optimizer = opt,
            metrics = ['mse'])
        minimodel.fit(x_train[:, mask_log[:, i]], x_train[:, mask_log[:, i]],
            batch_size = params.batch_size,
            epochs = params.epochs_prime,
            verbose = 1, validation_data = (x_test[:, mask_log[:, i]], x_test[:, mask_log[:, i]]), workers = jobs)
        w = minimodel.get_weights()
        if params.prepath_neurons:
            mini_weights.append([w[0], w[1], w[2], w[3]])
        else:
            mini_weights.append([w[0], w[1]])

    ## Convert weights to matrix
    if params.prepath_neurons:
        new_w1  = np.zeros((ngenes, npaths*params.prepath_neurons), float)
        new_b1  = np.zeros(npaths*params.prepath_neurons, float)
        new_w2  = np.zeros((npaths*params.prepath_neurons, npaths), float)
        new_b2  = np.zeros(npaths, float)
        for i in range(len(mini_weights)):
            for j in range(params.prepath_neurons):
                new_w1[mask_log[:, i], i*params.prepath_neurons + j] = mini_weights[i][0][:, j]
            new_b1[i*params.prepath_neurons:i*params.prepath_neurons + params.prepath_neurons] = mini_weights[i][1]
            new_w2[i*params.prepath_neurons:i*params.prepath_neurons + params.prepath_neurons, i:(i+1)] = mini_weights[i][2]
            new_b2[i] = mini_weights[i][3]
        ## Add modified weights
        w = mod.get_weights()
        w[0] = new_w1
        w[1] = new_b1
        w[5] = new_w2
        w[6] = new_b2
        mod.set_weights(w)
        if params.epochs_prime2:
            mod.layers[0].trainable = False
            mod.layers[1].trainable = False
            mod.compile(loss='mse',
                optimizer = opt,
                metrics = ['mse'])
            mod.fit(x_train, y_train,
              batch_size = params.batch_size,
              epochs = params.epochs_prime2,
              verbose = 1, callbacks = [callbacks], validation_data = (x_test, y_test), workers = jobs)
            mod.layers[0].trainable = True
            mod.layers[1].trainable = True


    else:
        new_w  = np.zeros((ngenes, npaths), float)
        new_b  = np.zeros(npaths, float)
        pcol = 0
        for i in range(len(mini_weights)):
            new_w[mask_log[:, i], i] = mini_weights[i][0][:, 0]
            new_b[i] = mini_weights[i][1]
            pcol += 1
        ## Add modified weights
        w = mod.get_weights()
        w[0] = new_w
        w[1] = new_b
        mod.set_weights(w)
        if params.epochs_prime2:
            mod.layers[0].trainable = False
            mod.compile(loss='mse',
                optimizer = opt,
                metrics = ['mse'])
            mod.fit(x_train, y_train,
                batch_size = params.batch_size,
                epochs = params.epochs_prime2,
                verbose = 1, callbacks = [callbacks], validation_data = (x_test, y_test), workers = jobs)
            mod.layers[0].trainable = True

# Train full model ####
opt = Adam(learning_rate = params.alpha)
mod.compile(loss='mse',
        optimizer = opt,
        metrics = ['mse'])

history = mod.fit(x_train, y_train,
  batch_size = params.batch_size,
  epochs = params.epochs,
  verbose = 1, callbacks = [callbacks], validation_data = (x_test, y_test), workers = jobs)

history_dict = history.history
pickle.dump( history_dict, open( name + "_history_model.pb", "wb" ), protocol = 4 )
pickle.dump( proj_labels, open( name + "_labels.pb", "wb" ), protocol = 4 )
mod.save( name )
