#! /usr/local/bin/python

#'#################################################################################
#'#################################################################################
#' Run random search on tcga data
#'#################################################################################
#'#################################################################################

import h5py
import pickle
import sys
import csv

from sklearn.model_selection import RandomizedSearchCV, train_test_split
from sklearn.preprocessing import LabelEncoder, OneHotEncoder

prop = float(sys.argv[1])
auto = sys.argv[2]


f = h5py.File('./assays.h5', 'r')
labels = f['label'][...]
methy = f['methy'][...]
f.close()

## embedding labels
# binary encode
onehot_encoder = OneHotEncoder(sparse=False)
integer_encoded = labels.reshape(len(labels), 1)
onehot_labels = onehot_encoder.fit_transform(integer_encoded)
index = range(len(labels))

x_train, x_test, y_train, y_test, index_train, index_test = train_test_split(methy, onehot_labels, index,
                                                    stratify = onehot_labels,
                                                    test_size = prop, random_state = 42)

if auto == "autoencoder":
    pickle.dump( [x_train, x_train], open( "train.pb", "wb" ), protocol = 4 )
    pickle.dump( [x_test, x_test], open( "test.pb", "wb" ), protocol = 4 )
else:
    pickle.dump( [x_train, y_train], open( "train.pb", "wb" ), protocol = 4 )
    pickle.dump( [x_test, y_test], open( "test.pb", "wb" ), protocol = 4 )


with open('test_indices.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows([index_test[index]] for index in range(len(index_test)))
