from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

# Deterministic behavior
import os
os.environ['TF_CUDNN_USE_AUTOTUNE'] = '0'
prng_seed = 1009732533  # OEIS A002205
import numpy
numpy.random.seed(prng_seed)
import tensorflow
tensorflow.set_random_seed(prng_seed)

batch_size = 128
nb_classes = 2
nb_epoch = 12

ndense = 512
dropout = 0.1
nlayer = 4

model_loss = 'categorical_crossentropy'
model_optimizer = 'adadelta'

import sys, getopt

try:
    option, argument = getopt.getopt(
        sys.argv[1:], '',
        ['batch-size=', 'nb-epoch=', 'ndense=', 'dropout=',
         'nlayer=', 'loss=', 'optimizer='])
except getopt.GetoptError:
    sys.exit()
print(option, file = sys.stderr)
for o, a in option:
    if o == '--batch-size':
        batch_size = int(a)
    elif o == '--nb-epoch':
        nb_epoch = int(a)
    elif o == '--ndense':
        ndense = int(a)
    elif o == '--dropout':
        dropout = float(a)
    elif o == '--nlayer':
        nlayer = int(a)
    elif o == '--loss':
        model_loss = a
    elif o == '--optimizer':
        model_optimizer = a

# the data, shuffled and split between train and test sets
#(X_train, y_train), (X_test, y_test) = mnist.load_data()
import ml_out
(X_train, y_train), (X_test, y_test) \
    = (ml_out.X_train, ml_out.y_train), \
    (ml_out.X_test, ml_out.y_test)
# Transform 1 vs. 2 photons into [0, 1]
y_train -= 1
y_test -= 1

nfeature = X_train.shape[1]

import keras.backend

if keras.backend.image_dim_ordering() == 'th':
    X_train = X_train.reshape(X_train.shape[0], nfeature)
    X_test = X_test.reshape(X_test.shape[0], nfeature)
else:
    X_train = X_train.reshape(X_train.shape[0], nfeature)
    X_test = X_test.reshape(X_test.shape[0], nfeature)

mean = X_train.mean(0)
X_train -= mean
X_test -= mean
std = X_train.std(0)
X_train /= std
X_test /= std
#print(list(X_train.mean(0)))
#sys.exit(0)

X_train = X_train.astype('float32')
X_test = X_test.astype('float32')
print('X_train shape:', X_train.shape, file = sys.stderr)
print(X_train.shape[0], 'train samples', file = sys.stderr)
print(X_test.shape[0], 'test samples', file = sys.stderr)

import keras.models
import keras.layers
import keras.utils

# convert class vectors to binary class matrices
Y_train = keras.utils.np_utils.to_categorical(y_train, nb_classes)
Y_test = keras.utils.np_utils.to_categorical(y_test, nb_classes)

model = keras.models.Sequential()

model.add(keras.layers.Dense(ndense, input_shape=(nfeature,)))
model.add(keras.layers.Activation('relu'))
model.add(keras.layers.Dropout(dropout))
for i in range(nlayer):
    model.add(keras.layers.Dense(ndense))
    model.add(keras.layers.Activation('relu'))
    model.add(keras.layers.Dropout(dropout))
model.add(keras.layers.Dense(2))
model.add(keras.layers.Activation('softmax'))

model.compile(loss=model_loss,
              optimizer=model_optimizer,
              metrics=['accuracy'])

model.fit(X_train, Y_train, batch_size=batch_size, nb_epoch=nb_epoch,
          verbose=1, validation_data=(X_test, Y_test))
score = model.evaluate(X_test, Y_test, verbose=0)
print('Test score:', score[0], file = sys.stderr)
print('Test accuracy:', score[1], file = sys.stderr)

f = open('ml_predict.py', 'w')
f.write('y_predict = ' +
        str(map(list, list(model.predict(X_test)))) + '\n')
f.close()
