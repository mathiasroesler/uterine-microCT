#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# unetSegmentation.py: Script to train the unet model
# Author: Mathias Roesler
# Last modified: 04/24

import argparse
import os

import numpy as np
import utils.utils as utils

from sklearn.model_selection import train_test_split

from keras_unet.models import custom_unet
from keras_unet.metrics import iou, iou_thresholded
from keras.optimizers import Adam
from keras.callbacks import ModelCheckpoint

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Trains the unet model"
    )

    parser.add_argument(
        "dir_path", type=str, metavar="dir-path",
        help="path from BASE to the dataset"
    )
    parser.add_argument(
        "base_name", type=str, metavar="base-name", help="name of the dataset"
    )
    parser.add_argument(
        "epochs", type=int, help="number of epochs to run"
    )
    parser.add_argument(
        "steps", type=int, help="number of steps per epoch"
    )

    # Parse input arguments
    args = parser.parse_args()

    load_directory = os.path.join(utils.HOME, utils.BASE, args.dir_path,
                                  args.base_name)

    imgs = utils.loadImageStack(load_directory + "/imgs")
    masks = utils.loadImageStack(load_directory + "/masks")

    # Convert to floats between 0 and 1
    imgs = np.asarray(imgs, dtype=np.float32) / imgs.max()
    masks = np.asarray(masks, dtype=np.float32) / masks.max()

    # Ensure shape are correct
    if len(imgs.shape) < 4:
        imgs = imgs.reshape(imgs.shape[0],
                            imgs.shape[1],
                            imgs.shape[2],
                            1)
    if len(masks.shape) < 4:
        masks = masks.reshape(masks.shape[0],
                              masks.shape[1],
                              masks.shape[2],
                              1)

imgs_val, imgs_train, masks_val, masks_train = train_test_split(
    imgs, masks, test_size=0.9, random_state=0)

print("imgs_train: ", imgs_train.shape)
print("masks_train: ", masks_train.shape)
print("imgs_val: ", imgs_val.shape)
print("masks_val: ", masks_val.shape)

input_shape = imgs_train[0].shape

model = custom_unet(
    input_shape,
    filters=32,
    use_batch_norm=True,
    dropout=0.3,
    dropout_change_per_layer=0.0,
    num_layers=4
)


model_filename = 'unet-model.keras'
callback_checkpoint = ModelCheckpoint(
    model_filename,
    verbose=1,
    monitor='val_loss',
    save_best_only=True,
)


model.compile(
    optimizer=Adam(),
    loss='binary_crossentropy',
    metrics=[iou, iou_thresholded]
)


history = model.fit(
    imgs_train, masks_train,
    steps_per_epoch=args.steps,
    epochs=args.epochs,
    validation_data=(imgs_val, masks_val),
    callbacks=[callback_checkpoint]
)
