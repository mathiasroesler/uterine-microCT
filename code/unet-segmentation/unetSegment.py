#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# unetSegmentation.py: Script to segment microCT using unet model
# Author: Mathias Roesler
# Last modified: 04/24

import os
import argparse

import numpy as np
import utils.utils as utils
import skimage.io as skio

import tf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Segment microCT dataset with unet model"
    )

    parser.add_argument(
        "dir_path", type=str, metavar="dir-path",
        help="path from BASE to the dataset"
    )
    parser.add_argument(
        "base_name", type=str, metavar="base-name", help="name of the dataset"
    )
    parser.add_argument(
        "model", type=str, metavar="model",
        help="Saved model to use"
    )
    parser.add_argument(
        "-e",
        "--extension",
        type=str,
        metavar="extension",
        help="extension for the saved images, default png",
        default="png",
    )

    # Parse input arguments
    args = parser.parse_args()

    load_directory = os.path.join(utils.HOME, utils.BASE, args.dir_path,
                                  args.base_name)
    save_directory = load_directory + "/masks"

    imgs = utils.loadImageStack(load_directory + "/imgs")

    # Convert to floats between 0 and 1
    imgs = np.asarray(imgs, dtype=np.float32) / imgs.max()

    # Ensure shape are correct
    if len(imgs.shape) < 4:
        imgs = imgs.reshape(imgs.shape[0],
                            imgs.shape[1],
                            imgs.shape[2],
                            1)

    # Load trained model
    model = tf.keras.models.load(args.model)

    for i, img in enumerate(imgs):
        mask = model.predict(img)
        skio.imsave("{}/{}_{}.{}".format(
            save_directory, args.base_name, i, args.extension
        ), mask, check_contrast=False)
