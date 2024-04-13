#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# cropImages.py: Script to resize images to for unet model
# Author: Mathias Roesler
# Last modified: 04/24

import argparse
import glob
import sys
import os

import numpy as np
import skimage.io as skio
import utils.utils as utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Resizes images and masks to be square"
    )

    parser.add_argument(
        "dir_path", type=str, metavar="dir-path",
        help="path from BASE to the dataset"
    )
    parser.add_argument(
        "base_name", type=str, metavar="base-name", help="name of the dataset"
    )
    parser.add_argument(
        "new_size", type=int, metavar="new-size",
        help="size of the new training image"
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

    new_size = args.new_size
    load_directory = os.path.join(utils.HOME, utils.BASE, args.dir_path,
                                  args.base_name)
    save_directory = os.path.join(load_directory, "resized/")

    img_list = sorted(glob.glob("*.{}".format(args.extension),
                                root_dir=load_directory))

    for img_name in img_list:
        path = os.path.join(load_directory, img_name)

        if not os.path.isfile(path):
            sys.stderr.write("Error: {} is not a file.\n".format(path))
            exit()

        img = skio.imread(path, as_gray=True)
        height, width = img.shape

        if height <= new_size and width <= new_size:
            # Make sure image block is the correct size
            img_block = np.pad(img, ((0, new_size-img.shape[0]),
                                     (0, new_size-img.shape[1])
                                     )
                               )

            skio.imsave("{}/{}_block_{}_{}.{}".format(
                save_directory, img_name, 0, 0, args.extension
            ), img_block, check_contrast=False)

        else:
            img_block = img[:512, :512]
            skio.imsave("{}/{}_block_{}_{}.{}".format(
                save_directory, img_name, 0, 0, args.extension
            ), img_block, check_contrast=False)
