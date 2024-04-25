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

from skimage.util import img_as_ubyte

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

    # Check directory exists
    if not os.path.isdir(save_directory):
        os.mkdir(save_directory)

    img_list = sorted(glob.glob("*.{}".format(args.extension),
                                root_dir=load_directory))

    for img_name in img_list:
        path = os.path.join(load_directory, img_name)

        if not os.path.isfile(path):
            sys.stderr.write("Error: {} is not a file.\n".format(path))
            exit()

        img = skio.imread(path, as_gray=True)
        height, width = img.shape

        if height % new_size:
            # Height padding
            pad_h = new_size * ((height // new_size) + 1) - height
        else:
            pad_h = 0

        if width % new_size:
            # Width padding
            pad_w = new_size * ((width // new_size) + 1) - width
        else:
            pad_w = 0

        # Pad image
        padded_img = np.pad(img, ((0, pad_h), (0, pad_w)))

        for i in range((height + pad_h) // new_size):
            for j in range((width + pad_w) // new_size):
                # Get each image block
                img_block = padded_img[
                    i*new_size:(i+1)*new_size,
                    j*new_size:(j+1)*new_size
                ]

                # Save the image block
                skio.imsave("{}/{}_block_{}_{}.{}".format(
                    save_directory, img_name, i, j, args.extension
                ), img_as_ubyte(img_block), check_contrast=False)