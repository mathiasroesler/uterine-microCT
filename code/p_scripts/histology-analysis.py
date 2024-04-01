#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# histology-analysis.py: Script to generate the thickness based on histology
# Author: Mathias Roesler
# Last modified: 06/23

import argparse
import os
import pickle

import numpy as np
import scipy.io

import thickness_analysis.plots as plots
import thickness_analysis.projection as projection
import utils.utils as utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Specific script to generate the thickness plot "
        "based on histology"
    )

    parser.add_argument(
        "dir_path", type=str, metavar="dir-path",
        help="path from BASE to the dataset"
    )
    parser.add_argument(
        "base_name", type=str, metavar="base-name", help="name of the dataset"
    )
    parser.add_argument(
        "-e",
        "--extension",
        type=str,
        metavar="extension",
        help="extension for the saved images, default png",
        default="png",
    )
    parser.add_argument(
        "-p",
        "--points",
        type=int,
        help="number of points to use for the projection, default 128",
        default=128,
    )
    parser.add_argument(
        "-P",
        "--polar",
        action="store_true",
        help="flag used to plot the angular thickness in polar projection, "
        "default False",
    )

    # Parse input arguments
    args = parser.parse_args()

    load_directory = os.path.join(utils.HOME, utils.BASE, args.dir_path,
                                  args.base_name)

    param_file = os.path.join(load_directory, args.base_name + ".toml")
    params = utils.parseTOML(param_file)
    params = params["thickness"]  # Extract the thickness parameters

    # Add the muscle segmentation to the load directory
    load_directory = os.path.join(load_directory, "muscle_segmentation")

    # Dicts for results of both horns
    avg_thickness = dict()
    avg_slice_thickness = dict()
    errors = dict()

    # Use only the right horn
    horn = "right"

    print("Processing {} horn".format(horn))
    print("   Loading mask stack")
    mask_stack = utils.loadImageStack(
        os.path.join(load_directory, "{}".format(horn)),
        extension=args.extension
    )

    circular_win_size = round(0.04 * args.points)

    print("   Loading centreline")
    centreline_dict = scipy.io.loadmat(
        load_directory + "/{}/centreline.mat".format(horn)
    )
    centreline = np.transpose(centreline_dict["centreline"])

    print("   Estimating muscle thickness")
    muscle_thickness, slice_thickness, _ = projection.estimateMuscleThickness(
        mask_stack, centreline, args.points, params[horn]["slice_nbs"], horn
    )

    # Rescale the thickness to mm
    muscle_thickness *= params["scaling_factor"]
    slice_thickness *= params["scaling_factor"]

    print(
        "{} horn muscle thickness: {:.2f} \u00B1 {:.2f}".format(
            horn, np.mean(muscle_thickness), np.std(muscle_thickness)
        )
    )

    avg_slice_thickness[horn] = utils.circularAverage(
        slice_thickness, circular_win_size
    ).round(5)

    # Plot everything
    plots.plotAngularThickness(avg_slice_thickness, projection=args.polar)

    # Save angular thickness
    with open(load_directory + "/angular_thickness.pkl", "wb") as f:
        pickle.dump(avg_slice_thickness, f)
