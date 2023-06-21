#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# histology-analysis.py: Script to generate the thickness based on histology
# Author: Mathias Roesler
# Last modified: 06/23

import os
import pickle
import scipy.io
import argparse
import numpy as np
import utils.utils as utils
import thickness_analysis.plots as plots
import thickness_analysis.projection as projection


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=
		"Specific script to generate the thickness plot based on histology")

	parser.add_argument("dir_path", type=str, metavar="dir-path",
		help="path from BASE to the dataset")
	parser.add_argument("base_name", type=str, metavar="base-name",
		help="name of the dataset")
	parser.add_argument("-e", "--extension", type=str, metavar="extension",
		help="extension for the saved images", default="png")
	parser.add_argument("--points", type=int, 
		help="number of points to use for the projection", default=128)

	# Parse input arguments
	args = parser.parse_args()
	
	load_directory = os.path.join(utils.HOME, utils.BASE, args.dir_path, 
		args.base_name)

	param_file = os.path.join(load_directory, args.base_name + ".toml")
	params = utils.parseTOML(param_file)
	params = params['thickness'] # Extract the thickness parameters

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
	mask_stack = utils.loadImageStack(os.path.join(
		load_directory, "{}".format(horn)), extension=args.extension)

	circular_win_size = round(0.04 * args.points)

	print("   Finding centreline")
	centreline_dict = scipy.io.loadmat(load_directory + 
		"/{}/centreline.mat".format(horn))
	centreline = np.transpose(centreline_dict["centreline"])

	print("   Estimating muscle thickness")
	muscle_thickness, slice_thickness = projection.estimateMuscleThickness(
		mask_stack, centreline, args.points, params[horn]["slice_nbs"], horn)  

	# Rescale the thickness to mm
	muscle_thickness *= params["scaling_factor"]
	slice_thickness *= params["scaling_factor"]

	print(u"{} horn muscle thickness: {:.2f} \u00B1 {:.2f}".format(horn, 
		np.mean(muscle_thickness), np.std(muscle_thickness)))
	
	avg_slice_thickness[horn] = utils.circularAverage(slice_thickness,
		circular_win_size)

	# Plot everything
	plots.plotAngularThickness(avg_slice_thickness)


	# Save angular thickness
	with open(load_directory + "/angular_thickness.pkl".format(
		horn), 'wb') as f:
		pickle.dump(avg_slice_thickness, f)
