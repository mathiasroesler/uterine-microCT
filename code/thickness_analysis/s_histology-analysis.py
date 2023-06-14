#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# s_histology-analysis.py: Script to generate the thickness based on histology
# Author: Mathias Roesler
# Last modified: 03/23

import os
import utils
import plots
import pickle
import scipy.io
import argparse
import projection
import numpy as np


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=
		"Specific script to generate the thickness plot based on histology")

	parser.add_argument("input_folder", type=str, metavar="input-folder",
		help="name of the folder containing the images")
	parser.add_argument("-e", "--extension", type=str, metavar="extension",
		help="extension for the saved images", default="png")
	parser.add_argument("--points", type=int, 
		help="number of points to use for the projection", default=128)

	# Parse input arguments
	args = parser.parse_args()
	
	full_path = os.path.join(utils.HOME, utils.BASE, args.input_folder, 
		"muscle_segmentation")

	param_file = full_path + "/analysis.toml"
	params = utils.parseTOML(param_file)

	# Dicts for results of both horns
	avg_thickness = dict()
	avg_slice_thickness = dict()
	errors = dict()

	# Use only the right horn
	horn = "right"

	print("Processing {} horn".format(horn))
	print("   Loading mask stack")
	mask_stack = utils.loadImageStack(os.path.join(
		full_path, "{}".format(horn)), extension=args.extension)

	circular_win_size = round(0.05 * args.points)

	print("   Finding centreline")
	centreline_dict = scipy.io.loadmat(full_path + 
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
	with open(full_path + "/{}/angular_thickness.pkl".format(
		horn), 'wb') as f:
		pickle.dump(avg_slice_thickness, f)
