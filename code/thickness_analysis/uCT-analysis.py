#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# uCT-analysis.py: Script to analyse muscle thickness in uterine horns
# Author: Mathias Roesler
# Last modified: 02/23

import os
import utils
import plots
import argparse
import projection
import numpy as np


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=
		"Determines the thickness of the muscle layers from a uCT dataset")

	parser.add_argument("input_folder", type=str, metavar="input-folder",
		help="name of the folder containing the images")
	parser.add_argument("-e", "--extension", type=str, metavar="extension",
		help="extension for the saved images", default="png")
	parser.add_argument("--horn", type=str, choices={"left", "right", "both"},
		help="horn to process", default="both")
	parser.add_argument("-p", "--points", type=int, metavar="points",
		help="number of points to use for the projection", default=32)
	parser.add_argument("-s", "--switch", action='store_true',
		help="switches the labels of the left and right horn")

	# Parse input arguments
	args = parser.parse_args()

	full_path = os.path.join(utils.HOME, utils.BASE, args.input_folder, 
		utils.MASK_FOLDER)

	param_file = full_path + "/analysis.toml"
	params = utils.parseTOML(param_file)

	if args.horn == "both":
		horns = ["left", "right"]

	else:
		horns = [args.horn]

	# Dicts for results of both horns
	avg_thickness = dict()
	avg_slice_thickness = dict()
	errors = dict()

	for i, horn in enumerate(horns):
		print("Processing {} horn".format(horn))
		print("   Loading mask stack")
		mask_stack = utils.loadImageStack(os.path.join(
			full_path, "{}_horn".format(horn)), extension=args.extension)

		param_dict = params[horn]
		short_mask_stack = mask_stack[
			param_dict['rot_start']:param_dict['rot_end']]
		nb_imgs = len(short_mask_stack)
		slice_nbs = [30, nb_imgs // 2, nb_imgs-52]

		print("   Finding centreline")
		centreline = projection.findCenterline(short_mask_stack, horn=horn)
		
		print("   Estimating muscle thickness")
		muscle_thickness, slice_thickness = projection.estimateMuscleThickness(
			short_mask_stack, centreline, args.points, slice_nbs)  

		# Rescale the thickness to mm
		muscle_thickness *= params["scaling_factor"]
		slice_thickness *= params["scaling_factor"]

		print(u"{} horn muscle thickness: {:.2f} \u00B1 {:.2f}".format(horn, 
			np.mean(muscle_thickness), np.std(muscle_thickness)))
		
		if args.switch:
			avg_thickness[horns[i-1]] = utils.movingAverage(
				muscle_thickness, 15)
			avg_slice_thickness[horns[i-1]] = utils.circularAverage(
				slice_thickness, 9)
			errors[horns[i-1]] = utils.movingStd(muscle_thickness, 15)

		else:
			avg_thickness[horn] = utils.movingAverage(muscle_thickness, 15)
			avg_slice_thickness[horn] = utils.circularAverage(slice_thickness, 9)
			errors[horn] = utils.movingStd(muscle_thickness, 15)

	# Plot everything
	plots.plotMuscleThickness(avg_thickness, errors)
	plots.plotAngularThickness(avg_slice_thickness)
