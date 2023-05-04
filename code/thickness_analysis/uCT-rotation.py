#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# uCT-rotation.py: Script to rotate the uterine uCT dataset
# Author: Mathias Roesler
# Last modified: 02/23

import os
import utils
import argparse
import projection
import numpy as np


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=
		"Rotates a uCT dataset of the uterus along both horns")

	parser.add_argument("input_folder", type=str, metavar="input-folder",
		help="name of the folder containing the images")
	parser.add_argument("-e", "--extension", type=str, metavar="extension",
		help="extension for the saved images", default="png")
	parser.add_argument("--horn", type=str, choices={"left", "right", "both"},
		help="horn to process", default="both")
	parser.add_argument("--padding", type=int, metavar="padding",
		help="padding value for the rotation", default=300)

	# Parse input arguments
	args = parser.parse_args()

	full_path = os.path.join(utils.HOME, utils.BASE, args.input_folder, 
		utils.MASK_FOLDER)

	param_file = full_path + "/analysis.toml"
	params = utils.parseTOML(param_file)

	print("Loading mask stack")
	mask_stack = utils.loadImageStack(full_path, extension=args.extension)

	if args.horn == "both":
		horns = ["left", "right"]

	else:
		horns = [args.horn]

	for horn in horns:
		print("Processing {} horn".format(horn))
		param_dict = params[horn]

		print("   Finding centreline")
		centreline = projection.findCenterline(
			mask_stack[params["start_nb"]:param_dict["end"]], horn=horn)

		origin = np.array([params["start_nb"], 
			centreline[0][0],
			centreline[0][1]])

		tip = np.array([centreline.shape[0],
			 centreline[-1][0], 
			centreline[-1][1]])

		theta = np.arccos(np.dot(tip, origin) / (
			np.linalg.norm(tip) * np.linalg.norm(origin)))

		if horn == "right":
			# Rotate in the other direction for the right horn
			theta = -theta

		print("   Rotating stack")
		rotated_mask_stack = projection.rotateImageStack(mask_stack, theta, 
			args.padding)

		print("   Saving stack")
		utils.saveImageStack(rotated_mask_stack, 
			os.path.join(full_path, "{}_horn".format(horn)), 
			params['prefix'], extension=args.extension)
