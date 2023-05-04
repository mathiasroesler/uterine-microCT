#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# s_projection.py: Script to generate the projection point figure
# Author: Mathias Roesler
# Last modified: 03/23

import os
import utils
import plots
import argparse
import projection
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=
		"Specific script to generate the projection point plot")

	parser.add_argument("input_folder", type=str, metavar="input-folder",
		help="name of the folder containing the images")
	parser.add_argument("img_nb", type=str, metavar="img-nb", 
		help="number of the image to use")
	parser.add_argument("-e", "--extension", type=str, metavar="extension",
		help="extension for the saved images", default="png")
	parser.add_argument("--horn", type=str, choices={"left", "right"},
		help="horn to process", default="right")
	parser.add_argument("-p", "--points", type=int, metavar="points",
		help="number of points to use for the projection", default=8)

	# Parse input arguments
	args = parser.parse_args()
	
	full_path = os.path.join(utils.HOME, utils.BASE, args.input_folder, 
		utils.MASK_FOLDER)

	param_file = full_path + "/analysis.toml"
	params = utils.parseTOML(param_file)

	img_file = full_path + "/{}_horn/".format(args.horn) + params["prefix"] + args.img_nb + "." + args.extension
	img = plt.imread(img_file)

	# 3D array needed for it to work
	centreline = projection.findCenterline(np.array([img, img]), 
		horn=args.horn)
	projection_points = projection.findProjectionPoints(img, centreline[0], args.points)

	plots.plotProjectionPoints(img, centreline[0], projection_points)
