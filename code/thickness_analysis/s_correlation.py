#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# s_correlation.py: Script to estimate correlation between histology and uCT
# Author: Mathias Roesler
# Last modified: 06/23

import os
import sys
import utils
import argparse
import scipy.io
import projection
import numpy as np


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=
		"Specific script to estimate the correlation between histology and uCT")

	parser.add_argument("uCT_folder", type=str, metavar="uCT-folder",
		help="name of the folder containing the uCT data")
	parser.add_argument("histo_folder", type=str, metavar="histo-folder",
		help="name of the folder containing the histology data")
	parser.add_argument("--horn", type=str, choices={"left", "right"},
		help="horn to process", default="right")

	# Parse input arguments
	args = parser.parse_args()

	# Set up variables
	uCT_path = os.path.join(utils.HOME, utils.BASE, args.uCT_folder, 
		utils.DATA_FOLDER)
	histo_path = os.path.join(utils.HOME, utils.BASE, args.histo_folder, 
		"muscle_segmentation")
	regions = ["cervical", "central", "ovarian"]
	horn = args.horn

	# Read data 
	uCT_data = np.load(uCT_path + "/angular_thickness.pkl", 
		allow_pickle=True)[horn]
	histo_data = np.load(histo_path + "/angular_thickness.pkl",
		 allow_pickle=True)[horn]

	try:
		assert(uCT_data.shape == histo_data.shape)

	except AssertionError:
		sys.stderr.write("Error: uCT and histology data show have same shape.\n")

	breakpoint()
	for i in range(uCT_data.shape[1]):
		correl_matrix = np.corrcoef(uCT_data[:, i], histo_data[:, i])
		print(u"{} section correlation factor: {:.2f}".format(
			regions[i].capitalize(), correl_matrix[0, 1]))
