#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Programs that converts a stack of images into a
# nifti volume. Requires SimpleITK package.
# Author: Mathias Roesler
# Last modified: 12/22

import os
import sys
import glob
import argparse
import SimpleITK as sitk

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=
		"Converts uCT dataset into the nifti format")

	parser.add_argument("input_folder", type=str, metavar="input-folder",
		help="name of the folder containing the images")
	parser.add_argument("save_name", type=str, metavar="save-name",
		help="name of the nifti file to save to")
	parser.add_argument("-e", "--extension", type=str, 
		help="extension of the images", default="png")

	# Parse input arguments
	args = parser.parse_args()

	input_folder_path = os.path.join(os.path.expanduser('~'), "Documents/phd",
		 args.input_folder)

	img_path = os.path.join(input_folder_path, "*." + args.extension)

	if not os.path.exists(input_folder_path):
		sys.stderr.write("\nError: the input folder {} does not exist \
			\n".format(args.input_folder))
		exit()

	img_list = sorted(glob.glob(img_path))

	if len(img_list) == 0:
		sys.stderr.write("Error: the folder {} does not contain any images \
			\n".format(args.input_folder))
		exit()

	if len(args.save_name.split('.')) == 1:
		# Add the .nii.gz extension if it does not exist
		save_name = args.save_name + ".nii.gz"	

	else:
		save_name = args.save_name

	reader = sitk.ImageSeriesReader()
	reader.SetFileNames(img_list)
	volume = reader.Execute()

	sitk.WriteImage(volume, os.path.join(input_folder_path, save_name))
