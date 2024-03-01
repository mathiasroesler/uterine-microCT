#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# testSeparateHorns.py: Test function for the horn separation algorithm
# Author: Mathias Roesler
# Last modified: 11/23

import sys
import scipy.io
import numpy as np
import utils.utils as utils
import matplotlib.pyplot as plt
import thickness_analysis.projection as projection


def separateHornsTest():
	""" Tests the horn separation algorithm on test images

	Arguments:

	Return:

	"""
	_dir = utils.HOME + '/' + utils.BASE + "/microCT/data/tests/"
	param_file = _dir + "test.toml"
	params = utils.parseTOML(param_file)

	for dataset in params["sets"]:
		print("Testing set {}".format(dataset))
		test_dir = _dir + dataset + "/muscle_segmentation"
		img_stack = utils.loadImageStack(test_dir) # Load test images
		centreline_dict = scipy.io.loadmat(test_dir + "/centreline.mat")
		centreline = np.transpose(centreline_dict["centreline"])
		centreline = np.round(centreline).astype(int) # Round and convert to int

		for i in range(2):
			centre_point = centreline[i]
			img = img_stack[i]

			if (centre_point[2:4] != np.array([0, 0])).all():
				# The horns are not clearly separated and three points are given
				# Create a vector between the left and right points

				n = np.array([
					centre_point[1] - centre_point[5], 
					centre_point[4] - centre_point[0]])
				n = n / np.linalg.norm(n)

				line_x, line_y = projection.separateHorns(
					img.shape, centre_point[2:4], n)

				for j in range(len(line_x)):
					# Clear half of the image based on the horn
					if params["horn"][i] == "left":
						img[j, line_x[j]:] = 0

					elif params["horn"][i] == "right":
						img[j, :line_x[j]] = 0;

				# Plot
				plt.imshow(img, cmap='gray') 
				plt.plot(line_x, line_y)
				plt.plot(centre_point[2], centre_point[3], '.r')
				plt.show()


separateHornsTest()
