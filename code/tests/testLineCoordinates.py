#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# testLineCoordinates.py: Test function for the line coordinates algorithm
# Author: Mathias Roesler
# Last modified: 11/23

import sys
import scipy.io
import numpy as np
import utils.utils as utils
import matplotlib.pyplot as plt
import thickness_analysis.projection as projection

_dir = utils.HOME + '/' + utils.BASE + "/microCT/data/tests/"
_nb_points = 10
_sets = ["AWA015", "AWA030"]
_horn = ["left", "right"]
_theta = [np.pi/4, np.pi/2, 3*np.pi/4, np.pi]


def lineCoordinatesTest():
	""" Tests the line coordinate algorithm on test images

	Arguments:

	Return:

	"""
	cnt = 0 # Counter for theta

	for dataset in _sets:
		test_dir = _dir + dataset + "/muscle_segmentation"
		img_stack = utils.loadImageStack(test_dir) # Load test images
		centreline_dict = scipy.io.loadmat(test_dir + "/centreline.mat")
		centreline = np.transpose(centreline_dict["centreline"])
		centreline = np.round(centreline).astype(int) # Round and convert to int

		for i in range(2):
			if _horn[i] == "left":
				centre_point = centreline[i, 0:2]

			else:
				centre_point = centreline[i, 4:6]

			img = img_stack[i]
			line_x, line_y = projection.findLineCoordinates(img.shape, 
				centre_point, _theta[cnt])

			cnt += 1 # Increment counter

			plt.imshow(img, cmap='gray') 
			plt.plot(line_x, line_y)
			plt.plot(centre_point[0], centre_point[1], '.r')
			plt.show()


lineCoordinatesTest()

