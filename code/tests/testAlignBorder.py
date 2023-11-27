#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# testAlignBorder.py: Test function for the border alignment alogrithm
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


def alignBorderTest():
	""" Tests the alignment alogrithm

	Arguments:

	Return:

	"""
	for dataset in _sets:
		test_dir = _dir + dataset + "/muscle_segmentation"
		img_stack = utils.loadImageStack(test_dir) # Load test images
		centreline_dict = scipy.io.loadmat(test_dir + "/centreline.mat")
		centreline = np.transpose(centreline_dict["centreline"])
		centreline = np.round(centreline).astype(int) # Round and convert to int

		for i, img in enumerate(img_stack):
			projection_points = projection.findProjectionPoints(img, centreline[i, :], 
				_nb_points, _horn[i])
			diff = np.diff(projection_points, axis=0)
			norm = np.linalg.norm(diff, axis=1)	
			thickness = norm[np.arange(0, projection_points.shape[0], 2)]

			# Find the two halves 
			indices = np.arange(len(projection_points))
			right_half = indices[(indices % 4 == 0) | ((indices - 1) % 4 == 0)]
			left_half = np.setdiff1d(indices, right_half)

			# Plot first half
			plt.imshow(img, cmap='gray') 
			plt.plot(projection_points[right_half][:, 0], 
				projection_points[right_half][:, 1], '.b')

			# Plot second half
			plt.plot(projection_points[left_half][:, 0], 
				projection_points[left_half][:, 1], '.g')
			plt.show()


alignBorderTest()

