#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# testCentrePoints.py: Test function for the centre point location
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


def centrePointsTest():
	""" Tests the location of the centre points in test images

	Arguments:

	Return:

	"""
	for dataset in _sets:
		test_dir = _dir + dataset + "/muscle_segmentation"
		img_stack = utils.loadImageStack(test_dir) # Load test images
		centreline_dict = scipy.io.loadmat(test_dir + "/centreline.mat")
		centreline = np.transpose(centreline_dict["centreline"])
		centreline = np.round(centreline).astype(int) # Round and convert to int

		for i in range(2):
			centre_point = centreline[i]
			img = img_stack[i]

			plt.imshow(img, cmap='gray') 

			for j in range(3):
				plt.plot(centre_point[j*2], centre_point[j*2+1], '.r')

			plt.show()


centrePointsTest()
