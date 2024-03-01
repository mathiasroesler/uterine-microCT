#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# testProjectionPoints.py: Test function for the projection points algorithm
# Author: Mathias Roesler
# Last modified: 11/23

import sys
import scipy.io
import numpy as np
import utils.utils as utils
import thickness_analysis.plots as plots
import thickness_analysis.projection as projection


def findProjectionPointsTest():
	""" Tests the projection point algorithm on test images

	Arguments:

	Return:

	"""
	_dir = utils.HOME + '/' + utils.BASE + "/microCT/data/tests/"
	param_file = _dir + "test.toml"
	params = utils.parseTOML(param_file)

	for dataset in params["sets"]:
		test_dir = _dir + dataset + "/muscle_segmentation"
		img_stack = utils.loadImageStack(test_dir) # Load test images
		centreline_dict = scipy.io.loadmat(test_dir + "/centreline.mat")
		centreline = np.transpose(centreline_dict["centreline"])
		centreline = np.round(centreline).astype(int) # Round and convert to int

		for i in range(2):
			projection_points = projection.findProjectionPoints(
				img_stack[i], centreline[i], params["nb_points"], 
				params["horn"][i]) 
			plots.plotProjectionPoints(img_stack[i], 
				centreline[i, i*4:i*4+2], projection_points)


findProjectionPointsTest()
