#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# testExcludePoints.py: Test function for excluding central points
# Author: Mathias Roesler
# Last modified: 04/24

import matplotlib.pyplot as plt
import numpy as np
import scipy.io

import thickness_analysis.projection as projection
import utils.utils as utils


def excludeCentralPointsTest():
    """Tests the exclusion algorithm

    Arguments:

    Return:

    """
    _dir = utils.HOME + "/" + utils.BASE + "/microCT/data/tests/"
    param_file = _dir + "test.toml"
    params = utils.parseTOML(param_file)

    for dataset in params["sets"]:
        print("Testing set {}".format(dataset))
        test_dir = _dir + dataset + "/muscle_segmentation"
        img_stack = utils.loadImageStack(test_dir)  # Load test images
        centreline_dict = scipy.io.loadmat(test_dir + "/centreline.mat")
        centreline = np.transpose(centreline_dict["centreline"])
        centreline = np.round(centreline).astype(int)  # Convert to int

        for i, img in enumerate(img_stack):
            projection_points = projection.findProjectionPoints(
                img, centreline[i, :], params["nb_points"], params["horn"][i]
            )

            if (centreline[i, 2:4] != np.array([0, 0])).all():
                projection_points = projection.excludeCentralPoints(
                    img, centreline[i, :], projection_points,
                    params["horn"][i])

            diff = np.diff(projection_points, axis=0)
            norm = np.linalg.norm(diff, axis=1)
            thickness = norm[np.arange(0, projection_points.shape[0], 2)]

            # Plot projection points
            plt.imshow(img, cmap="gray")
            plt.plot(
                projection_points[:, 0],
                projection_points[:, 1],
                ".b",
            )

            plt.show()

            # Angular muscle thickness
            right_half = np.arange(0, len(thickness), 2)
            left_half = np.arange(1, len(thickness), 2)

            # Order thickness to go from 0 to 2pi
            ordered_thickness = np.concatenate(
                (thickness[right_half], thickness[left_half])
            )

            x_values = np.linspace(0, 2 * np.pi, len(thickness),
                                   endpoint=False)
            plt.plot(x_values, ordered_thickness)
            plt.show()


excludeCentralPointsTest()
