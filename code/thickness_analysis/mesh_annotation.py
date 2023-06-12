#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# mesh_annotation.py: Script that adds the thickness data to a vtu mesh
# Author: Mathias Roesler
# Last modified: 06/23

import os
import sys
import utils
import meshio
import argparse
import numpy as np

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=
		"Annotates a .vtu mesh with the thickness values of each slice")

	# Parse input arguments
	parser.add_argument("dataset", type=str, metavar="dataset",
		help="name of the dataset that will be used")
	parser.add_argument("--mesh-dir", type=str, default="mesh/",
		help="path to the directory containing the mesh")
	parser.add_argument("--data-dir", type=str, default="microCT/data",
		help="path to the directory containing the thickness")
	parser.add_argument("--horn", type=str, choices={"left", "right", "both"},
		help="horn to annotate", default="both")

	args = parser.parse_args()

	# Set arguments
	mesh_name = args.dataset + "_volumetric_mesh.vtu"
	mesh_path = os.path.join(utils.HOME, utils.BASE, args.mesh_dir,
		mesh_name)

	if args.horn == "both":
		horns = ["left", "right"]

	else:
		horns = [args.horn]

	# Read the mesh file
	mesh = meshio.read(mesh_path)
	z_coords = mesh.points[:, 2]
	y_coords = mesh.points[:, 0]
	nb_points = len(z_coords)
	nb_slices = round(z_coords.max())
	
	# TODO read thickness data from npy file

	# Thickness data dictionary to be added to the mesh 	
	point_data_dict = dict()
	point_data_array = np.zeros((nb_points, 1))
	point_data_name = "thickness"

	# TODO Extract horn data
	for i in range(nb_slices):
		# Get the indices of the points on the slice
		slice_idx_list = np.where((z_coords >= i) * (z_coords < (i+1)))[0]

		slice_y_points = y_coords[slice_idx_list]

		# TODO Take into account the fact that both horns are not the same size
		# Find the middle of the slice
		slice_y_max = slice_y_points.max()
		slice_y_min = slice_y_points.min()
		mid_point = (slice_y_max + slice_y_min) / 2

		left_idx_list = np.where((slice_y_points <= mid_point))[0]
		right_idx_list = np.where((slice_y_points > mid_point))[0]

		point_data_array[slice_idx_list[left_idx_list]] = 0.05*i
		point_data_array[slice_idx_list[right_idx_list]] = 0.01*i

	# Add the data dictionary to the mesh
	point_data_dict[point_data_name] = point_data_array
	mesh.point_data = point_data_dict
	
	# Save new mesh
	mesh.write("test.vtu")
