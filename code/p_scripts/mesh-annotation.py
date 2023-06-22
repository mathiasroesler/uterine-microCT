#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# mesh-annotation.py: Script that adds the thickness data to a vtu mesh
# Author: Mathias Roesler
# Last modified: 06/23

import os
import sys
import meshio
import argparse
import numpy as np
import utils.utils as utils

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=
		"Annotates a .vtu mesh with the thickness values of each slice")

	# Parse input arguments
	parser.add_argument("base_name", type=str, metavar="base-name",
		help="name of the dataset")
	parser.add_argument("--mesh-dir", type=str, default="mesh/",
		help="path from BASE to the mesh, default mesh/")
	parser.add_argument("--data-dir", type=str, default="microCT/data",
		help="path form BASE to the thickness data, default microCT/data")
	parser.add_argument("--horn", type=str, choices={"left", "right", "both"},
		help="horn to annotate, default both", default="both")
	parser.add_argument("--not-d", action='store_true',
		help="flag used if the dataset is not downsampled, default False")

	args = parser.parse_args()

	# Set arguments
	mesh_directory = os.path.join(utils.HOME, utils.BASE, args.mesh_dir)
	thickness_directory = os.path.join(utils.HOME, utils.BASE, args.data_dir, 
		args.base_name)	
	mesh_name = os.path.join(mesh_directory, 
		args.base_name + "_volumetric_mesh")
	
	if not args.not_d:
		# If the dataset is downsampled
		thickness_directory = os.path.join(thickness_directory, "downsampled")

	# Add the muscle segmentation to the load directory
	thickness_directory = os.path.join(thickness_directory, 
		"muscle_segmentation")
	horns = ["left", "right"]
	
	# Read the mesh file
	mesh = meshio.read(mesh_name + ".vtu")
	z_coords = mesh.points[:, 2]
	y_coords = mesh.points[:, 0]
	nb_points = len(z_coords)
	nb_slices = round(z_coords.max())
	
	# Read thickness data
	thickness_data = np.load(thickness_directory + "/muscle_thickness.pkl",
		allow_pickle=True)
	left_thickness = thickness_data['left']
	right_thickness = thickness_data['right']

	# Thickness data dictionary to be added to the mesh 	
	point_data_dict = dict()
	point_data_array = np.zeros((nb_points, 1))
	point_data_name = "thickness"

	for i in range(nb_slices):
		# Get the indices of the points on the slice
		slice_idx_list = np.where((z_coords >= i) * (z_coords < (i+1)))[0]
		
		# Deal with edge cases
		if i == 0:
			slice_idx_list = np.where((z_coords >= (i-1)) *
				(z_coords < (i+1)))[0]
		
		elif i == nb_slices-1:
			slice_idx_list = np.where((z_coords >= i) *
				(z_coords < (i+2)))[0]
				
		# Create the thickness array 
		thickness_array = np.zeros((2, 1))

		if i < len(left_thickness):
			thickness_array[0] = left_thickness[i]

		if i < len(right_thickness):
			thickness_array[1] = right_thickness[i]
		
		# Remove 0 values if the horns are not the same length
		thickness_array = np.delete(thickness_array, 
			np.where(thickness_array == 0)[0])

		# If two horns
		if thickness_array.shape[0] == 2:
			slice_y_points = y_coords[slice_idx_list]

			# Find the middle of the slice
			slice_y_max = slice_y_points.max()
			slice_y_min = slice_y_points.min()
			mid_point = (slice_y_max + slice_y_min) / 2

			left_idx_list = np.where((slice_y_points <= mid_point))[0]
			right_idx_list = np.where((slice_y_points > mid_point))[0]

			point_data_array[slice_idx_list[left_idx_list]] = thickness_array[0]
			point_data_array[slice_idx_list[right_idx_list]] = thickness_array[1]

		# If only one horn left
		elif thickness_array.shape[0] == 1:
			point_data_array[slice_idx_list] = thickness_array[0]

	# Add the data dictionary to the mesh
	point_data_dict[point_data_name] = point_data_array
	mesh.point_data = point_data_dict
	
	# Save new mesh
	mesh.write(mesh_name + "_annotated.vtu")
