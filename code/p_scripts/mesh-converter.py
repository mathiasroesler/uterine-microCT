#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# mesh-converter.py: Script that converts a vtu mesh to cmgui format
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
		"Converts an annotated vtu mesh to a cmgui format")

	# Parse input arguments
	parser.add_argument("mesh_name", type=str, metavar="mesh-name",
		help="name of the mesh to convert")
	parser.add_argument("--mesh-dir", type=str, default="mesh/",
		help="path from BASE to the mesh, default mesh/")

	args = parser.parse_args()

	# Set arguments
	mesh_path = os.path.join(utils.HOME, utils.BASE, args.mesh_dir)
	mesh_file = mesh_path + '/' + args.mesh_name
	
	# Read the mesh file
	mesh = meshio.read(mesh_file + ".vtu")
	
	# Extract information
	nodes = mesh.points
	thickness = mesh.point_data["thickness"]
	elements = mesh.cells_dict["tetra"]

	# Write EX files
	print("Writing exnode file")
	utils.writeExNode(mesh_file + ".exnode", nodes, thickness)

	print("Writing exelem file")
	utils.writeExElem(mesh_file + ".exelem", elements)

