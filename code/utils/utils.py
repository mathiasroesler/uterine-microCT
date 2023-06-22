#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# utils.py: Utility functions for the quantative analysis
# Author: Mathias Roesler
# Last modified: 02/23

import os
import sys
import glob
import tomli
import numpy as np
import skimage.io as skio

HOME = os.path.expanduser('~')
BASE = "Documents/phd"


def loadImageStack(dir_path, extension="png"):
	""" Loads the images found in the directory
	
	Arguments:
	dir_path -- str, path to the folder containing the images.	
	extension -- str, image extension, default value png.

	Return:
	img_stack -- ndarray, image stack.
	
	"""
	if not os.path.isdir(dir_path):
		sys.stderr.write("Error: the directory {} does not exists.\n".format(
			dir_path))
		exit()	

	img_list = sorted(glob.glob("*.{}".format(extension), root_dir=dir_path))

	# Load first image to get size of image
	img = skio.imread(os.path.join(dir_path, img_list[0]), as_gray=True)
	nb_x_pixels, nb_y_pixels = img.shape
	stack_size = len(img_list)

	# Pre-allocate and read images
	img_stack = np.zeros([stack_size, nb_x_pixels, nb_y_pixels], dtype=np.uint8)

	for i, img_name in enumerate(img_list):
		path = os.path.join(dir_path, img_name)

		if not os.path.isfile(path):
			sys.stderr.write("Error: {} is not a file.\n".format(path))
			exit()

		img_stack[i, :, :] = skio.imread(path, as_gray=True)
	
	return img_stack


def saveImageStack(img_stack, save_path, img_prefix, start_nb=0, 
		extension="png"):
	""" Saves the images in the stack to the save directory

	Arguments:
	img_stack -- ndarray, stack of images to save.
	save_path -- str, path of the directory in which to save images.
	img_prefix -- str, prefix for the images in the stack.
	start_nb -- int, number at which to start saving the images, 
		default value 0.
	extension -- str, image extension, default value png.

	Return:
	
	"""
	if not os.path.isdir(save_path):
		sys.stderr.write("Error: the directory {} does not exists.\n".format(
			save_path))
		exit()	
	
	i = 0

	for img in img_stack:
		if not np.sum(img) == 0:
			img_path = os.path.join(save_path, img_prefix + str(
				"%03d" % (i + start_nb)))

			if not img.dtype == np.dtype(np.uint8):
				img = img.astype(np.uint8)

			skio.imsave("{}.{}".format(img_path, extension), img, 
				check_contrast=False)
			
			i += 1

def parseTOML(toml_file):
	""" Parse a toml file
	
	Arguments:
	toml_file -- str, path to the toml file.
	
	Return:
	data -- dict, dictonary containing the parameters.

	"""	
	with open(toml_file, 'rb') as f:
		data = tomli.load(f)
				
	return data	


def movingAverage(array, window_size):
	""" Computes the moving average of an array
		
	Arguments:
	array -- ndarray, array over which to compute the average.
	window_size -- int, size of the window over which to compute.

	Return:
	averaged_array -- ndarray, values of the moving average.

	"""
	array_size = len(array)
	averaged_array = np.zeros(array.shape)
	half_window = window_size // 2
	
	if window_size > len(array):
		sys.stderr.write("Error: window size is greater than array size.\n")
		exit(1)

	elif window_size == len(array):
		return np.mean(array)

	for i in range(array_size):
		win_start = max(0, i-half_window)
		win_end = min(array_size, i+half_window+1)
		
		window = array[win_start:win_end]
		averaged_array[i] = np.mean(window, axis=0)

	return averaged_array


def circularAverage(array, window_size):
	""" Computes the average for a circular array

	Arguments:
	array -- ndarray, array over which to compute the standard deviation.
	window_size -- int, size of the window over which to compute.

	Return:
	mean_array -- ndarray, values of the moving standard deviation.

	"""
	array_size = len(array)
	half_window_size = window_size // 2
	mean_array = np.zeros(array.shape)

	if window_size > len(array):
		sys.stderr.write("Error: window size is greater than array size.\n")
		exit(1)

	elif window_size == len(array):
		return [np.mean(array)]

	for i in range(len(array)):
		win_start = i - half_window_size
		win_end = i + half_window_size

		if abs(win_start) + abs(win_end) < window_size:
			# If window size is odd add 1
			win_end += 1	
		
		# Calculate mean with wrapped edges
		mean_array[i] = np.mean(array[
			np.arange(win_start, win_end) % array_size], axis=0)

	return mean_array


def movingStd(array, window_size):
	""" Computes the standard deviation for each window
		
	Arguments:
	array -- ndarray, array over which to compute the standard deviation.
	window_size -- int, size of the window over which to compute.

	Return:
	std_array -- ndarray, values of the moving standard deviation.

	"""
	array_size = len(array)
	std_array = np.zeros(array_size)
	
	if window_size > len(array):
		sys.stderr.write("Error: window size is greater than array size.\n")
		exit(1)

	elif window_size == len(array):
		return [np.std(array)]

	for i in range(0, array_size // window_size):
		window = array[i*window_size:(i+1)*window_size]
		std_array[i*window_size + window_size // 5] = np.std(window)

	return std_array
	

def writeExElem(file_path, elements):
	""" Writes out the data from a mesh to a exnode file

	Arguments:
	file_path -- str, path to the file to save to.
	elements -- ndarray, list of nodes associated with each tetrahedra,
		size = Nx4

	Return:

	"""
	try:
		assert(elements.shape[1] == 4)

	except AssertionError:
		sys.stderr.write("Error: elements should contain 4 nodes\n")

	with open(file_path, "w") as f:
		# Write the exnode file header
		f.write("Group name: uterus\n")
		f.write("Shape.  Dimension=1\n")
		f.write("#Scale factor sets=1\n")
		f.write("l.lagrange, #Scale factors=2\n")
		f.write("#Nodes=4\n")
		f.write("#Fields=1\n")
		f.write("1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
		f.write(" x. l.Lagrange, no modify, standard node based.\n")
		f.write("  #Nodes=4\n")
		f.write("  1. #Values=1\n")
		f.write("	Value indices: 1\n")
		f.write("	Scale factor indices: 1\n")
		f.write("  2. #Values=1\n")
		f.write("	Value indices: 1\n")
		f.write("	Scale factor indices: 2\n")
		f.write("  3. #Values=1\n")
		f.write("	Value indices: 1\n")
		f.write("	Scale factor indices: 3\n")
		f.write("  4. #Values=1\n")
		f.write("	Value indices: 1\n")
		f.write("	Scale factor indices: 4\n")
		f.write(" y. l.Lagrange, no modify, standard node based.\n")
		f.write("  #Nodes=4\n")
		f.write("  1. #Values=1\n")
		f.write("	Value indices: 1\n")
		f.write("	Scale factor indices: 1\n")
		f.write("  2. #Values=1\n")
		f.write("	Value indices: 1\n")
		f.write("	Scale factor indices: 2\n")
		f.write("  3. #Values=1\n")
		f.write("	Value indices: 1\n")
		f.write("	Scale factor indices: 3\n")
		f.write("  4. #Values=1\n")
		f.write("	Value indices: 1\n")
		f.write("	Scale factor indices: 4\n")
		f.write(" z. l.Lagrange, no modify, standard node based.\n")
		f.write("  #Nodes=4\n")
		f.write("  1. #Values=1\n")
		f.write("	Value indices: 1\n")
		f.write("	Scale factor indices: 1\n")
		f.write("  2. #Values=1\n")
		f.write("	Value indices: 1\n")
		f.write("	Scale factor indices: 2\n")
		f.write("  3. #Values=1\n")
		f.write("	Value indices: 1\n")
		f.write("	Scale factor indices: 3\n")
		f.write("  4. #Values=1\n")
		f.write("	Value indices: 1\n")
		f.write("	Scale factor indices: 4\n")

		for i, nodes in enumerate(elements):
			f.write("Element: {} 0 0\n".format(i+1))
			f.write(" Nodes: \n")
			f.write("  {} {} {} {}\n".format(
				nodes[0], nodes[1], nodes[2], nodes[3])) 
			f.write(" Scale factor: \n")
			f.write("  1 1 1 1\n")


def writeExNode(file_path, nodes, thickness):
	""" Writes out the nodes and the thickness from a mesh 
		to a exnode file

	Arguments:
	file_path -- str, path to the file to save to.
	nodes -- ndarray, list of coordinates for each node.
		size = Nx3
	thickness -- ndarray, list of thickness value for each node.
		size = Nx1

	Return:

	"""
	try:
		# Check for number of coordinates
		assert(nodes.shape[1] == 3)

	except AssertionError:
		sys.stderr.write("Error: nodes should have three coordinates\n")
		exit()

	try:
		# Check that thickness and nodes have the same dimension
		assert(nodes.shape[0] == thickness.shape[0])

	except AssertionError:
		sys.stderr.write("Error: nodes and thickness should have the same " \
		"number of elements\n")
		exit()

	with open(file_path, "w") as f:
		# Write exnode file header
		f.write("Group name: uterus\n")
		f.write("#Fields=2\n")
		f.write("1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
		f.write(" x. Value index=1, #Derivatives=0\n")
		f.write(" y. Value index=2, #Derivatives=0\n")
		f.write(" z. Value index=3, #Derivatives=0\n")
		f.write("2) thickness, field, rectangular cartesian, #Components=1\n")
		f.write(" t. Value index=4, #Derivatives=0\n") 

		for i in range(len(nodes)):
			f.write("Node: {}\n".format(i+1))
			f.write(" {} {} {}\n".format(
				nodes[i][0], nodes[i][1], nodes[i][2]))
			f.write(" {}\n".format(thickness[i, 0]))
