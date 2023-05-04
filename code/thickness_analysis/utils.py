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
MASK_FOLDER = "downsampled/muscle_segmentation"


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

	for horn in ["left", "right"]:
		param_dict = data[horn]
		
		for item in param_dict:
			if param_dict[item] == -1:
				param_dict[item] = None

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
