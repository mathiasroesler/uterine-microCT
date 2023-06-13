#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# projection.py: Functions used to find the projections points
# Author: Mathias Roesler
# Last modified: 03/23

import sys
import plots
import numpy as np
import skimage.io as skio
import skimage.draw as skd
import skimage.measure as skme
import skimage.transform as skt
import skimage.morphology as skmo


def findLineCoordinates(img_shape, centre_point, theta):
	""" Finds the line in the image given an angle and a point.

	Arguments:
	img_shape -- ndarray, shape of the image that is being analysed.
	centre_point -- ndarray, coordinates of the centre point (XY).
	theta -- angle in radians.
	
	Return:
	line_x -- list[int], x coordinates of the line that belong to the
		image.
	line_y -- list[int], y coordinates of the line that belong to the
		image.

	"""
	# Ensure that coordinates are integers
	x_centre = int(np.round(centre_point[0]))
	y_centre = int(np.round(centre_point[1]))

	if theta != np.pi:
		y_points = np.arange(img_shape[0], dtype=int)
		x_points = (((y_centre - y_points) * np.cos(theta)) /
			 np.sin(theta) + x_centre)
		x_points = x_points.astype(int) # Convert to int
		
		# Find the points that are in the image
		intersection = np.intersect1d(np.where(x_points >= 0),
			np.where(x_points < img_shape[1]))
		x_points = x_points[intersection]
		y_points = y_points[intersection]

	else:
		# Case of the vertical line that goes through the centre
		x_points = np.arange(img_shape[1], dtype=int)
		y_points = np.ones(img_shape[1], dtype=int) * y_centre

	# Create the line and get the pixel values that it cuts through
	line_y, line_x = skd.line(y_points[0], x_points[0], y_points[-1],
		x_points[-1])

	return line_x, line_y


def findProjectionPoints(img, centre_point, nb_points, horn):
	""" Find the projection points from the centre point onto the muscle
	layers given the desired number of points.

	Arguments:
	img -- ndarray, image to analyse.
	centre_point -- list[int], coordinates of the centre point (XY).
	nb_points -- int, number of desired projection points, must be a 
		multiple of 2.
	horn -- str {left, right}, horn that is been analysed

	Return:
	projection_points -- ndarray, list of the coordinates of the
		projection points.

	"""
	try:
		assert(nb_points % 2 == 0)

	except AssertionError:
		sys.stderr.write("Error: the number of points needs to be even.\n")
		exit(1)

	angles = np.arange(1, 1+(nb_points/2)) * (np.pi / (nb_points / 2))
	projection_points = np.zeros((nb_points*2, 2), dtype=int)
	
	if (centre_point[2:4] != np.array([0, 0])).all():
		# The horns are not clearly separated and three points are given
		# Create a vector between the left and right points
		v = centre_point[1:2] - centre_point[4:6]
		v = v / np.linalg.norm(v)

		# Use the middle point to draw a line
		line_x, _ = findLineCoordinates(img.shape, centre_point[2:4], 
			np.arccos(np.dot(v, [0, -1])))

		for i in range(len(line_x)):
			# Clear half of the image based on the horn
			if horn == "left":
				img[i, line_x[i]:] = 0

			elif horn == "right":
				img[i, :line_x[i]] = 0;

		# Select the correct centre point based on the horn
		if horn == "left":
			centre_point = centre_point[0:2]

		elif horn == "right":
			centre_point = centre_point[4:6]
	
	for i, theta in enumerate(angles):
		# Find the line for the given angle
		line_x, line_y = findLineCoordinates(img.shape, centre_point, theta)
		line = img[line_y, line_x]

		# Find the indices of rising and falling edges
		coords = np.where(line[:-1] != line[1:])[0] 
		for j in np.arange(0, len(coords), 2):
			coords[j] += 1

		x_coords = line_x[coords]
		y_coords = line_y[coords]

		projection_points[i*4:(i+1)*4] = createProjectionPointCoords(
			x_coords, y_coords, centre_point, theta, img)

	return projection_points


def createProjectionPointCoords(x_coords, y_coords, centre_point, theta, img):
	""" Creates the (x, y) pairs of coordinates for the projection points

	The points that are to the right of the centre point are placed first
	in the array. In the case of the vertical line, the top points are 
	placed first.

	Arguments:
	x_coords -- list[int], list of coordinates of the projection points on
		the x axis.
	y_coords -- list[int], list of coordinates of the projection points on
		the y axis.
	centre_point -- ndarray, coordinates of the centre point.
	theta -- float, angle at which the projection line is on.
	
	Return:
	projection_points -- ndarray, list of coordinates of the 
		projection points.

	"""
	try:
		assert(len(x_coords) == len(y_coords))

	except AssertionError:
		sys.stderr.write("Error: x_coords and y_coords should have" \
			" the same size.\n")
		exit(1)

	point_list = np.transpose([x_coords, y_coords])
	
	if theta == np.pi:
		# In the case of the horizontal line flip x and y
		point_list = np.transpose([y_coords, x_coords])
		centre_point = np.flip(centre_point)

	# Find the indices of the two nearest points to the centre
	diff = point_list - centre_point
	distances_neg = np.linalg.norm(diff[diff[:, 1] < 0], axis=1)
	distances_pos = np.linalg.norm(diff[diff[:, 1] >= 0], axis=1)
	point_list_neg = point_list[diff[:, 1] < 0]
	point_list_pos = point_list[diff[:, 1] >= 0]

	if len(point_list) == 4:
		# If 4 points were found the indices are obvious
		min_idx = 1
		max_idx = 0
	
	else:
		nearest_points_idx_neg = np.argmin(distances_neg)
		nearest_points_idx_pos = np.argmin(distances_pos)
		min_idx = np.min(nearest_points_idx_neg)
		max_idx = np.max(nearest_points_idx_pos)
		
		if min_idx == 0:
			# If only one point is found
			min_idx = 1


	if theta == np.pi:
		# Flip point list back to XY
		point_list_neg = np.fliplr(point_list_neg)
		point_list_pos = np.fliplr(point_list_pos)

	# Create the sets of points on the inner and outer edges
	first_set = point_list_neg[min_idx-1:min_idx+1]
	second_set = point_list_pos[max_idx:max_idx+2]

	if diff[min_idx, 1] < 0:
		projection_points = np.concatenate((first_set, second_set))

	elif diff[min_idx, 1] > 0:
		projection_points = np.concatenate((second_set, first_set))

	return projection_points


def estimateMuscleThickness(img_stack, centreline, nb_points, slice_nbs, horn):
	""" Estimates the muscle thickness of each slice
	
	Arguments:
	img_stack -- ndarray, stack of images to process (ZXY dimensions).
	centreline -- ndarray, coordinates of the centre point on each slice.
	nb_points -- int, number of points to use for projection.
	slice_nbs -- list[int], number of the slices at which to save
		angular thickness.
	horn -- str {left, right}, horn that is been analysed

	Return:
	muscle_thickness_array -- ndarray, muscle thickness of each slice.
	slice_thickness_array -- ndarray, muscle thickness at different angles
		for three slices.
	
	"""
	try:
		assert(len(img_stack) == len(centreline))

	except AssertionError:
		sys.stderr.write("Error: the image stack and the centreline do not " \
			"have the same size.\n")
		exit()

	nb_imgs = len(img_stack)
	muscle_thickness_array = np.zeros(nb_imgs)
	slice_thickness_array = list()
	idx_removed_slices = list()

	for i, img in enumerate(img_stack):
		try:
			projection_points = findProjectionPoints(img, centreline[i, :], 
				nb_points, horn)
			diff = np.diff(projection_points, axis=0)
			norm = np.linalg.norm(diff, axis=1)	
			thickness = norm[np.arange(0, projection_points.shape[0], 2)]
			muscle_thickness_array[i] = np.mean(thickness)

			if i in slice_nbs:
				# Find the four quadrants
				quad_1 = thickness[np.arange(0, nb_points // 2, 2)]
				quad_2 = thickness[np.arange(nb_points // 2, nb_points-2, 2)]
				quad_3 = thickness[np.arange(1, nb_points // 2, 2)]
				quad_4 = thickness[np.arange(1+nb_points // 2, nb_points-1, 2)]

				# Add two last points to correct quadrants
				quad_1 = np.concatenate(([thickness[nb_points-1]], quad_1))
				quad_2 = np.concatenate((quad_2, [thickness[nb_points-2]]))

				# Order thickness to go from 0 to 2pi
				ordered_thickness = np.concatenate((
					quad_1, quad_2, quad_3, quad_4))

				# Roll array to line up 0 with anti-mesometrial border
				max_idx = np.argmax(ordered_thickness)
				ordered_thickness = np.roll(ordered_thickness, 
					nb_points - max_idx)
				slice_thickness_array.append(ordered_thickness)

		except:
			# Write to stderr and add slice number to the remove list
			sys.stderr.write("Warning: unable to process image number {}\n".format(
				i))
			idx_removed_slices.append(i)

	# Remove the slices the values of the slices that were not processed
	muscle_thickness_array = np.delete(muscle_thickness_array, 
		idx_removed_slices)

	return muscle_thickness_array, np.transpose(slice_thickness_array)
