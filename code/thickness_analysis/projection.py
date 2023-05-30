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


def findCenterline(img_stack, horn="left"):
	""" Find the centreline for the selected horn
	
	Arguments:
	img_stack -- ndarray, stack of images to process (ZXY dimensions)
	horn -- str {left, right}, horn to process, default value left.

	Return:
	centreline -- tuple(ndarry), tuple containing the coordinates of the
		centreline.

	"""
	if not horn in {"left", "right"}:
		sys.stderr.write("Error: {} invalid horn selection.\n".format(horn))
		exit(1)

	nb_images = img_stack.shape[0]
	centreline = np.zeros([nb_images, 2], dtype=int)

	for i, img in enumerate(img_stack):
		label_img = skme.label(img)
		regions = skme.regionprops(label_img)
		label = findHornRegion(regions, horn)
		centreline[i, :] = np.round(regions[label]['centroid'])

		# Refine the centroid location by filling the horn
		filled_img = skmo.flood(img, tuple(centreline[i, :]))
		regions = skme.regionprops(filled_img.astype(np.uint8))
		label = findHornRegion(regions, horn)
		
		centreline[i, :] = np.round(regions[label]['centroid'])
			
	return centreline	


def findHornRegion(regions, horn="left"):
	""" Finds the label associated with the desired horn

	Arguments:
	regions -- list of RegionProperties, properties of the regions.
	horn -- str {left, right}, horn to process, default value left.

	Return:
	label -- int, label of the region corresponding to the horn.	

	"""
	extrema = None
	label = 0

	# Find correct regions associated with uterine horns
	for region in regions:
		if region['area'] > 300:
			tmp_value = region['centroid'][1]
			
			if extrema == None:
				extrema = tmp_value

			if horn == "left" and extrema >= tmp_value:
				extrema = tmp_value
				label = region['label'] - 1

			elif horn == "right" and extrema <= tmp_value:
				extrema = tmp_value						
				label = region['label'] - 1

	return label	


def findProjectionPoints(img, centre_point, nb_points):
	""" Find the projection points from the centre point onto the muscle
	layers in four directions (top, bottom, left, right)

	Arguments:
	img -- ndarray, image to analyse.
	centre_point -- list[int], coordinates of the centre point.
	nb_points -- int, number of desired projection points, must be a 
		multiple of 2.

	Return:
	projection_points -- ndarray, list of the coordinates of the
		projection points.

	"""
	if not nb_points % 2 == 0:
		sys.stderr.write("Error: the number of points needs to be even.\n")
		exit(1)

	# Ensure that coordinates are integers
	x_centre = int(np.round(centre_point[1]))
	y_centre = int(np.round(centre_point[0]))

	angles = np.arange(1, 1+(nb_points/2)) * (np.pi / (nb_points / 2))
	projection_points = np.zeros((nb_points*2, 2), dtype=int)

	for i, theta in enumerate(angles):
		if theta != np.pi:
			y_points = np.arange(img.shape[0], dtype=int)
			x_points = (((y_centre - y_points) * np.cos(theta)) /
				 np.sin(theta) + x_centre)
			x_points = x_points.astype(int) # Convert to int
			
			# Find the points that are in the image
			intersection = np.intersect1d(np.where(x_points >= 0),
				np.where(x_points < img.shape[1]))
			x_points = x_points[intersection]
			y_points = y_points[intersection]

		else:
			x_points = np.arange(img.shape[1], dtype=int)
			y_points = np.ones(img.shape[1], dtype=int) * y_centre

		line_y, line_x = skd.line(y_points[0], x_points[0], y_points[-1],
			x_points[-1])
		line = img[line_y, line_x]
		coords = np.where(line[:-1] != line[1:])[0]

		for j in np.arange(0, len(coords), 2):
			coords[j] += 1

		projection_points[i*4:(i+1)*4] = createProjectionPointCoords(
			line_x[coords], line_y[coords], centre_point)

	return projection_points


def createProjectionPointCoords(x_coords, y_coords, centre_point):
	""" Creates the (x, y) pairs of coordinates for the projection points

	Arguments:
	x_coords -- list[int], list of coordinates of the projection points on
		the x axis.
	y_coords -- list[int], list of coordinates of the projection points on
		the y axis.
	centre_point -- ndarray, coordinates of the centre point.
	
	Return:
	projection_points -- ndarray, list of coordinates of the 
		projection points.

	"""
	if len(x_coords) != len(y_coords):
		sys.stderr.write("Error: x_coords and y_coords should have" \
			" the same size.\n")
		exit(1)

	point_list = np.transpose([y_coords, x_coords])

	if len(x_coords) == 4:
		return point_list

	diff = point_list - centre_point
	points_below = np.where(diff[:, 1] < 0)[0] # Before x centre
	points_above = np.where(diff[:, 1] > 0)[0] # After x centre

	if points_below.size == 0:
		points_below = np.where(diff[:, 0] > 0)[0]
		points_above = np.where(diff[:, 0] < 0)[0]

	lower_idx = points_below[
		np.argmin(np.linalg.norm(diff[points_below], axis=1))]
	upper_idx = points_above[
		np.argmin(np.linalg.norm(diff[points_above], axis=1))]

	if diff[upper_idx, 0] < 0:
		upper_points = point_list[upper_idx-1:upper_idx+1]
		lower_points = point_list[lower_idx:lower_idx+2]
		projection_points = np.concatenate((
			upper_points, lower_points))	

	elif diff[upper_idx, 0] >= 0:
		upper_points = point_list[upper_idx:upper_idx+2]
		lower_points = point_list[lower_idx-1:lower_idx+1]
		projection_points = np.concatenate((
			lower_points, upper_points))
	

	return projection_points


def estimateMuscleThickness(img_stack, centreline, nb_points, slice_nbs):
	""" Estimates the muscle thickness of each slice
	
	Arguments:
	img_stack -- ndarray, stack of images to process (ZXY dimensions).
	centreline -- ndarray, coordinates of the centre point on each slice.
	nb_points -- int, number of points to use for projection.
	slice_nbs -- list[int], number of the slices at which to save
		angular thickness.

	Return:
	muscle_thickness_array -- ndarray, muscle thickness of each slice.
	slice_thickness_array -- ndarray, muscle thickness at different angles
		for three slices.
	
	"""
	nb_imgs = len(img_stack)
	muscle_thickness_array = np.zeros(nb_imgs)
	slice_thickness_array = list()

	for i, img in enumerate(img_stack):
		projection_points = findProjectionPoints(img, centreline[i, :], 
			nb_points)

		diff = np.diff(projection_points, axis=0)
		norm = np.linalg.norm(diff, axis=1)	
		thickness = norm[np.arange(0, projection_points.shape[0], 2)]
		muscle_thickness_array[i] = np.mean(thickness)
		
		if i in slice_nbs:
			# Order thickness to go from 0 to 2pi
			ordered_thickness = np.append(
				thickness[np.arange(0, len(thickness), 2)],
				thickness[np.arange(1, len(thickness), 2)])

			# Roll array to line up 0 with anti-mesometrial border
			max_idx = np.argmax(ordered_thickness)
			ordered_thickness = np.roll(ordered_thickness, nb_points - max_idx)
			slice_thickness_array.append(ordered_thickness)

	return muscle_thickness_array, np.transpose(slice_thickness_array)
