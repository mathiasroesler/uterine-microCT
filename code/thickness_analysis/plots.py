#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# plots.py: Functions used for plotting
# Author: Mathias Roesler
# Last modified: 03/23

import sys
import sys
import numpy as np
import matplotlib.pyplot as plt



def plotProjectionPoints(img, centre, projection_points):
	""" Plots the centre point and projection points of an image

	Arguments:
	img -- ndarray, image to display.
	centre -- ndarray, coordinates of the centre point.
	projection_points -- ndarray, coordinates of the projection points.

	Return:
	
	"""
	fig, ax = plt.subplots()

	plt.imshow(img, cmap='gray') # Plot image
	ax.plot(centre[1], centre[0], '.r', markersize=18) # Plot centre

	for point in projection_points:
		ax.plot(point[1], point[0], '.b', markersize=18)
	
	plt.show()
	

def plotMuscleThickness(muscle_thickness, errors):
	""" Plots the muscle thickness of both horns on the same plot

	The number of points is normalised so that both sets are shown
	between 0 and 1

	Arguments:
	muscle_thickness -- dict(ndarray), thickness of the horns.
	errors -- dict(ndarray), errors for the muscle thickness.

	Return:
	
	"""
	fig, ax = plt.subplots()
	colors = {"left": "black", "right": "silver"} 

	for horn in muscle_thickness.keys():
		horn_thickness = muscle_thickness[horn]
		error_bars = errors[horn]
		horn_length = len(horn_thickness)

		ax.errorbar(np.linspace(0, 1, horn_length), horn_thickness, 
			yerr=error_bars,
			label="{} horn".format(horn.capitalize()), linewidth=4, 
			color=colors[horn])

	ax.tick_params(length=12, width=4, labelsize=22)
	
	# Reset x-axis ticks
	plt.xticks(ticks=[0, 0.5, 1], 
		labels=["Cervical end", "Centre", "Ovarian end"])

	plt.ylim([0, 0.45])
	plt.xlim([0, 1])
	plt.xlabel("Location", fontsize=22)
	plt.ylabel("Muscle thickness (in mm)", fontsize=22)
	plt.legend(fontsize=22)

	plt.show()


def plotAngularThickness(slice_thickness):
	""" Plots the muscle thickness of one slice as a function of the 
	angle theta

	Arguments:
	slice_thickness -- dict(ndarray), array containing the angluar thickness
		of three slices for each horn.

	Return:
	
	"""
	fig, ax = plt.subplots()
	colors = {"left": "black", "right": "silver"} 

	for horn in slice_thickness.keys():
		# Add last value to top of each slice
		y_values = slice_thickness[horn]

		# Create x-axis values
		nb_points = y_values.shape[0]
		x_values = np.arange(nb_points)

		ax.plot(x_values, y_values[:, 0], 
			linestyle='solid', color=colors[horn],
			label="{} horn cervical end".format(horn.capitalize()), 
			linewidth=4)
		ax.plot(x_values, y_values[:, 1],
			linestyle='dashed', color=colors[horn],
			label="{}  horn centre".format(horn.capitalize()), 
			linewidth=4)
		ax.plot(x_values, y_values[:, 2],
			linestyle='dotted', color=colors[horn],
			label="{} horn ovarian end".format(horn.capitalize()),
			linewidth=4)

	# Change tick parameters
	ax.tick_params(length=12, width=4, labelsize=22)
	plt.xticks(ticks=
		[nb_points // 4, nb_points // 2,  3*nb_points // 4, nb_points-1],
		labels=[r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{4}$" ,r"$2\pi$"])

	# Set labels and legends
	plt.xlabel(r"Angle $\theta$ (in rad)", fontsize=22)
	plt.ylabel("Muscle thickness (in mm)", fontsize=22)
	plt.legend(fontsize=22)

	plt.ylim([0, 0.65])
	plt.xlim([0, nb_points-1])
	plt.show()
