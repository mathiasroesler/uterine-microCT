#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# plots.py: Functions used for plotting
# Author: Mathias Roesler
# Last modified: 06/23

import sys
import numpy as np
import matplotlib.pyplot as plt


def plotProjectionPoints(img, centre, projection_points):
	""" Plots the centre point and projection points of an image

	Arguments:
	img -- ndarray, image to display.
	centre -- ndarray, coordinates of the centre point (XY).
	projection_points -- ndarray, coordinates of the projection points (XY).

	Return:
	
	"""
	fig, ax = plt.subplots()

	plt.imshow(img, cmap='gray') # Plot image
	ax.plot(centre[0], centre[1], '.r', markersize=18) # Plot centre

	for point in projection_points:
		ax.plot(point[0], point[1], '.b', markersize=18)
	
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

	ax.tick_params(length=12, width=4, labelsize=24)
	
	# Reset x-axis ticks
	plt.xticks(ticks=[0, 0.2, 0.6, 1], 
		labels=["Cervix", "Cervical end", "Centre", "Ovarian end"])

	plt.ylim([0, 0.7])
	plt.xlim([0, 1])
	plt.xlabel("Locations", fontsize=24)
	plt.ylabel("Muscle thickness (in mm)", fontsize=24)
	plt.legend(fontsize=24)

	plt.show()


def plotAngularThickness(slice_thickness, projection=False):
	""" Plots the muscle thickness of one slice as a function of the 
	angle theta

	Arguments:
	slice_thickness -- dict(ndarray), array containing the angluar thickness
		of four slices for each horn.
	projection -- str, projection type of the plot, default value False.

	Return:
	
	"""
	fig, ax = plt.subplots(len(slice_thickness.keys()), 1, 
		subplot_kw={"polar": projection})
	colors = {"left": "black", "right": "silver"} 

	if not hasattr(ax, "__len__"):
		# If only one subplot is created
		ax = [ax] # Convert to list for rest of code to work
	
	for i, horn in enumerate(slice_thickness.keys()):
		y_values = slice_thickness[horn]

		# Create x-axis values so that everything is normalised
		nb_points = y_values.shape[0]
		x_values = np.linspace(0, 2*np.pi, nb_points, endpoint=False)

		ax[i].plot(x_values, y_values[:, 0], 
			linestyle='dashdot', color=colors[horn],
			label="Cervix", linewidth=4)
		ax[i].plot(x_values, y_values[:, 1], 
			linestyle='solid', color=colors[horn],
			label="Cervical end", linewidth=4)
		ax[i].plot(x_values, y_values[:, 2],
			linestyle='dashed', color=colors[horn],
			label="Centre",	linewidth=4)
		ax[i].plot(x_values, y_values[:, 3],
			linestyle='dotted', color=colors[horn],
			label="Ovarian end", linewidth=4)
		ax[i].set_title("{} horn muscle thickness (in mm)".format(
			horn.capitalize()), fontsize=24)

		# Change tick parameters
		ax[i].tick_params(length=12, width=4, labelsize=24)
		
		if projection:
			ax[i].set_rlabel_position(-22.5)  # Move radial labels

			ax[i].set_rmax(1.1) # Set radial max
			ticks = plt.xticks()[0]

			# Set labels and legends
			angle = np.deg2rad(35)
			plt.legend(loc="lower left", fontsize=24,
				bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))

			plt.xticks(ticks=ticks, labels=['0',r'$\frac{\pi}{4}$',\
					r'$\frac{\pi}{2}$',r'$\frac{3\pi}{4}$', r'$\pi$',\
					r'$\frac{5\pi}{4}$',r'$\frac{3\pi}{2}$',\
					r'$\frac{7\pi}{4}$'])

		else:
			plt.xlim([0, 2*np.pi])
			plt.ylim([0, 1.1])
			ticks = np.linspace(0, 2*np.pi, 9)
			plt.legend(fontsize=24)

			plt.xticks(ticks=ticks, labels=['0',r'$\frac{\pi}{4}$',\
					r'$\frac{\pi}{2}$',r'$\frac{3\pi}{4}$', r'$\pi$',\
					r'$\frac{5\pi}{4}$',r'$\frac{3\pi}{2}$',\
					r'$\frac{7\pi}{4}$', r'2$\pi$'])


	plt.show()
