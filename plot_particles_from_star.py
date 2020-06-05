#!/usr/bin/env python
# Simple script to plot particles coordinate from star file

import os, sys, argparse, shutil, os.path, glob, string
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from shutil import copyfile

def plot_2d
	""" Plotting 2D array """
	# Prepare the data
	x = np.linspace(0, 10, 100)
	# Plot the data
	plt.plot(x, x, label='linear')

	# Add a legend
	plt.legend()

	# Show the plot
	plt.show()
	
if __name__=='__main__':
	# get name of input starfile, output starfile, output stack file
	parser = argparse.ArgumentParser(description='Perform transform and average of microtubule segment and generate a new star file for the average')
	parser.add_argument('--istar', help='Input particle star file',required=True)
	#parser.add_argument('--odir', help='Output dir for transformed and average files',required=True)
	#parser.add_argument('--ostar', help='Output star file for the average segments',required=True)
	#parser.add_argument('--avgclass', help='Average Class or not',required=False,default=0)
	#parser.add_argument('--nomicro', help='Test mode for only this number of micrographs',required=False)

	args = parser.parse_args()
	
	
	
	instar = open(args.istar, 'r')
	starlabels = learnstarheader(instar)
	coorxcol = starcol_exact_label(starlabels, '_rlnCoordinateX')
	coorycol = starcol_exact_label(starlabels, '_rlnCoordinateY')
	firstline = 1

	for line in instar:
		record = line.split()
		if len(record)==len(starlabels): # if line looks valid
			# Extract coordinate
			if firstline == 1:
				origin = np.array([float(record[coorxcol]) float(record[coorycol])])
			else:
				origin = np.append(origin, [float(record[coorxcol]) float(record[coorycol])], 0)
			firstline = 0
		
	# 
	print(origin)
	instar.close()
