#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Script to read align star files of filament, recenter and fit new curve 
# using RANSAC and interpolate new origins
# Recenter newCoordinateX = CoordinateX + OriginX*binFactor
# TODO: How to eliminate very bad one already
"""
Created on Sat Jun  6 17:35:42 2020

@author: kbui2
"""

import os, sys, argparse, os.path, glob, math
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import RANSACRegressor
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline

	
def calculatepsi(pts):
	"""Calculate the psi from the tangent"""
	x = pts[:,0]
	y = pts[:,1]
	xd = np.diff(x)
	yd = np.diff(y)
	psirad = np.arctan2(yd,xd)*-1
	psi = np.array(psirad)*180/math.pi
	psi = np.hstack([psi, [psi.max()]])
	return psi
	#print(psi)
	
def ransac_polyfit(x, y, step, filename):
	"""RANSAC polyfit"""
	if x[1] > x[-1]:
		flipud = 1	
	else:
		flipud = 0
	x = x[:, np.newaxis]
	estimator =  RANSACRegressor()
	model = make_pipeline(PolynomialFeatures(2), estimator)
	model.fit(x, y)
	x_plot = np.linspace(x.min(), x.max())
	y_plot = model.predict(x_plot[:, np.newaxis])
	
	# calculate new origin
	xd = np.diff(x_plot)
	yd = np.diff(y_plot)
	dist = np.sqrt(xd**2 + yd**2)
	u = np.cumsum(dist)
	u = np.hstack([[0], u])
	t = np.arange(0, u.max(), step)
	xn = np.interp(t, u, x_plot)
	yn = np.interp(t, u, y_plot)
	neworigin = np.vstack([xn, yn])
	plt.plot(x_plot, y_plot, color='lightgreen', linestyle='--', linewidth=3)
	plt.plot(xn, yn, 'ro')
	plt.plot(x[:, 0], y, 'b+')
	plt.savefig(filename)
	plt.close()
	if flipud == 1:
		neworigin = np.fliplr(neworigin)
	return neworigin.T
	
 
def ransac_polyfit2(x, y, order=1, n=5, k=100, t=0.1, d=5, f=0.8):
	# Thanks https://en.wikipedia.org/wiki/Random_sample_consensus
	
	# n – minimum number of data points required to fit the model
	# k – maximum number of iterations allowed in the algorithm
	# t – threshold value to determine when a data point fits a model
	# d – number of close data points required to assert that a model fits well to data
	# f – fraction of close data points required
	besterr = np.inf
	bestfit = None
	for kk in range(k):
		print(kk)
		maybeinliers = np.random.randint(len(x), size=n)
		maybemodel = np.polyfit(x[maybeinliers], y[maybeinliers], order)
		print(maybemodel)
		alsoinliers = np.abs(np.polyval(maybemodel, x)-y) < t
		if sum(alsoinliers) > d and sum(alsoinliers) > len(x)*f:
			bettermodel = np.polyfit(x[alsoinliers], y[alsoinliers], order)
			thiserr = np.sum(np.abs(np.polyval(bettermodel, x[alsoinliers])-y[alsoinliers]))
			if thiserr < besterr:
				bestfit = bettermodel
				besterr = thiserr
	return bestfit
 

def learnstarheader(infile):
	"""Learn which column contains which information from an already open starfile"""
	infile.seek(0) # Go to the beginning of the starfile
	doneheader = False
	doneprelabels = False
	headerlabels = []
	while not doneprelabels:
		line=infile.readline()
		if line.startswith('loop_'):
			doneprelabels = True # read until 'loop_'
	while not doneheader:
		line=infile.readline()
		if not line.startswith('_'): # read all lines the start with '_'
			doneheader = True
		else:
			headerlabels += [line] 
	infile.seek(0) # return to beginning of starfile before return
	return headerlabels

def writestarheader(outfile,headerlabels):		  
	"""With an already opened starfile write a header"""
	outfile.write('\ndata_\n\nloop_\n')
	for label in headerlabels:
		outfile.write(label)

def readstarline(infile):
	"""Read a record (line) from an already open starfile and return XXX"""
	line=infile.readline()
	records = line.split()
	return records

def writestarline(outfile,records):
	"""Write a record (line) to an already open starfile"""
	for item in records:
		outfile.write(item+'  ')
	outfile.write('\n')

def starcol_exact_label(starlabels, label):
	"""New function to do exact match of relion label such as _rlnImageCol"""
	for i, s in enumerate(starlabels):
		record=s.split()
		if label == record[0]:
			return i
	return -1
	
if __name__=='__main__':
	# get name of input starfile, output starfile, output stack file
	parser = argparse.ArgumentParser(description='Plot coordinate of star file')
	parser.add_argument('--istarpattern', help='Input particle star file',required=True)
	parser.add_argument('--odir', help='Output dir for plots',required=True)
	#parser.add_argument('--ostar', help='Output star file for the average segments',required=True)
	parser.add_argument('--recenter', help='Do recenter',required=False,default=1)
	parser.add_argument('--step', help='Step for interpolation (pixel)',required=True)
	parser.add_argument('--bin', help='Bin factor for current alignment (can be non inetger)',required=True)

	#parser.add_argument('--nomicro', help='Test mode for only this number of micrographs',required=False)

	args = parser.parse_args()
	
	# Hardcode step
	# binfactr = 4
	step = float(args.step)
	recenter = int(args.recenter)
	bin = float(args.bin)
	
	liststar=glob.glob(args.istarpattern)
	outdir = args.odir
	
	# Check if the outstar file exists
	if os.path.exists(args.odir) == 0:
		print ('Output dir does not exists')
		sys.exit(0)
	
	shiftX = 0
	shiftY = 0
	
	
	for starfile in liststar:
		basename = os.path.basename(starfile)
		basename = str.replace(basename, '.star', '')
		print(basename)
		instar = open(starfile, 'r')
		starlabels = learnstarheader(instar)
		coorxcol = starcol_exact_label(starlabels, '_rlnCoordinateX')
		coorycol = starcol_exact_label(starlabels, '_rlnCoordinateY')
		imagecol = starcol_exact_label(starlabels, '_rlnImageName')
		helicaltracklengthcol = starcol_exact_label(starlabels, '_rlnHelicalTrackLength')
		shiftxcol = starcol_exact_label(starlabels, '_rlnOriginX')
		shiftycol = starcol_exact_label(starlabels, '_rlnOriginY')
		psipriorcol = starcol_exact_label(starlabels, '_rlnAnglePsiPrior')
		psicol = starcol_exact_label(starlabels, '_rlnAnglePsi')
		firstline = 1
		for line in instar:
			record = line.split()
			if len(record)==len(starlabels): # if line looks valid
				#print(record[psipriorcol])
				# Extract coordinate
				if recenter == 1:
					shiftX = float(record[shiftxcol])*bin
					shiftY = float(record[shiftycol])*bin
				if firstline == 1:
					origin = np.array([float(record[coorxcol]), float(record[coorycol])])
					singledataline = record
					singleline = record
					origin = np.array([float(record[coorxcol]) + shiftX, float(record[coorycol]) + shiftY])
				else:
					if firstline == 2:
						singletracklength = float(record[helicaltracklengthcol])
					origin = np.vstack((origin, [float(record[coorxcol]) + shiftX, float(record[coorycol]) + shiftY]))
				firstline +=1
				
		
		instar.close()
		#print(origin)
		neworigin = ransac_polyfit(origin[:,0], origin[:,1], step, outdir + "/" + basename + ".png")
		psi = calculatepsi(neworigin)
		#print(psi)
		
		# Writing new star file
		outstar = open(outdir + '/' + basename + '.star', 'w')
		writestarheader(outstar, starlabels)
		medpsi = np.median(psi)
		partandstack=singledataline[imagecol].split('@')
		# Reset shift
		singledataline[shiftxcol] = "{:.6f}".format(0) 
		singledataline[shiftycol] = "{:.6f}".format(0) 
		for npart in range(len(neworigin)):
			# New origin
			singledataline[coorxcol] = "{:.6f}".format(neworigin[npart,0]) 
			singledataline[coorycol] = "{:.6f}".format(neworigin[npart,1]) 
			# New tracklength, image & psi
			singledataline[helicaltracklengthcol] = "{:.6f}".format(singletracklength*npart) 
			singledataline[psicol] = "{:.6f}".format(psi[npart] - medpsi) 
			singledataline[psipriorcol] = "{:.6f}".format(psi[npart]) 
			singledataline[imagecol]=str(npart + 1).zfill(6)+'@'+partandstack[1]
			writestarline(outstar, singledataline)
			#print(singledataline)

		outstar.close()
		
