#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Script to read align star files of filament, recenter and fit new curve 
# using RANSAC and interpolate new origins
# Recenter newCoordinateX = CoordinateX + OriginX*binFactor
# v0.1 Take care of case where constant line & switch X, Y if spread of X is too small
"""
Created on Sat Jun  6 17:35:42 2020

@author: kbui2
"""

import os, sys, argparse, os.path, glob, math
import matplotlib.pyplot as plt
import numpy as np
from sklearn import linear_model
from sklearn.linear_model import RANSACRegressor
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline

class PolynomialRegression(object):
    def __init__(self, degree=3, coeffs=None):
        self.degree = degree
        self.coeffs = coeffs

    def fit(self, X, y):
        self.coeffs = np.polyfit(X.ravel(), y, self.degree)

    def get_params(self, deep=False):
        return {'coeffs': self.coeffs}

    def set_params(self, coeffs=None, random_state=None):
        self.coeffs = coeffs

    def predict(self, X):
        poly_eqn = np.poly1d(self.coeffs)
        y_hat = poly_eqn(X.ravel())
        return y_hat

    def score(self, X, y):
        return mean_squared_error(y, self.predict(X))
	
def calculatepsi(pts):
	"""Calculate the psi from the tangent"""
	x = pts[:,0]
	y = pts[:,1]
	threshold = 3
	# Y axis constant
	q1y = np.quantile(y, 0.25)
	q3y = np.quantile(y, 0.75)
	yq = y[(y>= q1y) & (y <= q3y)]
	if np.std(yq) < threshold:
		psi = np.array(y)*0 - 180
		return psi
	# X axis constant
	q1x = np.quantile(x, 0.25)
	q3x = np.quantile(x, 0.75)
	xq = x[(x>= q1x) & (x <= q3x)]
	if np.std(xq) < threshold:
		psi = np.array(y)*0 - 90
		return psi
	# Normal case
	xd = np.diff(x)
	yd = np.diff(y)
	psirad = np.arctan2(yd,xd)*-1
	psi = np.array(psirad)*180/math.pi
	psi = np.hstack([psi, [psi.max()]])
	return psi
	#print(psi)
	
def ransac_polyfit(x, y, step, filename, threshold = 5, outlierthreshold = 10, switchXYthresh = 20):
	"""RANSAC polyfit"""
	#%%Taking care of y = Constant first
	q1y = np.quantile(y, 0.25)
	q3y = np.quantile(y, 0.75)
	yq = y[(y>= q1y) & (y <= q3y)]
	# threshold = Threshold for spread (pixel)
	#outlierthreshold = 10 # Threshold for eliminate outlier
	# switchXYthresh = 20 if std of X is < 30 pixel, switch XY for robust estimation
	if np.std(yq) < threshold:
		print("WARNING: Y is almost constant")
		xn = np.arange(x.min(), x.max(), step)
		yn = np.array(xn)*0 + np.median(y) 
		neworigin = np.vstack([xn, yn])
		plt.plot(xn, yn, color='lightgreen', linestyle='--', linewidth=3)
		plt.plot(xn, yn, 'ro')
		plt.plot(x, y, 'b+')
		plt.axis('equal')
		plt.savefig(str.replace(filename, '.png', '_WARNING.png'))
		plt.close()
		return neworigin.T
	
	# X = constant
	q1x = np.quantile(x, 0.25)
	q3x = np.quantile(x, 0.75)
	xq = x[(x>= q1x) & (x <= q3x)]
	if np.std(xq) < threshold:
		print("WARNING: X is almost constant")
		yn = np.arange(y.min(), y.max(), step)
		xn = np.array(yn)*0 + np.median(x) 
		neworigin = np.vstack([xn, yn])
		plt.plot(xn, yn, color='lightgreen', linestyle='--', linewidth=3)
		plt.plot(xn, yn, 'ro')
		plt.plot(x, y, 'b+')
		plt.axis('equal')
		plt.savefig(str.replace(filename, '.png', '_WARNING.png'))
		plt.close()
		return neworigin.T
	
	# Normal case
	# Avoid X spread two small, switch X & Y std(X)  < 10
	fliplr = 0
	if x[1] > x[-1]:
		flipud = 1	
	else:
		flipud = 0
	if np.std(xq) < switchXYthresh:
		print('WARNING: Switch XY due to small X spread')
		fliplr = 1
		x2 = y[:, np.newaxis]
		xorig  = x
		x = y
		y = xorig

	else:
		x2 = x[:, np.newaxis]
		
	estimator =  RANSACRegressor()
	model = make_pipeline(PolynomialFeatures(2), estimator)
	model.fit(x2, y)
	#%% Filterting outliner. Doesnt seem necessary
	x_plot = x2[:,0]
	y_plot = model.predict(x2)
	#y_error = np.abs(y - y_plot)
	#x_plot = x[y_error < 	outlierthreshold*np.median(y_error)]
	#y_plot = y_plot[y_error < 	outlierthreshold*np.median(y_error)]

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

	if flipud == 1:
		neworigin = np.fliplr(neworigin)
	if fliplr == 1:
		neworigin = np.flipud(neworigin)
		plt.plot(y_plot, x_plot, color='lightgreen', linestyle='--', linewidth=3)
		plt.plot(yn, xn, 'ro')
		plt.plot(y, x, 'b+')
	else:
		plt.plot(x_plot, y_plot, color='lightgreen', linestyle='--', linewidth=3)
		plt.plot(xn, yn, 'ro')
		plt.plot(x, y, 'b+')
	plt.axis('equal')
	plt.savefig(filename)
	plt.close()
	return neworigin.T

#%%	
 
def ransac_polyfit2(x, y, step, filename):
	ransac = RANSACRegressor(PolynomialRegression(degree=2),
                         residual_threshold=2 * np.std(y),
                         random_state=0)
	ransac.fit(np.expand_dims(x, axis=1), y)
	inlier_mask = ransac.inlier_mask_
	y_hat = ransac.predict(np.expand_dims(x, axis=1))
	plt.plot(x, y, 'bx', label='input samples')
	plt.plot(x[inlier_mask], y[inlier_mask], 'go', label='inliers (2*STD)')
	plt.plot(x, y_hat, 'r-', label='estimated curve')
	neworigin = np.vstack([x, y_hat])
	plt.axis('equal')
	plt.savefig(filename)
	plt.close()
	return neworigin.T

 

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
		
