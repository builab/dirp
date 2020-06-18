#!/usr/bin/env python3
# Plot MT origin from Frealign par file

import os, sys, argparse, os.path, glob, math
import matplotlib.pyplot as plt
import numpy as np

def plotorigins(pts, shift, filename):
	""" Plotting 2D array """
	# Plot the data
	ptscorr = np.array(pts - shift)
	plt.figure()
	plt.plot(ptscorr[:, 0], ptscorr[:,1], color='lightgreen', linestyle='--', linewidth=3)
	plt.plot(pts[:, 0], pts[:, 1],  'r+', label='origin')
	plt.plot(ptscorr[:, 0], ptscorr[:,1], 'b.', label='shift adjusted origin')
	plt.xlabel('X')
	plt.ylabel('Y')
	plt.legend(loc='upper right', shadow=True, fontsize='x-large')
	plt.axis('equal')
	plt.savefig(filename)
	plt.close()
	
def readparfile(infile):
	"""The whole spider file content into an np array of float"""
	data = np.loadtxt(infile, dtype='float', comments='C')
	return data
	
if __name__=='__main__':
	# get name of input starfile, output starfile, output stack file
	parser = argparse.ArgumentParser(description='Plot coordinate of Frealign file. Need each MT a unique number in column 8')
	parser.add_argument('--i', help='Input par file',required=True)
	parser.add_argument('--origin', help='Origin file in txt two column',required=True)
	parser.add_argument('--odir', help='Output dir for plots',required=True)
	parser.add_argument('--angpix', help='Pixel Size in Angstrom',required=True)

	args = parser.parse_args()
	# Hardcode step
	parfile = args.i
	angpix = float(args.angpix)
	
	outdir = args.odir
	
	# Check if the outstar file exists
	if os.path.exists(args.odir) == 0:
		print ('Output dir does not exists')
		sys.exit(0)
	
	"""The whole par file content into an np array of float"""
	data = readparfile(parfile)
	origins = readparfile(args.origin)
	if len(data) != len(origins):
		print('Par file and origin file do not have the same number of record')
		sys.exit(0)
		
	mtlist = np.unique(data[:,7])
	#for mtid in mtlist:
	for mtid in mtlist[:3]:
		mtorigin = origins[data[:,7] == mtid, :]
		mtshift = data[data[:,7] == mtid, 4:5]
		mtshift = np.array(mtshift)/angpix
		plotorigins(mtorigin, mtshift, args.odir + '/MT_' + "{}".format(int(mtid)) + '.png')
		