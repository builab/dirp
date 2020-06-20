#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Script to read align star files of filament, recenter and fit new curve 
# using RANSAC and interpolate new origins
# Recenter newCoordinateX = CoordinateX + OriginX*binFactor
# v0.3 Take care of case where constant line & switch X, Y if spread of X is too small
# Take care of segID unique in spider
"""
Created on Sat Jun  6 17:35:42 2020

@author: kbui2
"""

import os, sys, argparse, os.path, glob, math, string
import matplotlib.pyplot as plt
import numpy as np



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
	parser.add_argument('--istar', help='Input particle star file',required=True)
	parser.add_argument('--o', help='Output spider',required=True)
	#parser.add_argument('--nomicro', help='Test mode for only this number of micrographs',required=False)

	args = parser.parse_args()
	
	
	instar = open(args.istar, 'r')
	outspider = open(args.o, 'w')
		
	
	starlabels = learnstarheader(instar)
	coorxcol = starcol_exact_label(starlabels, '_rlnCoordinateX')
	coorycol = starcol_exact_label(starlabels, '_rlnCoordinateY')
	microcol = starcol_exact_label(starlabels, '_rlnMicrographName')
	helicalidcol = starcol_exact_label(starlabels, '_rlnHelicalTubeID')
	psipriorcol = starcol_exact_label(starlabels, '_rlnAnglePsiPrior')
	dfucol = starcol_exact_label(starlabels, '_rlnDefocusU')
	dfvcol = starcol_exact_label(starlabels, '_rlnDefocusV')
	dfacol = starcol_exact_label(starlabels, '_rlnDefocusAngle')
	ctfmeritcol = starcol_exact_label(starlabels, '_rlnCtfFigureOfMerit')

	# Stupid
	tmp = '1 2 3 4 5 6 7 8 9 10 11 12'
	line2 = tmp.split()
	count=1
	helicalid = 0
	prevhelicalid = 0
	prevmicroname = ''
	for line in instar:
		record = line.split()
		if len(record)==len(starlabels): # if line looks valid
			microname=record[microcol]
			micronum = os.path.basename(microname)
			micronum = str.replace(micronum, '.mrc', '')
			micronum = str.replace(micronum, 'tt_wt_dmt_', '')
			micronum = int(micronum)
			#print(record[psipriorcol])
			if microname != prevmicroname:
				prevmicroname = microname
				helicalid += 1
				prevhelicalid = record[helicalidcol]	
			if prevhelicalid != record[helicalidcol]:
				helicalid += 1
				prevhelicalid = record[helicalidcol]
				
			
			# Extract coordinate
			line2[1] = str(10)
			line2[0] = str(count)
			line2[2] = str(micronum)
			line2[3] = helicalid
			line2[4] = record[coorxcol]
			line2[5] = record[coorycol]
			line2[6] = "{:.6f}".format(-float(record[psipriorcol]))
			line2[7] = record[dfucol]
			line2[8] = record[dfvcol]
			line2[9] = record[dfacol]
			line2[10] = "{:.6f}".format(1)
			line2[11] = record[ctfmeritcol]
			writestarline(outspider,line2)
				
	instar.close()
	outspider.close()
