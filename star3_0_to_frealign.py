#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Convert star file to frealign with helical id column
"""
Created on Sat Jun  6 17:35:42 2020

@author: kbui2
"""

import os, argparse, os.path


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

def writefrealignline(outfile, records):
	"""Write a record (line) to an already open starfile"""
	FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	outfile.write(FORMAT%(records[0],records[1],records[2],records[3],records[4],records[5],records[6],records[7],records[8],records[9],records[10],records[11],records[12],records[13],records[14],records[15]))
	

	
if __name__=='__main__':
	# get name of input starfile, output starfile, output stack file
	parser = argparse.ArgumentParser(description='Plot coordinate of star file')
	parser.add_argument('--istar', help='Input particle star file',required=True)
	parser.add_argument('--opar', help='Output frealign',required=True)
	#parser.add_argument('--nomicro', help='Test mode for only this number of micrographs',required=False)

	args = parser.parse_args()
	
	
	instar = open(args.istar, 'r')
	outfreali= open(args.opar, 'w')
		
	
	starlabels = learnstarheader(instar)
	originxcol = starcol_exact_label(starlabels, '_rlnOriginX')
	originycol = starcol_exact_label(starlabels, '_rlnOriginY')
	microcol = starcol_exact_label(starlabels, '_rlnMicrographName')
	helicalidcol = starcol_exact_label(starlabels, '_rlnHelicalTubeID')
	psicol = starcol_exact_label(starlabels, '_rlnAnglePsi')
	tiltcol = starcol_exact_label(starlabels, '_rlnAngleTilt')
	rotcol = starcol_exact_label(starlabels, '_rlnAngleRot')
	dfucol = starcol_exact_label(starlabels, '_rlnDefocusU')
	dfvcol = starcol_exact_label(starlabels, '_rlnDefocusV')
	dfacol = starcol_exact_label(starlabels, '_rlnDefocusAngle')
	magcol = starcol_exact_label(starlabels, '_rlnMagnification')
	detpixelsizecol = starcol_exact_label(starlabels, '_rlnDetectorPixelSize')

	linefr = [0]*16
	mag = 10000
	occ = "100"
	sigma = 0.5

	npart=1
	helicalid = 0
	microlist ={}
	prevhelicalid = 0
	prevmicroname = ''
	micronum = 1
	for line in instar:
		record = line.split()
		if len(record)==len(starlabels): # if line looks valid
			microname=record[microcol]
			microname = os.path.basename(microname)
			# Create a dictionary, if microname not exist as key, insert
			if microlist.get(microname):
				if prevhelicalid != record[helicalidcol]:
				helicalid += 1
				prevhelicalid = record[helicalidcol]
			else:
				microlist[microname] = micronum
				micronum =+ 1
				helicalid += 1
				prevhelicalid = record[helicalidcol]
			


			# Extract coordinate			
			mag = float(record[magcol])
			detpixelsize = float(record[detpixelsizecol])
			# Convert to Angstrom
			pixelsize = detpixelsize/mag*10000
			shx = -float(record[originxcol])*pixelsize
			shy = -float(record[originycol])*pixelsize
			phi = float(record[rotcol])
			theta = float(record[tiltcol])
			psi = float(record[psicol])
			
			
			
			linefr[0] = npart
			linefr[1] = psi
			linefr[2] = theta
			linefr[3] = phi
			linefr[4] = shx
			linefr[5] = shy
			linefr[6] = int(mag)
			linefr[7] = helicalid
			linefr[8] = record[dfucol]
			linefr[9] = record[dfvcol]
			linefr[10] = record[dfacol]
			linefr[11] = occ
			linefr[12] = 0
			linefr[13] = sigma
			linefr[14] = 0
			linefr[15] = 0
	
			writefrealignline(outfreali,linefr)
			
			npart += 1  
			
			
	instar.close()
	outfreali.close()
