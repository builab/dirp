#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Script to read aligned spider files of filament to star file
# Cannot deal with star file with # in it
# Convert spider into array of float for efficient memory
# This script is not very robust but very efficient


import os, sys, argparse, math
import numpy as np


def learnstarheader(infile):
	"""Learn which column contains which information from an already open starfile"""
	infile.seek(0) # Go to the beginning of the starfile
	doneheader = False
	doneprelabels = False
	headerlabels = []
	headeroptics = []
	doneoptics = True
	while not doneprelabels:
		line=infile.readline()
		if line.startswith('#'):
			continue
		# Check if star 3.1 format
		if line.startswith('data_optics'):
			doneoptics = False
		if line.startswith('data_particles'):
			doneoptics = True
		if line.startswith('loop_') & doneoptics == True:
			doneprelabels = True # read until 'loop_'
		headeroptics += [line]
	while not doneheader:
		line=infile.readline()
		if line.startswith('#'):
			continue
		if not line.startswith('_'): # read all lines the start with '_'
			doneheader = True
		else:
			headerlabels += [line] 
	infile.seek(0) # return to beginning of starfile before return
	return headeroptics, headerlabels

def is_star3_1(infile):
	"""Learn starfile is 3.1 or not"""
	infile.seek(0) # Go to the beginning of the starfile
	is_star3_1 = False
	doneheader = False
	while not doneheader:
		line=infile.readline()
		# Check if star 3.1 format
		if line.startswith('data_optics'):
			doneheader = True
			is_star3_1 = True
		
	infile.seek(0) # return to beginning of starfile before return
	return is_star3_1

def writestarheader(outfile, headeroptics, headerlabels):		  
	"""With an already opened starfile write a header"""
	for label in headeroptics:
		outfile.write(label)
	for label in headerlabels:
		outfile.write(label)

def appendstarlabel(headerlabels, label):
	newheaderlabels = headerlabels + [label + ' \n']
	return newheaderlabels

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

def starcountparticles(starfile):
	"""Count how many particles inside the star file"""
	nopart = 0
	instar = open(starfile, 'r')
	starlabels = learnstarheader(instar)
	for line in instar:
		record = line.split()
		if len(record)==len(starlabels): # if line looks valid
			nopart += 1			
	return nopart

def readspiderfile(infile):
	"""The whole spider file content into an np array of float"""
	data = np.loadtxt(infile, dtype='float', comments=';')
	return data	
	
if __name__=='__main__':
	# get name of input starfile, output starfile, output stack file
	parser = argparse.ArgumentParser(description='Convert spider alignment to frealign and star')
	parser.add_argument('--ispider', help='Input spider file',required=True)
	parser.add_argument('--istar', help='Input segment star file',required=True)
	parser.add_argument('--angpix', help='Input particle pixel size',required=True)
	parser.add_argument('--ostar', help='Corresponding output star file',required=True)
	#parser.add_argument('--nomicro', help='Test mode for only this number of micrographs',required=False)

	args = parser.parse_args()
	

	
	instar = open(args.istar, 'r')
	outstar = open(args.ostar, 'w')		
	
	data = readspiderfile(args.ispider)
	
	if is_star3_1(instar) == False:
		print("WARNING: Script will fail since star file is not in 3.1 format")
		exit(0)
	
	# Reading star file header
	[staroptics, starlabels] = learnstarheader(instar)
	
	#print(staroptics)

	
	if len(starlabels) < 19:
		starlabels = appendstarlabel(starlabels, '_rlnGroupNumber #19')
		starlabels = appendstarlabel(starlabels, '_rlnOriginXAngst #20')
		starlabels = appendstarlabel(starlabels, '_rlnOriginYAngst #21')
		starlabels = appendstarlabel(starlabels, '_rlnAngleRot #22')
		starlabels = appendstarlabel(starlabels, '_rlnAngleTilt #23')
		starlabels = appendstarlabel(starlabels, '_rlnAnglePsi #24')
	
	#print(starlabels)
	
	coorxcol = starcol_exact_label(starlabels, '_rlnCoordinateX')
	coorycol = starcol_exact_label(starlabels, '_rlnCoordinateY')
	orixcol = starcol_exact_label(starlabels, '_rlnOriginXAngst')
	oriycol = starcol_exact_label(starlabels, '_rlnOriginYAngst')
	tiltpriorcol = starcol_exact_label(starlabels, '_rlnAngleTiltPrior')
	psipriorcol = starcol_exact_label(starlabels, '_rlnAnglePsiPrior')
	psicol = starcol_exact_label(starlabels, '_rlnAnglePsi')
	rotcol = starcol_exact_label(starlabels, '_rlnAngleRot')
	tiltcol = starcol_exact_label(starlabels, '_rlnAngleTilt')
	helicalidcol = starcol_exact_label(starlabels, '_rlnHelicalTubeID')
	dfucol = starcol_exact_label(starlabels, '_rlnDefocusU')
	dfvcol = starcol_exact_label(starlabels, '_rlnDefocusV')
	dfacol = starcol_exact_label(starlabels, '_rlnDefocusAngle')
	ctfmeritcol = starcol_exact_label(starlabels, '_rlnCtfFigureOfMerit')
	imagecol = starcol_exact_label(starlabels, '_rlnImageName')
	groupcol = starcol_exact_label(starlabels, '_rlnGroupNumber')
	classcol = starcol_exact_label(starlabels, '_rlnClassNumber')
	microcol = starcol_exact_label(starlabels, '_rlnMicrographName')

	writestarheader(outstar, staroptics, starlabels)
	
	npart=0
	groupnumber = 0
	helicalid = 0
	prevhelicalid = 0
	prevmicroname = ''
	pixelsize = float(args.angpix) #3.1

	for line in instar:
		record = line.split()
		if len(record) > 10: # if line looks valid
			#print(record)
			shx_s = data[npart, 14]
			shy_s = data[npart, 15]
			psi = 360 - data[npart, 13]
			psirad = psi*math.pi/180
			theta = data[npart, 3]
			phi =  data[npart, 4]
			
			# Angular conversion
			shx = -(shx_s*math.cos(psirad) + shy_s*math.sin(psirad))
			shy = -(-shx_s*math.sin(psirad) + shy_s*math.cos(psirad))

			shx = shx*pixelsize
			shy = shy*pixelsize
			
			microname = record[microcol]
			if microname != prevmicroname:
				groupnumber +=1
				prevmicroname = microname
				helicalid += 1
				prevhelicalid = record[helicalidcol]
				
			if prevhelicalid != record[helicalidcol]:
				helicalid += 1
				prevhelicalid = record[helicalidcol]
						
			record[tiltpriorcol] = "{:.6f}".format(theta)
			record[psipriorcol] = "{:.6f}".format(psi)

			
			# Check this for 3_1
			if len(record) < 19:
				record += ["{:5d}".format(groupnumber)]
				record += ["{:.6f}".format(-shx)]
				record += ["{:.6f}".format(-shy)]
				record += ["{:.6f}".format(phi)]
				record += ["{:.6f}".format(psi)]
				record += ["{:.6f}".format(theta)]
			else:
				record[groupcol] = "{:5d}".format(groupnumber)
				record[rotcol] = "{:.6f}".format(phi)
				record[psicol] = "{:.6f}".format(psi)
				record[tiltcol] = "{:.6f}".format(theta)
				record[orixcol] = "{:.6f}".format(-shx)
				record[oriycol] = "{:.6f}".format(-shy)
			
			writestarline(outstar, record)
			npart += 1
			

				
	instar.close()
	outstar.close()
