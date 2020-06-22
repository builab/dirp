#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Script to read aligned frealign par file & original star file
# Convert into new star file
# 2020/06/24 fix convert frealign to relion -shx/apix, -shy/apix


import argparse
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

def readparfile(infile):
	"""The whole spider file content into an np array of float"""
	data = np.loadtxt(infile, dtype='float', comments='C')
	return data

def writefrealignline(outfile, records):
	"""Write a record (line) to an already open starfile"""
	FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	outfile.write(FORMAT%(records[0],records[1],records[2],records[3],records[4],records[5],records[6],records[7],records[8],records[9],records[10],records[11],records[12],records[13],records[14],records[15]))
	
	
if __name__=='__main__':
	# get name of input starfile, output starfile, output stack file
	parser = argparse.ArgumentParser(description='Convert spider alignment to frealign and star')
	parser.add_argument('--ipar', help='Input par file',required=True)
	parser.add_argument('--istar', help='Input segment star file',required=True)
	parser.add_argument('--ostar', help='Corresponding output star file',required=True)

	args = parser.parse_args()
	

	
	instar = open(args.istar, 'r')
	outstar = open(args.ostar, 'w')		

	
	data = readparfile(args.ipar)
	
	# Reading star file header
	starlabels = learnstarheader(instar)

	
	orixcol = starcol_exact_label(starlabels, '_rlnOriginX')
	oriycol = starcol_exact_label(starlabels, '_rlnOriginY')
	tiltpriorcol = starcol_exact_label(starlabels, '_rlnAngleTiltPrior')
	psipriorcol = starcol_exact_label(starlabels, '_rlnAnglePsiPrior')
	psicol = starcol_exact_label(starlabels, '_rlnAnglePsi')
	rotcol = starcol_exact_label(starlabels, '_rlnAngleRot')
	tiltcol = starcol_exact_label(starlabels, '_rlnAngleTilt')
	pixelsizecol = starcol_exact_label(starlabels, '_rlnDetectorPixelSize')

	
	writestarheader(outstar, starlabels)

	npart = 0
	for line in instar:
		record = line.split()
		if len(record)==len(starlabels): # if line looks valid
				# Star file
				angpix = float(record[pixelsizecol])
				record[orixcol] = "{:.6f}".format(-data[npart,4]/angpix)
				record[oriycol] = "{:.6f}".format(-data[npart,5]/angpix)
				record[psipriorcol] = "{:.6f}".format(data[npart,1])
				record[tiltpriorcol] = "{:.6f}".format(data[npart,2])
				record[rotcol] = "{:.6f}".format(data[npart,3])
				record[psicol] = "{:.6f}".format(data[npart,1])
				record[tiltcol] = "{:.6f}".format(data[npart,2])	
				writestarline(outstar,record)
				npart += 1

	outstar.close()
