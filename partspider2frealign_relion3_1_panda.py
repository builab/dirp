#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Script to read aligned spider files of filament to star file
# Convert spider into array of float for efficient memory
# This can be improved using Panda to write to star file
# To be fixed in the future


import os, sys, argparse, math
import pandas as pd
import starfile
import numpy as np

def readspiderfile(infile):
	"""The whole spider file content into an np array of float"""
	data = np.loadtxt(infile, dtype='float', comments=';')
	return data

def writefrealignline(outfile, records):
	"""Write a record (line) to an already open starfile"""
	FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	outfile.write(FORMAT%(records[0],records[1],records[2],records[3],records[4],records[5],records[6],records[7],records[8],records[9],records[10],records[11],records[12],records[13],records[14],records[15]))
	
	
if __name__=='__main__':
	# get name of input starfile, output starfile, output stack file
	parser = argparse.ArgumentParser(description='Convert spider alignment to frealign and star')
	parser.add_argument('--ispider', help='Input spider file',required=True)
	parser.add_argument('--istar', help='Input segment star file',required=True)
	parser.add_argument('--angpix', help='Input particle pixel size',required=True)
	parser.add_argument('--ofreali', help='Output frealign file',required=True)
	parser.add_argument('--ostar', help='Corresponding output star file',required=True)
	#parser.add_argument('--nomicro', help='Test mode for only this number of micrographs',required=False)

	args = parser.parse_args()
	
	stardict = starfile.read(args.istar)
	
	# Check optics group
	if stardict['optics'] is not None:
		print("WARNING: Script will fail since star file is not in 3.1 format")
		exit(0)
	
	df_part = stardict['particles']
	pixelsize = float(stardict['optics'].loc[0, 'rlnImagePixelSize'])
	
	# Should be done as panda finally
	outfreali = open(args.ofreali, 'w')
	#outstar = open(args.ostar, 'w')		
	
	data = readspiderfile(args.ispider)
	
	firstline = 1
	
	# Stupid
	tmp = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16'
	linefr = tmp.split()
	mag = 10000
	occ = 100
	sigma = 0.5

	groupnumber = 0
	helicalid = 0
	prevhelicalid = 0
	prevmicroname = ''

	for npart in range(len(df_part)):
		partandstack=df_part.loc[npart, 'rlnImageName'].split('@')
		imagename=partandstack[1]
		basename = os.path.basename(imagename)
		# Take care of segavg or particles
		basename = str.replace(basename, '.mrcs', '');
			
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
				
		# Convert frealign
		linefr[0] = npart + 1
		linefr[1] = psi
		linefr[2] = theta
		linefr[3] = phi
		linefr[4] = shx
		linefr[5] = shy
		linefr[6] = mag
		linefr[7] = helicalid
		linefr[8] = float(df_part.loc[npart, 'rlnDefocusU'])
		linefr[9] = float(df_part.loc[npart, 'rlnDefocusV'])
		linefr[10] = float(df_part.loc[npart, 'rlnDefocusAngle'])
		linefr[11] = occ
		linefr[12] = 0
		linefr[13] = sigma
		linefr[14] = 0
		linefr[15] = 0
			

				
		df_part.loc[npart, 'rlnAngleTiltPrior'] = theta
		df_part.loc[npart, 'rlnAnglePsiPrior'] = psi

			
		
		if len(record) < 22:
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
			
			

		writefrealignline(outfreali,linefr)
		writestarline(outstar, record)
			

				
	outfreali.close()
	outstar.close()
