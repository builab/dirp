#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Script to read aligned spider files of filament to star file
# Convert spider into array of float for efficient memory
# This is extremely slow


import os, sys, argparse, math
import pandas as pd
import starfile
import numpy as np

def readspiderfile(infile):
	"""The whole spider file content into an np array of float"""
	data = np.loadtxt(infile, dtype='float', comments=';')
	return data
	
if __name__=='__main__':
	# get name of input starfile, output starfile, output stack file
	parser = argparse.ArgumentParser(description='Convert spider alignment to frealign and star')
	parser.add_argument('--ispider', help='Input spider file',required=True)
	parser.add_argument('--istar', help='Input segment star file',required=True)
	#parser.add_argument('--angpix', help='Input particle pixel size',required=True)
	parser.add_argument('--ofreali', help='Output frealign file',required=True)
	parser.add_argument('--ostar', help='Corresponding output star file',required=True)
	#parser.add_argument('--nomicro', help='Test mode for only this number of micrographs',required=False)

	args = parser.parse_args()
	
	stardict = starfile.read(args.istar)
	
	# Check optics group
	if stardict['optics'] is None:
		print("WARNING: Script will fail since star file is not in 3.1 format")
		exit(0)
	
	df_part = stardict['particles']
	pixelsize = float(stardict['optics'].loc[0, 'rlnImagePixelSize'])
	#print(df_part)
	
	# Should be done as panda finally
	outfreali = open(args.ofreali, 'w')
	#outstar = open(args.ostar, 'w')		
	
	data = readspiderfile(args.ispider)
	#print(data)
	
	mag = 10000
	occ = 100
	sigma = 0.5
	
	# Extra field
	header_list = ["rlnOriginXAngst", "rlnOriginYAngst", "rlnAngleRot", "rlnAngleTilt", "rlnAnglePsi"]
	df_extra = pd.DataFrame(columns = header_list)	
	df_part_out = df_part.drop(columns=header_list, errors='ignore').copy()
	#df_part_out = df_part.copy()
	#if not 'rlnOriginXAngst' in df_part_out.columns:
		
		
	helicalid = 0
	prevhelicalid = 0
	prevmicroname = ''

	for npart in range(len(df_part_out)):
		# Calculate new helicalID
		microname=df_part_out.loc[npart, 'rlnMicrographName']
		if microname != prevmicroname:
			print(helicalid)
			prevmicroname = microname
			helicalid += 1
			prevhelicalid = df_part_out.loc[npart, 'rlnHelicalTubeID']
			print(helicalid)
		if prevhelicalid != df_part_out.loc[npart, 'rlnHelicalTubeID']:
			helicalid += 1
			prevhelicalid = df_part_out.loc[npart, 'rlnHelicalTubeID']
			print(helicalid)
	
		#partandstack=df_part.loc[npart, 'rlnImageName'].split('@')
		#imagename=partandstack[1]
		#basename = os.path.basename(imagename)
		# Take care of segavg or particles
		#basename = str.replace(basename, '.mrcs', '');
			
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
	
		df_part_out.loc[npart, 'rlnAngleTiltPrior'] = theta
		df_part_out.loc[npart, 'rlnAnglePsiPrior'] = psi
		df_part_out.loc[npart, 'rlnHelicalTubeID'] = helicalid
		
		df_extra.loc[npart, 'rlnOriginXAngst'] = -shx
		df_extra.loc[npart, 'rlnOriginYAngst'] = -shy
		df_extra.loc[npart, 'rlnAngleRot'] = phi
		df_extra.loc[npart, 'rlnAnglePsi'] = psi
		df_extra.loc[npart, 'rlnAngleTilt'] = theta
							

	
	try:
		os.remove(args.ostar)
	except OSError as e:
    		print ("File {:s} not exist. Good to go!".format(args.ostar))
	df_out = pd.concat([df_part_out, df_extra], axis=1)
	stardict['particles'] = df_part_out
	starfile.write(stardict, args.ostar)
		
