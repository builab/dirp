#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  26 22:56:14 2021
using starfile for robust conversion. Very quick
"""

import numpy as np
import starfile
import argparse, os

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Combine helical star files. Assuming same optic groups')
	parser.add_argument('--i', help='Input star file',required=True)
	parser.add_argument('--o', help='Output frealign file',required=True)

	
	args = parser.parse_args()
	
	# Default value
	mag = 10000
	occ = 100
	sigma = 0.5
	
	try:
		os.remove(args.o)
	except OSError as e:
    		print ("File {:s} not exist. Good to go!".format(args.o))
		
	stardict = starfile.read(args.i)
	df_optics = stardict['optics']	
	df_part = stardict['particles']
	
	npart = len(df_part)
	
	out = np.zeros((npart, 16))
	
	out[:, 0] = np.arange(1, npart + 1, 1)
	out[:, 1] = df_part.loc[:, 'rlnAnglePsi'].to_numpy()	
	out[:, 2] = df_part.loc[:, 'rlnAngleTilt'].to_numpy()	
	out[:, 3] = df_part.loc[:, 'rlnAngleRot'].to_numpy()	
	out[:, 4] = df_part.loc[:, 'rlnOriginXAngst'].to_numpy()*-1	
	out[:, 5] = df_part.loc[:, 'rlnOriginYAngst'].to_numpy()*-1	
	out[:, 6] = out[:, 6] + mag
	out[:, 7] = df_part.loc[:, 'rlnHelicalTubeID'].to_numpy()	
	out[:, 8] = df_part.loc[:, 'rlnDefocusU'].to_numpy()	
	out[:, 9] = df_part.loc[:, 'rlnDefocusV'].to_numpy()	
	out[:, 10] = df_part.loc[:, 'rlnDefocusAngle'].to_numpy()	
	out[:, 11] = out[:, 11] + occ
	out[:, 13] = out[:, 13] + sigma
	
	np.savetxt(args.o, out, fmt='%9.2f', delimiter=' ')

