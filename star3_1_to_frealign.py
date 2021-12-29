#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  26 22:56:14 2021
using starfile for robust conversion
"""

import numpy as np
import starfile
import argparse, os

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Combine helical star files. Assuming same optic groups')
	parser.add_argument('--i1', help='Input star file 1',required=True)
	parser.add_argument('--i2', help='Input star file 2',required=True)
	parser.add_argument('--o', help='Output star file',required=True)

	
	args = parser.parse_args()
	
	try:
		os.remove(args.o)
	except OSError as e:
    		print ("File {:s} not exist. Good to go!".format(args.o))
		
	stardict1 = starfile.read(args.i1)
	df_optics1 = stardict1['optics']	
	df_particles1 = stardict1['particles']
