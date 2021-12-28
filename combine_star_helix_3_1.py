#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  26 22:56:14 2021
combine_star_helix to avoid the same helicalID

"""


import numpy as np
import starfile
import argparse

def renumber_helicalID(df_part):
	count=1
	helicalid = 0
	prevhelicalid = 0
	prevmicroname = ''
	df_renumber = df_part.copy()
	for i in range(len(df_part)):
		microname=df_part.loc[i, 'rlnMicrographName']
		micronum = os.path.basename(microname)
		# This is highly dependent on the micrograph name, need to fix it in the future
		# Do a temporary fix for now. Asumming file name is xxxx_01234_1-3.mrc
		#micronum = re.sub('.*_(\d\d\d\d\d).mrc', '\\1', micronum)
		#print(micronum)
		#micronum = int(micronum)
		#print(record[psipriorcol])
		if microname != prevmicroname:
			prevmicroname = microname
			helicalid += 1
			prevhelicalid = df_part.loc[i, 'rlnHelicalTubeID']	
		if prevhelicalid != df_part.loc[i, 'rlnHelicalTubeID']:
			helicalid += 1
			prevhelicalid = df_part.loc[i, 'rlnHelicalTubeID']
		df_renumber.loc[i, 'rlnHelicalTubeID'] = helicalid
		
	return df_renumber

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Combine helical star files. Assuming same optic groups')
	parser.add_argument('--i1', help='Input star file 1',required=True)
	#parser.add_argument('--i2', help='Input star file 2',required=True)
	parser.add_argument('--o', help='Output star file',required=True)

	
	args = parser.parse_args()
		

	stardict1 = starfile.read(args.i1)
	
	df_optics1 = stardict1['optics']	
	
	df_particles1 = stardict1['particles']
	
	df_renumber1 = renumber_helicalID(df_particles1)
	
	stardict1['particles'] = df_renumber1
	
	# Offset to load in case many different object. Not use now
	starfile.write(stardict1, args.o)
	
