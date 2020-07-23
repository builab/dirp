#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Script to read spider files converted from star to ihrsr format
# It will read the coordinate and estimate the rot angle of IHRSR format


import os, sys, argparse, math
import numpy as np



def readspiderfile(infile):
	"""The whole spider file content into an np array of float"""
	data = np.loadtxt(infile, dtype='float', comments=';')
	return data

def writespiderline(outfile, records):
	"""Write a record (line) to an already open starfile"""
	FORMAT = "%7d %7d %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n"
	outfile.write(FORMAT%(records[0],records[1],records[2],records[3],records[4],records[5],records[6],records[7],records[8],records[9],records[10]))
		
	
if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Convert spider alignment to frealign and star')
	parser.add_argument('--ispider', help='Input spider file',required=True)
	parser.add_argument('--ospider', help='Output IHRSR spider file',required=True)
	#parser.add_argument('--nomicro', help='Test mode for only this number of micrographs',required=False)

	args = parser.parse_args()
	
	
	data = readspiderfile(args.ispider)
	outspider = open(args.ospider, 'w')		

	
	# Reading star file header
	prevseg = 0
	prevmicro = 0	
	for npart in len(data):
		if npart == len(data) - 1 :
			data[npart,6] = data[npart-1,6]
			continue
		micro = data[npart, 2]
		seg = data[npart, 3]
		if micro != prevmicro or seg != prevseg:
			rot=-math.atan2(data[npart, 4]-data[npart + 1, 4],data[npart, 5]-data[npart + 1, 5])
			rot = rot*180/math.pi
			if rot <0 :
				rot=rot+360
			if rot > 360:
				rot=rot-360
			data[npart, 6] = rot
			prevmicro = micro
			prevseg = seg
		else:
			data[npart, 6] = data[npart -1, 6]
	writespiderline(outspider,recrod)

							
	outspider.close()
