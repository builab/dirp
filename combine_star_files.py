#!/usr/bin/env python3
# Script to combine star files pattern into a single star file
# Main code coming from John Rubinstein createstackfromstar.py
# HB 2020/05/30 Tested & verified

import os, sys, argparse, shutil, os.path, glob, string

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
		

if __name__=='__main__':
	# get name of input starfile, output starfile, output stack file
	parser = argparse.ArgumentParser(description='Combine star files into one from a regular expression pattern')
	parser.add_argument('--istarpattern', help='Input particle star file pattern (Extract/MT*.star)',required=True)
	parser.add_argument('--ostar', help='Output combine star file (particles_cbn.star)',required=True)
	parser.add_argument('--nomicro', help='Test mode for only this number of micrograph',required=False)

	args = parser.parse_args()
	#print(args.istar)
	if args.nomicro is not None:
		testmode = 1
		nomicro=args.nomicro
		print("Operating in test mode for " + args.nomicro + " micrographs")
	else:
		testmode = 0
		print("Operating for the whole dataset")		
	liststar=glob.glob(args.istarpattern)
	outstar = open(args.ostar, 'w')
	
	count = 1
	for starfile in liststar:
		# Control for test mode
		if ( testmode == 1 and count > int(nomicro) ):
			print("Finish test mode")
			break	
		print("Combining " + starfile)
		instar = open(starfile, 'r')
		starlabels = learnstarheader(instar)
		if count == 1: # First star file
			writestarheader(outstar, starlabels)	
		count += 1
		for line in instar:
			record = line.split()
			if len(record)==len(starlabels): # if line looks valid			
				writestarline(outstar, record)
		instar.close()

	outstar.close()
