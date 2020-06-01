#!/usr/bin/env python
# Script to read the segment_average alignment file and then propagate the alignment to the particles
# This script still doesn't take into a bit of rotation between particles & average
# HB 2020/05/31

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

def starcol_exact_label(starlabels, label):
	"""New function to do exact match of relion label such as _rlnImageCol"""
	for i, s in enumerate(starlabels):
		record=s.split()
		if label == record[0]:
			return i
	return -1

if __name__=='__main__':
	# get name of input starfile, output starfile, output stack file
	parser = argparse.ArgumentParser(description='Create particle star files from the segment average')
	parser.add_argument('--istar', help='Input segment average star file',required=True)
	parser.add_argument('--imtstardir', help='Directory of microtubule star directory (E.g.: Avg2D)',required=True)
	parser.add_argument('--ostar', help='Output particle star file',required=True)
	parser.add_argument('--nomicro', help='Test mode for only this number of micrograph',required=False)
	
	if args.nomicro is not None:
		testmode = 1
		nomicro=args.nomicro
		print("Operating in test mode for " + args.nomicro + " micrographs")
	else:
		testmode = 0
		print("Operating for the whole dataset")
		
	avgstar = open(args.istar , 'r')
	indir = args.imtstardir
	partstar = open(args.ostar, 'w')
	
	starlabels = learnstarheader(avgstar)
	writestarheader(partstar, starlabels)
	# Column definition
	helicaltubeidcol = starcol_exact_label(starlabels, '_rlnHelicalTubeID')
	imagecol =  starcol_exact_label(starlabels, '_rlnImageName')
	anglerotcol = starcol_exact_label(starlabels, '_rlnAngleRot')
	angletiltcol = starcol_exact_label(starlabels, '_rlnAngleTilt')
	anglepsicol = starcol_exact_label(starlabels, '_rlnAnglePsi')
	originxcol = starcol_exact_label(starlabels, '_rlnOriginX')
	originycol = starcol_exact_label(starlabels, '_rlnOriginY')
	anglepsipriorcol = starcol_exact_label(starlabels, '_rlnAnglePsiPrior')
	angletiltpriorcol = starcol_exact_label(starlabels, '_rlnAngleTiltPrior')

	count = 1 # for every line in starfile
	for line in avgstar:
		record = line.split()
		if len(record)==len(starlabels): # if line looks valid
			# Control for test mode
			if ( testmode == 1 and count > int(nomicro) ):
				print("Finish test mode")
				break
			count += i
			partandstack=record[imagecol].split('@')
			helicaltubeid = record[helicaltubeidcol]
			basename = os.path.basename(partandstack[1])
			basename = str.replace(basename, "_avg.mrcs", "");
			print(basename)
			# write a record to the new starfile
			instar=open(indir + basename +  ".star", 'r')
			for item in instar:
				partrecord = item.split()
				if len(partrecord) == len(starlabels):
					partrecord[originxcol]=record[originxcol]
					partrecord[originycol]=record[originycol]
					partrecord[anglepsicol]=record[anglepsicol]
					partrecord[angletiltcol]=record[angletiltcol]
					partrecord[anglerotcol]=record[anglerotcol]
					partrecord[angletiltpriorcol]=record[angletiltpriorcol]
					partrecord[anglepsipriorcol]=record[anglepsipriorcol]
					writestarline(partstar,partrecord)
			instar.close()
	partstar.close()
	avgstar.close()
