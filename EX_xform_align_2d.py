#!/usr/bin/env python
# Python script to transform and average particles after reextraction
# Create also star file of the average microtubule
# Experimental script to create average of each class from 2d file

import os, sys, argparse, shutil, os.path, glob, string
import numpy as np
from numpy import *
from shutil import copyfile

def applytransformation(instar, outputrootname):		  
	"""Apply transformation to a stack"""
	transform2d = "relion_stack_create --i " + instar + " --o " + outputrootname + " --apply_transformation"
	print(transform2d)
	os.system(transform2d)


def splitstarclass(instar): # Very bad with big star file
	"""Split star file into sub-star file base on class number"""
	instarhandle = open(instar, 'r')
	starlabels = learnstarheader(instarhandle)
	classnocol = starcol_exact_label(starlabels, '_rlnClassNumber')
	classlist=np.array([])
	firstline = 1
	for line in instarhandle:
		record = line.split()
		if len(record)==len(starlabels): # if line looks valid
			if firstline == 1:
				data = np.array([record])
			else:
				data = np.append(data, [record], 0)
			firstline = 0
			classlist = np.append(classlist, [int(record[classnocol])])
	instarhandle.close()
	classno = np.unique(classlist)
	index = np.where(classlist == 1)
	print(classlist)
	print(len(index))
					
	
def averagestack(instack, outstack):
	"""Average stack"""
	average2d = "relion_image_handler --i " + instack + " --o " + outstack + " --average"
	print(average2d)
	os.system(average2d)
	
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
	

if __name__=='__main__':
	# get name of input starfile, output starfile, output stack file
	parser = argparse.ArgumentParser(description='Perform transform and average of microtubule segment and generate a new star file for the average')
	parser.add_argument('--istarpattern', help='Input particle star file pattern (Extract/MT*.star)',required=True)
	parser.add_argument('--odir', help='Output dir for transformed and average files',required=True)
	parser.add_argument('--ostar', help='Output star file for the average segments',required=True)
	parser.add_argument('--avgclass', help='Average Class or not',required=False,default=0)
	parser.add_argument('--nomicro', help='Test mode for only this number of micrographs',required=False)

	args = parser.parse_args()
	
	avgclass = int(args.avgclass)
	
	if args.nomicro is not None:
		testmode = 1
		nomicro=int(args.nomicro)
		print("Operating in test mode for " + args.nomicro + " micrographs")
	else:
		testmode = 0
		print("Operating for the whole dataset")	

	liststar=glob.glob(args.istarpattern)
	outdir = args.odir
	
	# Check if the outstar file exists
	if os.path.exists(args.ostar):
		os.remove(args.ostar)
	
	count = 1
	# Create average microtubule
	for starfile in liststar:
		if ( testmode == 1 and count > int(nomicro) ):
			print("Finish test mode")
			break
		count += 1
		basename = os.path.basename(starfile)
		basename = string.replace(basename, ".star", "")
		applytransformation(starfile, outdir + "/" + basename)	
		splitstarclass(outdir + "/" + basename + ".star")
		averagestack(outdir + "/" + basename + ".mrcs", outdir + "/" + basename + "_avg.mrcs")
		
	sys.exit(0)
	# Create average star file
	print ("Create output " + args.ostar + " from " + outdir + "/*.star")
	listxformstar=glob.glob(outdir + "/*.star")
	starfile = listxformstar[0]
	outstar = open(args.ostar, 'w')
	instar = open ( starfile, 'r')
	starlabels = learnstarheader(instar)
	instar.close()
	# write new trimmed (24 col) output starfile header
	writestarheader(outstar, starlabels[:24])   
	
	count = 1
	for starfile in listxformstar:
		if ( testmode == 1 and count > int(nomicro) ):
			print("Finish test mode")
			break	
		count += 1
		print("Reading " + starfile)
		instar = open ( starfile, 'r' )
		# learn the starfile header and column for image names
		starlabels=learnstarheader(instar)
		imagecol = starcol_exact_label(starlabels, '_rlnImageName')
		for line in instar:
			record = line.split()
			if len(record)==len(starlabels): # if line looks valid
				partandstack=record[imagecol].split('@')
				# write a record to the new starfile
				record[imagecol] =str(1).zfill(6)+'@'+ str.replace(starfile, ".star", "_avg.mrcs")
				# Write a trim record to 24 column only
				writestarline(outstar,record[:24])
				instar.close()
				break
	outstar.close()
