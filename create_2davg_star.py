#!/usr/bin/env python3
# Python script to generate star file of 2d average 
# Create also star file of the average microtubule
# Update 2020/06/15

import os, sys, argparse, shutil, os.path, glob, string
from shutil import copyfile
	
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
	parser = argparse.ArgumentParser(description='Generate average star file')
	parser.add_argument('--istar', help='Original input particle star file',required=True)
	parser.add_argument('--ostar', help='Output star file for the average segments',required=True)
	parser.add_argument('--avgdir', help='Directory for the average segments',required=True)
	parser.add_argument('--trim', help='Only print the first 24 records',required=False, default=0)
	parser.add_argument('--resetpsiprior', help='Reset angle psi prior to 0',required=False, default=0)
	parser.add_argument('--nomicro', help='Test mode for only this number of micrographs',required=False)

	args = parser.parse_args()
	
	dotrim = int(args.trim)
	resetpsiprior = int(args.resetpsiprior)
	
	indir = args.avgdir
	
	if args.nomicro is not None:
		testmode = 1
		nomicro=args.nomicro
		print("Operating in test mode for " + args.nomicro + " micrographs")
	else:
		testmode = 0
		print("Operating for the whole dataset")	
	
	# Check if the outstar file exists
	if os.path.exists(args.ostar):
		os.remove(args.ostar)
	

	# Create 2d average star file
	instar = open(args.istar, 'r')
	starlabels = learnstarheader(instar)
	helicaltubeidcol = starcol_exact_label(starlabels, '_rlnHelicalTubeID')
	imagecol = starcol_exact_label(starlabels, '_rlnImageName')
	psipriorcol = starcol_exact_label(starlabels, '_rlnAnglePsiPrior')
	outstar = open(args.ostar, 'w')
	# Write trim star header
	if trim > 0:
		lastcol = 24
	else:
		lastcol = len(starlabels)
		
	writestarheader(outstar, starlabels[:lastcol])   

	
	currid=''
	currimage=''
	newmt=0
	count = 1
	for line in instar:
		# Control for test mode
		if ( testmode == 1 and count > int(nomicro) ):
			print("Finish test mode")
			break		
		record = line.split()
		if len(record)==len(starlabels): # if line looks valid			
			partandstack=record[imagecol].split('@')
			imagename=partandstack[1]
			basename = os.path.basename(imagename)
			basename = str.replace(basename, ".mrcs", "");	
			#print(basename)
			helicaltubeid = record[helicaltubeidcol]
			if currimage == imagename:
				if currid == helicaltubeid:
					newmt = 0
				else:
					newmt = 1		
			else: 		
				newmt = 1
			currid = helicaltubeid
			currimage = imagename
			if newmt == 1:
				if os.path.exists(indir + '/' + basename + '_MT' + helicaltubeid + '_avg.mrcs'):
					print(basename + "_MT" + helicaltubeid)
					if dotrim == 1:
						record[psipriorcol] = "{:.6f}".format(0)
					record[imagecol] =str(1).zfill(6)+'@'+ indir + '/' + basename + '_MT' + helicaltubeid + '_avg.mrcs'
					# Write a trim record to 24 column only
					writestarline(outstar,record[:lastcol])
				count += 1
			newmt = 0

	instar.close()
	outstar.close()
