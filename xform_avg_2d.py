#!/usr/bin/env python3
# Python script to transform and average particles after reextraction
# Create also star file of the average microtubule
import os, sys, argparse, shutil, os.path, glob, string
from shutil import copyfile

def applytransformation(instar, outputrootname):		  
	"""Apply transformation to a stack"""
	transform2d = "relion_stack_create --i " + instar + " --o " + outputrootname + " --apply_transformation"
	print(transform2d)
	os.system(transform2d)
	
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

def readstarline(infile):
	"""Read a record (line) from an already open starfile and return XXX"""
	line=infile.readline()
	records = line.split()
	return records
	
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
	parser.add_argument('--nomicro', help='Test mode for only this number of micrographs',required=False)

	args = parser.parse_args()
	
	if args.nomicro is not None:
		testmode = 1
		nomicro=args.nomicro
		print("Operating in test mode for " + args.nomicro + " micrographs")
	else:
		testmode = 0
		print("Operating for the whole dataset")	

	liststar=glob.glob(args.istarpattern)
	outdir = args.odir
	outstar = open(args.ostar, 'w')
	
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
		averagestack(outdir + "/" + basename + ".mrcs", outdir + "/" + basename + "_avg.mrcs")
		
	# Create average star file
	listxformstar=glob.glob(outdir + "/*.star")
   	starfile = listxformstar[0]
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
