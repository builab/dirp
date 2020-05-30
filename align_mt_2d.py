#!/usr/bin/env python
# Python script to run 2d alignment and create a combine star file
# Realistically, only the combined star file is needed to keep
# Right now, a lot of hard code but will be refined later
# HB 2020/05/30

import os, sys, argparse, shutil, os.path, glob, string
from shutil import copyfile

# Global variable
iter=6
tau2_fudge=12
thread=4
offset_range=10
helical_outer_diameter=500
sigma_psi=5

# Still hard-code the relion command now
def align2d(starfile, outprefix):		  
	"""Performing align 2d of mt star file using relion"""
	align2d = "relion_refine --i " + starfile + " --o " + outprefix + " --dont_combine_weights_via_disc --no_parallel_disc_io --preread_images --pool 200 --pad 2 --ctf --iter 6 --tau2_fudge 12 --particle_diameter 600 --K 1 --flatten_solvent --oversampling 1 --psi_step 2 --offset_range 10 --offset_step 2 --helical_outer_diameter 500 --sigma_psi 5 --dont_check_norm --norm --scale --j 4"
	print(align2d)
	os.system(align2d)
	
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
	parser = argparse.ArgumentParser(description='Perform 2D alignment of microtubule segment')
	parser.add_argument('--idir', help='Input dir for microtubule star file',required=True)
	parser.add_argument('--odir', help='Output dir for MT star files',required=True)
	parser.add_argument('--min_particles', help='Minium number of particles',required=False,default=10)
	parser.add_argument('--nomicro', help='Test mode for only this number of micrographs',required=False)

	args = parser.parse_args()
	
	min_part = int(args.min_particles)

	if args.nomicro is not None:
		testmode = 1
		nomicro=args.nomicro
		print("Operating in test mode for " + args.nomicro + " micrographs")
	else:
		testmode = 0
		print("Operating for the whole dataset")	

	liststar=glob.glob(args.idir + "/*MT*.star")
	outdir = args.odir

	try:
		os.mkdir( outdir, 0755 );
	except:
		print("Output directory exist")
	
	count = 1
	for starfile in liststar:
		# Control for test mode
		if ( testmode == 1 and count > int(nomicro) ):
			print("Finish test mode")
			break	
		count += 1
		basename = os.path.basename(starfile)
        	basename = string.replace(basename, ".star", "")
		alndir = outdir + "/" + basename
		try:
			os.mkdir( alndir, 0755 );
		except:
			print( alndir + " exists")
		# Check if the min_particle statisfy, otherwise skip
		if starcountparticles(starfile) < min_part:
			print("Skip " + starfile)
			continue
		# Perform alignment
		align2d(starfile, alndir + "/" + basename)
		outstar =  alndir + "/" + basename + "_it006_data.star"
		outmrc = alndir + "/" + basename + "_it006_classes.mrcs"
		try:
			shutil.copyfile(outstar, outdir + "/" + basename + ".star")
		except:	
			print("Error copying file " + outstar)
		try:
			shutil.rmtree(alndir)
		except:
			print("Cannot remove " + alndir)




