#!/usr/bin/env python3
# Python script to run 2d alignment and create a combine star file
# Realistically, only the combined star file is needed to keep

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
	#os.system(align2d)

if __name__=='__main__':
	# get name of input starfile, output starfile, output stack file
	parser = argparse.ArgumentParser(description='Perform 2D alignment of microtubule segment')
	parser.add_argument('--idir', help='Input dir for microtubule star file',required=True)
	parser.add_argument('--odir', help='Output dir for MT star files',required=True)
	parser.add_argument('--ostar', help='Output combined star file',required=True)
	parser.add_argument('--nomicro', help='Test mode for only this number of micrograph',required=False)

	args = parser.parse_args()

	if args.nomicro is not None:
		testmode = 1
		nomicro=args.nomicro
		print("Operating in test mode for " + args.nomicro + " micrographs")
	else:
		testmode = 0
		print("Operating for the whole dataset")	

	liststar=glob.glob(args.idir + "/*MT*.star")
	outdir = args.odir
	outstar = args.outstar

	if os.path.exists(outdir):
		shutil.rmtree( outdir )
	os.mkdir( outdir, 0755 );

	for starfile in liststar:
		# Control for test mode
		if ( testmode == 1 and count > int(nomicro) ):
			print("Finish test mode")
			break	
		basename = os.path.basename(starfile)
        	basename = string.replace(basename, ".star", "")
		alndir = outdir + basename
		try:
			os.mkdir( alndir, 0755 );
		except:
			print( alndir + " exists")
		# Perform alignment
		align2d(starfile, alndir + "/" + basename)
		outstar =  alndir + "/" + basename + "_it006_data.star"
		outmrc = alndir + "/" + basename + "_it006_classes.mrcs"
		shutil.copyfile(outstar, outdir + "/" + basename + "_it006_data.star")
		shutil.copyfile(outmrc, outdir + "/" + basenaem + "_it006_classes.mrcs")
		shutil.rmtree(alndir)
		count += 1


# Combine star file
count=0
if os.path.exists(cbndir + cbnStar):
	os.remove(cbndir + cbnStar)

for starfile in liststar:
	basename = os.path.basename(starfile)
	basename = string.replace(basename, ".star", "")
	grep= "grep mrc " + cbndir + basename + "_it006_data.star" + " >> " + cbndir + cbnstar
	count=count+1
	if count==1:	
		shutil.copyfile(cbndir + basename + "_it006_data.star", cbndir + cbnstar)
	else:
		print(grep)
		os.system(grep)

