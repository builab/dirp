#!/usr/bin/env python3
# Python script to run 2d alignment and create a combine star file
# Realistically, only the combined star file is needed to keep
import os
import shutil
import os.path
import glob
import string
from shutil import copyfile

liststar=glob.glob("wf2/*MT*.star")
#count=1
alndir = "wf2/Class2D/"
outdir="wf2/2Daverage/"
cbnstar = "particles_2dalign.star"

if os.path.exists(alndir):
	shutil.rmtree(alndir)

if os.path.exists( outdir):
	shutil.rmtree( outdir )

os.mkdir( outdir, 0755 );

for starfile in liststar:
	os.mkdir( alndir, 0755 );
	basename = os.path.basename(starfile)
        basename = string.replace(basename, ".star", "")
	#print(basename)
	align2d = "relion_refine --i " + starfile + " --o " + alndir + basename + " --dont_combine_weights_via_disc --no_parallel_disc_io --preread_images --pool 200 --pad 2 --ctf --iter 6 --tau2_fudge 12 --particle_diameter 600 --K 1 --flatten_solvent --oversampling 1 --psi_step 2 --offset_range 10 --offset_step 2 --helical_outer_diameter 500 --sigma_psi 5 --dont_check_norm --norm --scale --j 4"
	print(align2d)
	os.system(align2d)
	outstar =  basename + "_it006_data.star"
	outmrc = basename + "_it006_classes.mrcs"
	shutil.copyfile(alndir + outstar, outdir + outstar)
	shutil.copyfile(alndir + outmrc, outdir + outmrc)
	shutil.rmtree(alndir)
	#count = count + 1
	#if count > 4:
	#	break

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

