#!/usr/bin/env python3
# Python script to transform and average particles after reextraction
# Create also star file of the average microtubule
import os
import shutil
import os.path
import glob
import string
from shutil import copyfile

outdir="wf2/2Daverage/"
reextractdir = "wf2/Xform/"
xformdir = "wf2/Xform2/"
liststar=glob.glob(reextractdir + "*MT*.star")

# Create average microtubule
for starfile in liststar:
	basename = os.path.basename(starfile)
        basename = string.replace(basename, ".star", "")
	#print(basename)
	outstar =  basename + ".star"
	transform2d = "relion_stack_create --i " + reextractdir  + outstar + " --o " + xformdir + basename + " --apply_transformation"
	average2d = "relion_image_handler --i " + xformdir + basename + ".mrcs " + " --o " + xformdir + basename + "_avg.mrcs --average"
	print(transform2d)
	os.system(transform2d)
	print(average2d)
	os.system(average2d)
	#shutil.copyfile(alndir + outstar, outdir + outstar)
	#shutil.copyfile(alndir + outmrc, outdir + outmrc)
	#shutil.rmtree(alndir)
	#count = count + 1
	#if count > 4:
	#	break
	#break
