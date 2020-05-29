#!/usr/bin/env python3
# Script to parse the particles from Extraction step into individual MT star file
# Main code coming from John Rubinstein createstackfromstar.py
# HB 2020/05/29

import os, sys, argparse, shutil, os.path, glob, string
import numpy as n



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

def starcol_containing_label(starlabels, substring):
    for i, s in enumerate(starlabels):
        if substring in s:
              return i
    return -1

# New function to do exact match of relion label such as _rlnImageCol
def starcol_exact_label(starlabels, label):
    for i, s enumerate(starlabels):
        if label == s:
            return i
    return -1
		

if __name__=='__main__':
    # get name of input starfile, output starfile, output stack file
    parser = argparse.ArgumentParser(description='Create multiple MT star files from a particle star file')
    parser.add_argument('--istar', help='Input particle star file from Extraction',required=True)
    parser.add_argument('--odir', help='Output dir for MT star files',required=True)
    args = parser.parse_args()
    # open input star, output star, output stack
    #with open(args.istar,'r') as instar, open (args.ostar,'w') as avgstar, open(args.ostack,'wb') as outstack:
        # prepare a temporary header for output stack
    instar = open(args.istar, 'r')
    outdir = args.odir
    starlabels = learnstarheader(instar)
    # Column definition
    helicaltubeidcol=starcol_exact_label(starlabels, '_rlnHelicalTubeID')
    print(helicaltubeidcol)
    imagecol=starcol_exact_label(starlabels, '_rlnImageName')
    print(imagecol)
    anglerotcol=25
    angletiltcol=26
    anglepsicol=27
    originxcol=22
    originycol=23
    anglepsipriorcol=4
    angletiltpriorcol=3

    nparts=0 # for every line in starfile
    #for line in avgstar:
    #    record = line.split()
    #    if len(record)==len(starlabels): # if line looks valid
    #        partandstack=record[imagecol].split('@')
    #        helicaltubeid = record[helicaltubeidcol]
    #        basename = os.path.basename(partandstack[1])
	#    basename = str.replace(basename, "_avg.mrcs", "");
	#    print(basename)
            # write a record to the new starfile
    #        instar=open(indir + basename +  ".star", 'r')
    #        for item in instar:
    #            partrecord = item.split()
     #           if len(partrecord) == len(starlabels):
      #              partrecord[originxcol]=record[originxcol]
      #              partrecord[originycol]=record[originycol]
       #            partrecord[angletiltcol]=record[angletiltcol]
        #            partrecord[anglerotcol]=record[anglerotcol]
         #           partrecord[angletiltpriorcol]=record[angletiltpriorcol]
          #          partrecord[anglepsipriorcol]=record[anglepsipriorcol]
           #         writestarline(partstar,partrecord)
           # instar.close()
            #record[imagecol] =str(1).zfill(6)+'@'+ str.replace(starfile, ".star", "_avg.mrcs")
            #print("{:.6f}".format(0))
            #record[psipriorcol] = "{:.6f}".format(0) 
            # Write a trim record to 24 column only
            #writestarline(avgstar,record[:24])
            #instar.close()
            #break

    #partstar.close()
    #avgstar.close()
