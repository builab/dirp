#!/usr/bin/env python3
# Script originally from John Rubinstein
# Temporary modify to work with Star3_1

import os, sys, argparse
import numpy as n

class LazyMRC:
    def __init__(self, fname, shape, dtype, idx):
        self.fname = fname
        self.shape = (int(shape[0]), int(shape[1]))
        self.idx = idx
        self.dtype = dtype
        self.length = (n.dtype(dtype).itemsize) * shape[0] * shape[1]
        self.offset = 1024 + idx * self.length
    def get(self):
        with open(self.fname) as f:
            f.seek(self.offset)
            data = n.reshape(n.fromfile(f, dtype=self.dtype, count= n.prod(self.shape)), self.shape, order='F')
        return data
    def view(self):
        return self.get()

def get_dtype(hdr):
    return {1:n.int16, 2:n.float32} [hdr['datatype']]

def readMRClazy (fname):
    hdr = readMRCheader(fname)
    shape = (hdr['nx'],hdr['ny'])
    num = hdr['nz']
    dtype = get_dtype(hdr)
    lazy_data = [LazyMRC(fname, shape, dtype, idx) for idx in range(num)]
    return lazy_data

def readMRCmemmap (fname, inc_header=False):
    """ Read a memory mapped MRC file and header """ 
    hdr = readMRCheader(fname)
    nx = hdr['nx']
    ny = hdr['ny']
    nz = hdr['nz']
    dtype = {0:n.int8, 1:n.int16, 2:n.float32, 6:n.uint16}[hdr['datatype']] 
    mm = n.memmap(fname, dtype, 'r',
                    offset=1024, shape=(nx, ny, nz), order='F')
    return (mm, hdr) if inc_header else mm



def writeMRCheader(filename,nxyz=0,dmin=0,dmax=0,dmean=0,mode=2,psize=1):
    """ Write an MRC header given key parameters about the MRC image or stack """
    filename.seek(0) # Go to the beginning of the file
    header = n.zeros(256, dtype=n.int32) # 1024 byte header (integer values)
    header_f = header.view(n.float32)    # floating point values
    header[:3] = nxyz
    header[3] = mode   # mode, 2 = float32 datatype
    header[7:10] = nxyz  # mx, my, mz (grid size)
    header_f[10:13] = [psize*i for i in nxyz] # xlen, ylen, zlen
    header_f [13:16] = 90.0 # CELLB
    header [16:19] = [1,2,3] # axis order
    header_f [19:22] = [dmin, dmax, dmean] # data stats
    header [52] = 542130509   # 'MAP ' chars
    header [53] = 16708
    header.tofile(filename)

def writeMRCsection(data,filename):
    """ Write a section of a stack to an MRC file """
    n.require(n.reshape(data, (-1,), order='F'), dtype=n.float32).tofile(filename)

def readMRCheader(fname):
    hdr = None
    with open(fname) as f:
        hdr = {}
        header = n.fromfile(f, dtype=n.int32, count=256)
        header_f = header.view(n.float32)
        [hdr['nx'], hdr['ny'], hdr['nz'], hdr['datatype']] = header[:4]
        [hdr['xlen'],hdr['ylen'],hdr['zlen']] = header_f[10:13]
        # print "Nx %d Ny %d Nz %d Type %d" % (nx, ny, nz, datatype)
    return hdr

def learnstarheader(infile):
    """Learn which column contains which information from an already open starfile"""
    infile.seek(0) # Go to the beginning of the starfile
    doneheader = False
    doneprelabels = False
    headerlabels = []
    headeroptics = []
    doneoptics = True
    while not doneprelabels:
        line=infile.readline()
        # Check if star 3.1 format
        if line.startswith('data_optics'):
            doneoptics = False
        if line.startswith('data_particles'):
            doneoptics = True
        if line.startswith('loop_') & doneoptics == True:
            doneprelabels = True # read until 'loop_'
        headeroptics += [line]
    while not doneheader:
        line=infile.readline()
        if not line.startswith('_'): # read all lines the start with '_'
            doneheader = True
        else:
            headerlabels += [line] 
    infile.seek(0) # return to beginning of starfile before return
    return headeroptics, headerlabels

def is_star3_1(infile):
    """Learn starfile is 3.1 or not"""
    infile.seek(0) # Go to the beginning of the starfile
    is_star3_1 = False
    doneheader = False
    while not doneheader:
        line=infile.readline()
        # Check if star 3.1 format
        if line.startswith('data_optics'):
            doneheader = True
            is_star3_1 = True
		
    infile.seek(0) # return to beginning of starfile before return
    return is_star3_1

def writestarheader(outfile, headeroptics, headerlabels):		  
    """With an already opened starfile write a header"""
    for label in headeroptics:
        outfile.write(label)
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

if __name__=='__main__':
    # get name of input starfile, output starfile, output stack file
    parser = argparse.ArgumentParser(description='Create a single stack from particles specified in a starfile')
    parser.add_argument('--istar', help='Input starfile from which to get particle information',required=True)
    parser.add_argument('--ostar', help='Output starfile with corrected file name for particle stack',required=True)
    parser.add_argument('--ostack', help='Output particle stack containing only particles in starfile',required=True)
    args = parser.parse_args()
    # open input star, output star, output stack
    print("WARNING: Only test compatibility with Relion 3.1")
    
    with open(args.istar,'r') as instar, open (args.ostar,'w') as outstar, open(args.ostack,'wb') as outstack:
        # prepare a temporary header for output stack
        nxyz=n.array([0,0,0])
        dmin=n.inf;dmax=-n.inf;dmean=0
        writeMRCheader(outstack,nxyz,dmin,dmax,dmean,mode=2)
        # learn the starfile header and column for image names
        staroptics, starlabels=learnstarheader(instar)
        imagecol = starcol_containing_label(starlabels, 'ImageName')
        if imagecol == -1: print 'The starfile does not contain a column with ImageName', exit()
        # write output starfile header
        writestarheader(outstar, staroptics, starlabels)
        nparts=0 # for every line in starfile
        for line in instar:
            record = line.split()
            if len(record)==len(starlabels): # if line looks valid
                nparts+=1
                partandstack=record[imagecol].split('@')
                print "reading particle",partandstack[0],"from stack",partandstack[1]
                # write the relevant particle image to the output stack
                # lazy MRC
                frame = readMRClazy(partandstack[1])
                cursection = frame[int(partandstack[0])-1].view()
                # write the particle image to the new stack
                writeMRCsection(cursection,outstack)
                dmean+=n.sum(cursection)
                dmin = min(n.min(cursection),dmin)
                dmax = max(n.max(cursection),dmax)
                # write a record to the new starfile
                record[imagecol]=str(nparts).zfill(6)+'@'+str(args.ostack)
                writestarline(outstar,record)
        # write final MRC header
        nxyz=n.array([cursection.shape[0],cursection.shape[1],nparts])
        dmean=dmean/(nxyz[0]*nxyz[1]*nxyz[2])
        writeMRCheader(outstack,nxyz,dmin,dmax,dmean,mode=2)

