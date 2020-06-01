#!/bin/bash
#SBATCH --ntasks=4                   # Number of processes
#SBATCH --partition=cpu
#SBATCH --job-name=xform_avg # Job name
#SBATCH --error=xform_avg.err
#SBATCH --output=xform_avg.out
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --mem-per-cpu=1GB

module load python-2.7.14-gcc-4.9.4-tzsdd67
module load relion/3.0b-jan19

python  dirp/xform_avg_2d_parallel.py --istarpattern "Xform/*.star" --odir AvgTest --ostar AvgTest/segment_average.star --j 4

