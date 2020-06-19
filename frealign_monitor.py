#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 13:19:02 2020
Script sequentially run through a list of mparameters file in Frealign
Doesn't work. Check frealign_kill script later
@author: kbui2
"""

import shutil, time, os.path, sys
import subprocess


def launchfrealign():
	''' Launch free align '''
	process = subprocess.Popen(['frealign_run_refine', ''])
	return process

def checkfrealignrunning():
	if os.path.exists('scratch/pid.log' ):
		return 1
	else:
		return 0
	
if __name__ == "__main__":
	# Checking interval
	interval = 300
	
	if len(sys.argv) < 2:
		print "Usage: python ./frealign_monitor.py mparameter1 mparameter2 mparametere3"
	
	mparamlist = sys.argv[1:]

	for i in range(len(mparamlist)):
		shutil.copy(mparamlist[i], 'mparameters')
		process = launchfrealign()
		while True:
			if checkprocessrunning() == 1:
				print('Frealign still alive ' + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
			else:
				print('Frealign done')
				break
			time.sleep(interval)
			
		
