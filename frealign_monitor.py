#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 13:19:02 2020
Script sequentially run through a list of mparameters file in Frealign
By monitoring the frealign.log file for normal termination
@author: kbui2
"""

import shutil, time, os.path, sys, re
import subprocess


def launchfrealign():
	''' Launch free align '''
	process = subprocess.Popen(['frealign_run_refine', ''])
	return process

def checkfrealignrunning():
	''' Check if frealign running by check the frealign.log file '''
	log = open('frealign.log', 'r')
	for line in log:
		if re.match('Normal termination of frealign run', line):
			return 0
	return 1		
		
	
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
			if checkfrealignrunning() == 1:
				print('Frealign still alive ' + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
			else:
				print('Frealign done')
				break
			time.sleep(interval)
			
		
		
