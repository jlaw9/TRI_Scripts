#! /usr/bin/env python2.7

import json
import os
import subprocess
import datetime
import time
import sys

# @param path the path of the json file
# @json_data the json dictionary to be written
def write_json(path, json_data):
	with open(path, 'w') as newJobFile:
		json.dump(json_data, newJobFile, sort_keys=True, indent=4)

# @param systemCall the command to run by bash
# @returns the status of either 0 or 1.
def runCommandLine(systemCall, get_output=False):
	# try running the command through
	if get_output:
		try:
			print "Running subprocess.check_output('%s') at %s" %(systemCall, getTimestamp())
			output = subprocess.check_output(systemCall, shell=True).strip()
			return output
		except AttributeError:
			sys.stderr.write("ERROR: need to use 'python2.7 %s ' to run subprocess.check_output('%s')"%(sys.argv[0], systemCall))
			return ''
	else:
		#run the call and return the status
		print 'Starting %s at %s' %(systemCall, getTimestamp())
		status = os.system(systemCall)
		return(status)

## simple method to get a time stamp
# @returns A timestamp
def getTimestamp():
	#get the time / timestamp
	currentTime = time.time()
	timeStamp = datetime.datetime.fromtimestamp(currentTime).strftime('%Y-%m-%d:%H-%M-%S')

	#return the stamp
	return(timeStamp)

def getDate():
	#get the time / timestamp
	currentTime = time.time()
	date = datetime.datetime.fromtimestamp(currentTime).strftime('%m%d%Y')

	#return the stamp
	return(date)


