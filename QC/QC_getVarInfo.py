#! /usr/bin/env python

# GOAL: print out the Variant info from a vcf and other metrics to the run's json file. 
# This should be called from QC_getRunInfo.sh as the other metrics are passed in from that script.

import sys
import os.path
import re
import json
from optparse import OptionParser

# ------------------------------------------
# ------------- FUNCTIONS ------------------
# ------------------------------------------

# (if flow-space values are present FAO/FRO values will be reported, if not then will use AO/RO values instead)
# based on FAO/FRO (or AO/RO) determine the total depth and allele frequency
# assign GT to each normal variant based on the determined thresholds: 
#@param Takes as input the line from a vcf file and the WT and HOM allele frequency cutoffs. (Could be different for Tumor / Normal pairs)
#@return Returns a tuple with 1. the GT, 2. alternate allele frequency, 3. the alt depth, and 4. the ref depth
def getGTInfo(line, WT_Cutoff, HOM_Cutoff):
	line = line.split("\t")
	alternates = line[4].split(",")
	GT = '.'
	alt_freq = '.'
	alt_depth = '.'
	ref_depth = '.'
	total_depth = 0
	# there should NOT be any multi-allelic calls in the input tumor/vcf files at this stage since these are filtered out upstream.
	# This "if" statement either could be removed, or replaced with an "if/else" to report any unexpected multi-allelic calls seen
	if len(alternates) != 1:
		GT = "ERROR: ", line[0], line[1], "has an alternate allele"
	else:
		info = dict(zip(line[8].split(":"), line[9].split(":"))) #Creates a dictionary with the description as the key, and the actual value as the value.
		try:
			if 'FAO' in info and info['FAO'] != '.': # Ozlem said that the flow space depths are more accurate than the base space depths, so use that if we can.
				alt_depth = info['FAO']
				ref_depth = info['FRO']
			elif 'AO' in info:
				alt_depth = info['AO']
				ref_depth = info['RO']
			elif 'AD' in info: # For older vcf files (3.2), the format is different.
				alt_depth = info['AD'].split(',')[1]
				ref_depth = info['AD'].split(',')[0]
			total_depth = int(alt_depth) + int(ref_depth)
			if total_depth != 0:
				alt_freq = float(alt_depth) / float(total_depth)
				if alt_freq < WT_Cutoff:
					GT = 'WT'
				elif alt_freq >= HOM_Cutoff:
					GT = 'HOM'
				else:
					GT = 'HET'
		except KeyError, ValueError:
			print "Leaving GT as . for this variant:", line
			pass
	return GT, str(alt_freq)[0:5], alt_depth, ref_depth # returns the GT, and first three decimal points of the alt frequencies, and the alt and ref depth.


# ---------------------------------------------------
# ------------- PRGRAM STARTS HERE ------------------
# ---------------------------------------------------
# set up the option parser
parser = OptionParser()

# add the options to parse
parser.add_option('-v', '--filtered_vcf', dest='vcf', help='The filetered.vcf file')
parser.add_option('-j', '--json', dest='json', help='The runs json file')
parser.add_option('-m', '--metrics', dest='metrics', help='The runs metrics to be added to the json file')
parser.add_option('-g', '--gt_cutoffs', dest='gt_cutoffs', nargs=2, type="float", help='WT_Cutoff, HOM_Cutoff. If a variant has FAO/FRO < WT_Cutoff, it is labelled as WT, >= HOM_Cutoff is HOM, in between is HET.')

(options, args) = parser.parse_args()

#check to make sure either ID or name was provided
if(not options.vcf or not options.json or not options.gt_cutoffs or not options.metrics):
	print "USAGE ERROR: all options are required"
	print "use -h for help"
	sys.exit(8)

vcf = open(options.vcf, 'r')
WT_Cutoff = options.gt_cutoffs[0]
HOM_Cutoff = options.gt_cutoffs[1]


if not os.path.isfile(options.json):
	# If the json file doesn't exist for some reason, then add these simple metrics to it
	sample = options.json.split("/")[-4]
	run_name = options.json.split("/")[-2]
	runData = {'sample': sample, 'name': run_name}
else:
	#read in the json file
	jsonData = open(options.json)
	runData = json.load(jsonData)
	# TEMP FOR WALES
	if 'project' in runData and runData['project'] == 'Wales':
		if options.json.split("/")[-4] == "":
			runData['name'] = "_".join(options.json.split("/")[-3:-1]) + "_Run1"
		else:
			runData['name'] = "_".join(options.json.split("/")[-4:])

numHET = 0
numHOM = 0
# Loop through the VCF file and GT all of the variants
for line in vcf:
	if line[0] != "#":
		gtInfo = getGTInfo(line, WT_Cutoff, HOM_Cutoff) 
		if gtInfo[0] == 'HET':
			numHET += 1
		elif gtInfo[0] == 'HOM':
			numHOM += 1
vcf.close()

# Now setup the json dictionary
het_hom = '.'
if numHOM != 0:
	het_hom = numHET / float(numHOM)

if 'run_data' not in runData:
	runData['run_data'] = {}

# add the run's variant info to the json dictionary
runData['run_data']['total_vars'] = numHET + numHOM   # total number of variants for this run
runData['run_data']['num_het'] = numHET  # number of HETs
runData['run_data']['num_hom'] = numHOM  # number of HOMs
runData['run_data']['het_hom'] = het_hom  # HET / HOM ratio

# put the metrics passed in by QC_getRunInfo.sh into a dictionary to add to the runData.
# they are set up as "key1:key2:key3;value1:value2:value3" so we can use dict(zip(keys, values))
metrics = dict(zip(options.metrics.split(';')[0].split(':'), options.metrics.split(";")[1].split(":")))

# Bash can't do float division, so do that division here
try:
	metrics['begin_amp_cov'] = float(metrics['begin_amp_cov']) / float(metrics['num_amps'])
	metrics['end_amp_cov'] = float(metrics['end_amp_cov']) / float(metrics['num_amps'])
	metrics['amp_cov'] = float(metrics['amp_cov']) / float(metrics['num_amps'])
	metrics['base_cov'] = float(metrics['base_cov']) / float(metrics['total_possible_bases'])
except ValueError:
	pass

# check the begin_amp_cov and end_amp_cov separately to not overwrite data incase something went wrong with getting the depths in QC_getRunINfo
if 'begin_amp_cov' in runData['run_data'] and runData['run_data']['begin_amp_cov'] != "" and metrics['begin_amp_cov'] == "":
	# don't overwrite what is already in the JSON file with nothing.
	pass
else:
	runData['run_data']['begin_amp_cov'] = metrics['begin_amp_cov']
if 'end_amp_cov' in runData['run_data'] and runData['run_data']['end_amp_cov'] != "" and metrics['end_amp_cov'] == "":
	# don't overwrite what is already in the JSON file with nothing.
	pass
else:
	runData['run_data']['end_amp_cov'] = metrics['end_amp_cov']

# add the other metrics to the json file.
for key in metrics:
	if key != 'num_amps' and key != 'begin_amp_cov' and key != 'end_amp_cov' and key != 'total_possible_bases':
		try:
			runData['run_data'][key] = float(metrics[key])
		except ValueError:
			runData['run_data'][key] = metrics[key]

#dump the json file
with open(options.json, 'w') as newJSONFile:
	json.dump(runData, newJSONFile, sort_keys=True, indent=4)

