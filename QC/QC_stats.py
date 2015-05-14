#! /usr/bin/env python

# GOAL: Take the information about each sample and generate the statistics for the QC metrics spreadsheet.
# GOAL: Find each runs json, and each samples QC json file to generate the statistics for the QC metrics spreadsheet.
# Writes the output as key:value so the QC_generateSheet.py will be able to use them.
import sys
import re
import math
import os
import fnmatch
import json
from optparse import OptionParser

# -----------------------------------------------------------------------------------------
# ------------------------------ FUNCTIONS DEFINED HERE -----------------------------------
# -----------------------------------------------------------------------------------------

# @param project_path the path to the project for which you want to generate the spreadsheet
# @param total_possible_bases the total_possible_bases for this QC project (either the bases in the project_bed file, or the project_bed intersected with the CDS bed. If exome data, only chr1 should be used
# @returns a json file ready for QC_generateSheets
def main(project_path, total_bases):
	# This variable is used to calculate the perc_avail_bases by calcStats. Make it global so you don't have to pass it around
	global total_possible_bases 
	total_possible_bases = total_bases
	
	# get both lists containing the json files.
	run_files, QC_files = findRunsJsons(project_path)

	all_samples_run_info = getAllRunInfo(run_files)
	all_samples_QC_info = getAllQCInfo(QC_files)
	# Now match the run and QC info
	all_samples_info = matchQCRunInfo(all_samples_run_info, all_samples_QC_info)
	# Now that the Run data and all of the QC_comparisons are pulled together for each sample, get the Reference Runs 
	all_samples_info = getAllStats(all_samples_info)
	# now write all of the metrics to a materJson file ready to be output by QC_generateSheets
	return generateMasterJson(all_samples_info)

# @param project_path the path to the project for which you want to generate the spreadsheet
# @returns a json file ready for QC_generateSheets
def main_runs_only(project_path):
	# find the QC_files, and return the dictionary of QC_info
	run_files, QC_files = findRunsJsons(project_path)
	return getAllRunInfo(run_files)

# @param project_path the path to the project for which you want to generate the spreadsheet
# @returns a json file ready for QC_generateSheets
def main_QC_only(project_path):
	# find the QC_files, and return the dictionary of QC_info
	run_files, QC_files = findRunsJsons(project_path)
	return getAllQCInfo(QC_files)

# Find and load info all of the run's json files
# @param the path to the project to look for the json files in
# @param the json file pattern used to find the .json files (i.e. *.json* or *_QC.json*)
# @return returns a dictionary containing each runs metrics
def findRunsJsons(project_path):
	QC_files = []
	run_files = []
	# recurse through the project_path and find the json files
	for root, dirnames, filenames in os.walk(project_path):
		for filename in fnmatch.filter(filenames, "*.json*"):
			# append all of the QC sample run comparison json files to this list
			if re.search("_QC.json", filename):
				QC_files.append(os.path.join(root, filename))
			else:
				run_files.append(os.path.join(root, filename))
	
	return run_files, QC_files	

# @param run_files a list of json files
# @return a dictionary containing the combined runs files
def getAllRunInfo(run_files):
	# the master dictionary containing each sample's runs and each run's info
	all_run_data = {}
	
	for jsonFile in run_files:
		# load the json file's data
		jsonData = json.load(open(jsonFile))
		# load this runs data into the master dictionary.
		# The dictionary is organized by sample, then by name, then by run data
		if 'sample' in jsonData:
		#	print "%s is being skipped. No 'sample' key found"%jsonFile
		#else:
			if jsonData['sample'] not in all_run_data:
				all_run_data[jsonData['sample']] = {'runs': {}}
			# the dictionary inside of run_data should hold all of this runs QC metrics such as % polyclonality and such
			# if this json file is not a run json, then don't load it.
			if 'json_type' in jsonData and jsonData['json_type'] != 'run' and jsonData['json_type'] != 'merged':
				pass
			elif 'run_data' not in jsonData:
				print "%s has no run_data"%jsonFile
			else:
				print "	%s run data is being loaded."%jsonFile
				name = 'name'
				if 'name' not in jsonData:
					name = 'run_name'
				all_run_data[jsonData['sample']]['runs'][jsonData[name]] = dict(jsonData.items() + jsonData['run_data'].items())
				if 'sample' not in all_run_data[jsonData['sample']]['runs'][jsonData[name]]:
					all_run_data[jsonData['sample']]['runs'][jsonData[name]]['sample'] = jsonData['sample']
					all_run_data[jsonData['sample']]['runs'][jsonData[name]]['run_num'] = jsonData[name]
#		except KeyError:
#			print jsonFile, jsonData
#			sys.exit(8)

	
	return all_run_data

# @param QC_files a list of _QC.json files
# @return a dictionary containing the combined samples QC_comparisons
def getAllQCInfo(QC_files):
	all_samples_QC_info = {}
	# For each sample, find the best ref run.  
	# A ref run must have > 95 mean read length, a high median read count, and above 50% available bases. 
	# Then choose the ref  based on best error rate
	for jsonFile in QC_files:
		# load the QC json file's data
		jsonData = json.load(open(jsonFile))
		# load this samples QC data into the master dictionary.
		# the dictionary inside of run_data should hold all of this runs QC metrics such as % polyclonality and such
		if 'QC_comparisons' in jsonData:
			if 'sample' in jsonData:
				all_samples_QC_info[jsonData['sample']] = jsonData['QC_comparisons']
			else:
				all_samples_QC_info[jsonData['sample_name']] = jsonData['QC_comparisons']
	
	return all_samples_QC_info


# @param all_samples_run_info a dictionary of each run's json file info
# @param all_samples_QC_info a dictionary of all_samples_QC_info from the _QC.json files found
# @returns a dictionary with the samples matched
def matchQCRunInfo(all_samples_run_info, all_samples_QC_info):
	for sample in all_samples_QC_info:
		# load this samples QC data into the master dictionary.
		# the dictionary inside of run_data should hold all of this runs QC metrics such as % polyclonality and such
		if sample in all_samples_run_info:
			all_samples_run_info[sample]['QC_comparisons'] = all_samples_QC_info[sample]
	
	return all_samples_run_info 

	
# get all of the error rate statistics for the project. Also find the ref runs for each sample
# @param all_samples_info a dictionary containing the QC_comparisons of each sample, and each run's QC metrics
# @returns the all_samples_info also containing the reference or best runs for each sample
def getAllStats(all_samples_info):
	tumor = normal = False
	# now get the REF runs.
	for sample in all_samples_info:
		# check to see if we're dealing with tumor normal or not in this sample
		tumor_dates = []
		normal_dates = []
		for run, run_data in all_samples_info[sample]['runs'].iteritems():
			if run != 'QC_comparisons':
				try:
					if 'run_type' not in run_data:
						if re.search("T-", run):
							all_samples_info[sample]['runs'][run]['run_type'] = "tumor"
							tumor = True
						else:
							all_samples_info[sample]['runs'][run]['run_type'] = "normal"
							normal = True
					elif run_data['run_type'] == "tumor":
						tumor = True
						tumor_dates.append(run_data['name'], run_data['run_date'])
					elif run_data['run_type'] != "tumor":
						normal = True
						normal_dates.append(run_data['name'], run_data['run_date'])
				except KeyError:
					print run_data
					print "%s has no run_type, or no run_date"%run
					sys.exit(8)
		# end of runs in sample for loop
	
		#if len(tumor_dates) > 0 and len(normal_dates) > 0:
		if tumor and normal:
			all_samples_info[sample]['tumor_normal'] = True
			tumor_ref = normal_ref = ''
			# this is a tumor normal sample. Check to see if we have same day pairs
			# same_day_pairs is a set of the same_day_pairs
			same_day_pairs = findSameDayPairs(all_samples_info[sample]['runs'], tumor_dates, normal_dates)
			if len(same_day_pairs) > 0:
				print same_day_pairs
				sys.exit()
				# get the QC_comparisons dictionaries from the sample's runs of these same_day_pairs so we can pass that to findBestRunComp
				same_day_pairs_runs = {'QC_comparisons': {}}
				for pair in same_day_pairs:
					same_day_pairs_runs['QC_comparisons'][pair] = runs['QC_comparisons'][pair]
				# Now find the best same day pair, and use those tumor and normal runs as the REF runs
				avail_runs, bestRunComp = findBestRunComp(same_day_pairs_runs)
				# Store the ref runs, and the same day pairs in the sample's dictionary
				all_samples_info[sample]['normal_ref'] = bestRunComp.split("vs")[0]
				all_samples_info[sample]['tumor_ref'] = bestRunComp.split("vs")[1]
				all_samples_info[sample]['same_day_paris'] = same_day_pairs
			else:
				print "I'm here in ", sample
				# if there arent any same day pairs, then find the best pairs and display those.
				avail_runs, bestRunComp = findBestRunComp(all_samples_info[sample])
				all_samples_info[sample]['normal_ref'] = bestRunComp.split("vs")[0]
				all_samples_info[sample]['tumor_ref'] = bestRunComp.split("vs")[1]
	
		else:
			# this is a germline or only tumor sample. We will treat them the same.
			all_samples_info[sample]['tumor_normal'] = False
			avail_runs, bestRunComp = findBestRunComp(all_samples_info[sample])
			all_samples_info[sample]['ref_run'] = getBetterRun(all_samples_info[sample]['runs'], avail_runs, bestRunComp)
	
	# Now get the error rate statistics
	if tumor and normal:
		# If there are tumor and normal runs in this project, then we will need to calculate the error rates for tn, nn, and tt combinations
		tn_error_rates = []
		nn_error_rates = []
		tt_error_rates = []
		# loop through all of the samples
		for sample in all_samples_info:
			# Loop through each Qc comparison in this sample, and get the error metrics for each one.
			for QC_comp, metrics in all_samples_info[sample]['QC_comparisons'].iteritems():
				if metrics['run1_type'] == 'normal' and metrics['run2_type'] == 'tumor':
					tn_error_rates.append(metrics['error_rate'])
				elif metrics['run1_type'] == 'normal':
					nn_error_rates.append(metrics['error_rate'])
				elif metrics['run1_type'] == 'tumor':
					tt_error_rates.append(metrics['error_rate'])
	
		all_samples_info['metrics'] = {}
		# Now get the mean error rate and Standard deviation for each type of comparison and store those in the all_sample_run_info dictionaries
		all_samples_info['metrics']['tn_mean_error_rate'], all_samples_info['metrics']['tn_SD'] = getStandardDeviation(tn_error_rates)
		#print nn_error_rates
		all_samples_info['metrics']['nn_mean_error_rate'], all_samples_info['metrics']['nn_SD'] = getStandardDeviation(nn_error_rates)
		#print tn_error_rates
		all_samples_info['metrics']['tt_mean_error_rate'], all_samples_info['metrics']['tt_SD'] = getStandardDeviation(tt_error_rates)
	
	else:
		error_rates = []
		# loop through all of the samples
		for sample in all_samples_info:
			# Loop through each Qc comparison in this sample, and get the error metrics for each one.
			for QC_comp, metrics in all_samples_info[sample]['QC_comparisons'].iteritems():
				error_rates.append(metrics['error_rate'])
	
		all_samples_info['mean_error_rate'], all_samples_info['SD'] = getStandardDeviation(error_rates)
	
	return all_samples_info

# @param runs a dictionary of the runs of the sample
# @param tumor_dates list of tumor runs dates
# @param normal_dates list of normal runs dates
# @returns a set of the same day pairs
def findSameDayPairs(runs, tumor_dates, normal_dates):
	same_day_pairs = set()
	#print runs
	for run, run_info in runs.iteritems():
		#print run, run_info
		# find the same day pairs, and add those to be written
		if isinstance(run_info, dict) and 'run_type' in runs[run]:
			if run_info['run_type'] == "tumor":
				for date in normal_dates:
					# if this date matches any of the normal dates, then add it to the same_day_pairs list
					if run_info['run_date'] == date[1]:
						# put the normal run first to follow the pattern we've been following. Because same_day_pairs is a set, we won't add the same normal / tumor pair twice
						same_day_pairs.add("%svs%s"%(date[0], run_info['name']))
			else:
				for date in tumor_dates:
					# if this date matches any of the normal dates, then add it to the same_day_pairs list
					if run_info['run_date'] == date[1]:
						same_day_pairs.add( "%svs%s"%(run_info['name'], date[0]))
	return same_day_pairs			


# Finds the best run out of all of the samle's runs
# @param runs takes a dictionary of the runs in a sample
# @returns the name of the best run.
def findBestRunComp(sample_dict):
	# first, filter the available runs by mean (or median) read length, and by % available read length
	read_length_filter = 95
	median_coverage_overall_filter = 100
	avail_runs = [] # a list of runs that pass the filters (available to be reference runs for the spreadsheet)
	has_run_data = True # variable to make sure this sample has run_data
	# keep looping until at least one run passes the filters.
	while len(avail_runs) == 0:
		has_run_data = False
		# loop through the runs in the sample, and put all of the runs that pass the filters in teh avail_runs list
		for run, run_data in sample_dict['runs'].iteritems():
			# get the read_length
#			if isinstance(runs[run], dict) and 'run_data' in runs[run]:
			has_run_data = True
			if 'median_read_length' in run_data:
				read_length = run_data['median_read_length']
			elif 'mean_read_length' in run_data:
				read_length = run_data['mean_read_length']
			elif 'median' in run_data:
				read_length = run_data['median']
			elif 'mean' in run_data:
				read_length = run_data['mean']
			
			# get the median_coverage_overall
			if 'median_coverage_overall' in run_data:
				median_coverage_overall = run_data['median_coverage_overall']
			elif 'median_read_coverage' in run_data:
				median_coverage_overall = run_data['median_read_coverage']

	
			# If this run passes the filters, then it is a candidate for the ref run
			if int(read_length) > read_length_filter and int(median_coverage_overall) > median_coverage_overall_filter:
				#avail_runs.append({run: run_info})
				avail_runs.append(run)
		# Lower the filters so if no runs passed the filters, a run could pass next pass through the while loop
		read_length_filter -= 5
		median_coverage_overall_filter -= 10
	
	if not has_run_data:
		for run, run_info in runs.iteritems():
			if isinstance(run_info, dict) and run != "QC_comparisons":
				print "%s has no run data"%run
				#avail_runs.append({run: run_info})
				avail_runs.append(run)
	# Now find the lowest error rate of the available run combos
	bestRunComp = ""
	bestErrorRate = 1.0
	# we also need to fileter for perc_avail_bases here
	perc_avail_bases_filter = .6
	# Keep looping until the bestRunComp is found (in case no runs have perc_avail_bases greater than .6)
	while bestRunComp == "":
		for run_comp, error_metrics in sample_dict['QC_comparisons'].iteritems():
			#print 'HERE YO', sample_dict
			for run in avail_runs:
				if re.search(run, run_comp):
					if error_metrics['error_rate'] < bestErrorRate and error_metrics['perc_avail_bases'] > perc_avail_bases_filter:
						bestRunComp = run_comp
						bestErrorRate = error_metrics['error_rate']
		perc_avail_bases_filter -= .05
	
	return avail_runs, bestRunComp
				

# @param runs the dictionary of error metrics for the sample
# @param avail_runs the list of available runs
# @param bestRunComp the name of the bestRunComp (i.e. Run1vsRun2)
# @returns the better of the two runs
def getBetterRun(runs, avail_runs, bestRunComp):
	# get the run names from the bestRunComparison (from Run1vsRun2)
	run1_name = bestRunComp.split('vs')[0]
	run2_name = bestRunComp.split('vs')[1]
	
	# find which of the best error rate runs are available
	if run1_name in avail_runs and run2_name in avail_runs:
		# if both are avail_runs, find the best run out of the two in the comparison by choosing the one with the least errors
		run1Errs = runs['QC_comparisons'][bestRunComp]['HET_WT'] + runs['QC_comparisons'][bestRunComp]['HOM_WT'] + runs['QC_comparisons'][bestRunComp]['HOM_HET']
		run2Errs = runs['QC_comparisons'][bestRunComp]['WT_HET'] + runs['QC_comparisons'][bestRunComp]['WT_HOM'] + runs['QC_comparisons'][bestRunComp]['HET_HOM']
		if run1Errs < run2Errs:
			return run1_name
		else:
			return run2_name
	elif run1_name in avail_runs:
		return run1_name
	else:
		return run2_name

# @param run the name of the current run
# @param runs the dictionary of runs for the current sample
# @param ref_type either "normal_ref", "tumor_ref", "ref_run", or the pair name (i.e. "N-1vsT-1") for same day pairs
# @param ref_action the action to take for writing the ref run. either "write", "-", or "REF"
# @param comp_type either 'same_' or "tn_" or ""
# @param mean_error_rate for this type of comparison. used to calc zscore
# @param SD used to calc zscore
# @returns a dictionary for this run containing the runmetrics
def getRefRunMetrics(run, sample_dict, ref_type, ref_action, comp_type, mean_error_rate, SD):
	ref_run = '' # the name of the ref_run. Will be left blank if this is a same_day tumor_normal pair
	QC_comp = '' # the name of the QC_comparison. (i.e. Run1vsRun2)
	pair_metrics = {}
	if ref_type not in sample_dict:
		# this must be a same_day_pair. I passed in the N-1vsT-1 (pair name) rather than the ref type
		QC_comp = ref_type
	else:
		ref_run = sample_dict[ref_type]
		# get nn, tt, or germline QC comparison
		QC_comp = run + 'vs' + ref_run
		# check to see if hte combination is correct
		if QC_comp not in sample_dict['QC_comparisons']:
			QC_comp = ref_run + 'vs' + run
	pair_metrics[comp_type + 'pair'] = QC_comp
	# If this isn't the ref run (same day pair will also return true here), or if this is the reference run and ref_action == "write", write the metrics
	if run != ref_run or (run == ref_run and ref_action == "write"):  
		# get the error metrics for this pair
		for key in sample_dict['QC_comparisons'][QC_comp]:
			# Don't add the GT info
			if not re.search("[A-Z]", key):
				# append the comp_type for printing later
				pair_metrics[comp_type + key] = sample_dict['QC_comparisons'][QC_comp][key]
		# get the zscore and the perc_avail_bases using the current error rate and bases, and the mean error rate of the project, and Standard deviation of the project
		error_rate = pair_metrics[comp_type + 'error_rate']
		total_eligible_bases = pair_metrics[comp_type + 'total_eligible_bases']
		pair_metrics[comp_type + 'zscore'], pair_metrics[comp_type + 'perc_avail_bases'] = calcStats(error_rate, total_eligible_bases, mean_error_rate, SD)
	else:
		# this is the ref_run. set either "-" or "REF" for each of the metrics
		for metric in ['pair', 'error_count', 'error_rate', 'total_eligible_bases', 'zscore', 'perc_avail_bases']:
			pair_metrics[comp_type + metric] = ref_action	
	return pair_metrics


# @params takes a list of error Rates and finds the standardDeviation of the error rates.
# @return returns the meanErrorRate and the standardDeviation
def getStandardDeviation(errorRates):
	# Calculate the mean Error rate and the standard Deviation in order to calculate the zscore statistic for each error rate.
	total_errorRates = 0
	for rate in errorRates:
		total_errorRates += rate
	meanErrorRate = float(total_errorRates) / len(errorRates)
	#print meanErrorRate
	total_deviation = 0
	for rate in errorRates:
		total_deviation += (rate - meanErrorRate) ** 2
	#print 'total_deviation:',total_deviation, 'total_rates:', len(errorRates)
	standardDeviation = math.sqrt(float(total_deviation) / float(len(errorRates)))
	#print standardDeviation
	return meanErrorRate, standardDeviation


# Function to get error Rate stats for each run.
def calcStats(errRate, total_eligible_bases, meanErrorRate, standardDeviation):
	# the total_possible_bases is a global variable setup by main
	percAvailable = float(total_eligible_bases) / total_possible_bases
	try:
		zscore = (errRate - meanErrorRate) / float(standardDeviation)
	except ZeroDivisionError:
		zscore = 0
	return zscore, percAvailable 


# Write the master output .json file, or better, return the .json file
# @param all_samples_info the dictionary containing the QC metrics for each run of each sample in the project
# @returns a json file ready to be passed to QC_generateSheets.py to be written into the QC_Spreadsheet and a list of schema tuples containing the output format for headers and for writing the output.
def generateMasterJson(all_samples_info):
	# a list containing the name of the header, and the format. the order of the schema is the order the headers will be written 
	schema = []
	# the master json file holding all of the metrics to be written for each run of each sample.
	out_json = {}
	
	# loop through the samples
	for sample in all_samples_info:
		# get the output for each run
		if sample != "metrics":
			out_json[sample] = {}
			for run, run_data in all_samples_info[sample]['runs'].iteritems():
				out_json[sample][run] = {'sample': sample, 'run_num': run}
				out_json[sample][run] = run_data
				if all_samples_info[sample]['tumor_normal']:
					# If this sample has tumor normal runs, then add each run's same (n-n or t-t) and diff (t-n) comparison
					if run_data['run_type'] == "normal":
						# get the nn comparison add this run's error metrics based on the ref_run found 
						pair_metrics = getRefRunMetrics(run, all_samples_info[sample], 'normal_ref', 'REF', 'same_', all_samples_info['metrics']['nn_mean_error_rate'], all_samples_info['metrics']['nn_SD'])
						out_json[sample][run] = dict(out_json[sample][run].items() + pair_metrics.items())
					else:
						# get t-t comparison add this run's error metrics based on the ref_run found 
						pair_metrics = getRefRunMetrics(run, all_samples_info[sample], 'tumor_ref', 'REF', 'same_', all_samples_info['metrics']['tt_mean_error_rate'], all_samples_info['metrics']['tt_SD'])
						out_json[sample][run] = dict(out_json[sample][run].items() + pair_metrics.items())
		
					# Now write either the same_day_pair, or the best tumor_normal pair if no same_day_pair is found.
					if 'same_day_pairs' in all_samples_info[sample]:
						for same_day_pair in all_samples_info[sample]['same_day_pairs']:
							if re.search(run, same_day_pair):
								pair_metrics = getRefRunMetrics(run, all_samples_info[sample], same_day_pair, 'REF', 'tn_', all_samples_info['metrics']['tn_mean_error_rate'], all_samples_info['metrics']['tn_SD'])
								out_json[sample][run] = dict(out_json[sample][run].items() + pair_metrics.items())
		
				else:
					# This is germline, then just set each runs metrics
					pair_metrics = getRefRunMetrics(run, all_samples_info[sample], 'ref_run', 'REF', '', all_samples_info['metrics']['mean_error_rate'], all_samples_info['metrics']['SD'])
					out_json[sample][run] = dict(out_json[sample][run].items() + pair_metrics.items())
	
	return out_json


# -----------------------------------------------------------------------
# ---------------------- PROGRAM STARTS HERE ----------------------------
# -----------------------------------------------------------------------

# I put everything in functions so QC_generateSheets.py can just run the functions in this script to get the dicitonary rather than writing it back out to a .json file
# This script can either be called alone, or imported. Just use main

# If this script is called from the command line, parse the arguments, and output the .json file containing all of the stats for this project
if (__name__ == "__main__"):
	# parse the arguments
	parser = OptionParser()

	# All of the arguments are specified here.
	parser.add_option('-p', '--project_path', dest='project_path', help='REQUIRED: /path/to/the/project_dir')
	parser.add_option('-o', '--out', dest='out_json', default='qc_stats.json', help='Specify the path/to/.json file that QC_generateSheets.py will use [default: %default]')
	parser.add_option('-r', '--run_info_only', dest='run_info_only', action="store_true", default=False, help="Get only the individual run's info")
	parser.add_option('-q', '--qc_info_only', dest='qc_info_only', action="store_true", default=False, help="Get only the QC comparison info")
	parser.add_option('-b', '--bases', dest='bases', type="int", nargs=3, help="Specify 1: the total_expected_bases, 2: the total_targeted_bases, 3: the total_possible_bases")

	# Gets all of the command line arguments specified and puts them into the dictionary args
	(options, args) = parser.parse_args()

	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(8)
	
	# Wales bases:
	#total_expected_bases = 83046
	#total_targeted_bases = 84447
	#total_possible_bases = 124490
	

	if options.run_info_only:
		out_json = main_runs_only(options.project_path)
	elif options.qc_info_only:
		out_json = main_QC_only(options.project_path)
	elif options.bases != None:
		out_json = main(options.project_path, options.bases[2]) #total_possible_bases
	else:
		print "USAGE ERROR: --bases are needed for combining QC and run info"
		sys.exit(8)
	
	
	# dump the json file
	with open(options.out_json, 'w') as statsOut:
		json.dump(out_json, statsOut, sort_keys=True, indent=4)
	
	print 'Finished generating %s QC stats'%options.out_json
