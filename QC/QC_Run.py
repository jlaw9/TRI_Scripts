#! /usr/bin/env python2.7

import os
import json
import sys
from alignment_stats import Align_Stats
# write_json and runCommandLine are imported from tools
from tools import *

# class to run the individual run parts and 3x3 tables of the QC scripts.
class QC_Run:
	# nothing to do here so far...
	def __init__(self, sample_json, recalc_3x3_tables):
		self.no_errors = True
		self.sample_json = sample_json
		self.recalc_3x3_tables = recalc_3x3_tables
		self.__QCDirectory = "%s/QC"%self.sample_json['analysis']['software_directory']

	# get the separate the runs of the sample by status
	def get_runs_status(self, runs):
		passing_runs = []
		pending_runs = []
		for run in runs:
			try:
				run_json = json.load(open(run))
				if run_json['pass_fail_status'] == 'pending':
					pending_runs.append(run)
				# only the runs that pass both of these cutoffs will be used to merge
				elif run_json['pass_fail_status'] == 'pass':
					if len(runs) > 1:
						if 'pass_fail_3x3_status' not in run_json or (run_json['pass_fail_3x3_status'] == 'pending' or run_json['pass_fail_3x3_status'] == 'pass'):
							passing_runs.append(run)
					# if there is only 1 run, then it won't have a pass_fail_3x3_status. it passes.
					else:
						passing_runs.append(run)
			except ValueError:
				pass
		print "pending: %s, passing: %s"%(pending_runs, passing_runs)
		return pending_runs, passing_runs


	# runTVC_COV will only run TVC or coverage analysis if a .vcf or .cov.xls file is not found in the sample dir.
	# @param run the run for which to run TVC and coverage analysis
	def runTVC_COV(self, run, pref):
		# load the run's json file
		run_json = json.load(open(run))
	   #default is to not flag dups
		dupFlag = '--remove_dup_flags'
	
	   #see if settings want to mark dups though
		if 'mark_dups' in self.sample_json['analysis']['settings']:
		   #if it is set to true, then we change the flag
		   if self.sample_json['analysis']['settings']['mark_dups'] == 'true':
			  dupFlag = '--flag_dups'
	
	   #default is AmpliSeq for coverage analysis
		coverageAnalysisFlag = '--ampliseq'
	
	   #see if the settings say targetseq
		if 'capture_type' in self.sample_json['analysis']['settings']:
		   #if it is set to true, then we change the flag
		   if self.sample_json['analysis']['settings']['capture_type'].lower() == 'targetseq' or self.sample_json['analysis']['settings']['capture_type'].lower() == 'target_seq':
			   coverageAnalysisFlag = '--targetseq'
	
		for file in run_json['analysis']['files']:
			command = 'bash %s/runTVC_COV.sh '%self.__QCDirectory + \
					'--ptrim PTRIM.bam ' + \
					'--cleanup %s %s '%(dupFlag, coverageAnalysisFlag) + \
					'--cov %s %s '%(self.sample_json['analysis']['settings']['qc_merged_bed'], self.sample_json['analysis']['settings']['qc_unmerged_bed']) + \
					'--tvc %s %s '%(self.sample_json['analysis']['settings']['project_bed'], self.sample_json['analysis']['settings']['%stvc_json'%pref]) + \
					'--output_dir %s %s/%s '%(run_json['run_folder'], run_json['run_folder'], file)
			# run TVC and Cov analysis on this sample.
			status = runCommandLine(command)
			if status != 0:
				sys.stderr.write("%s runTVC_COV.sh had an error!!\n"%run)
				self.no_errors = False


	# getRunInfo will only be run if one of the needed parameters is not found.
	# @param run the json file of the run
	def getRunInfo(self, run, pref):
		run_json = json.load(open(run))
		# if the following metrics have already been gathered, then skip running getRunInfo.sh
		if 'run_data' in run_json and 'ts_tv' in run_json['run_data'] and 'beg_amp_cov' in run_json['run_data'] and 'num_het' in run_json['run_data'] \
				and 'amp_cov' in run_json['run_data'] \
				and ('pools_total' in run_json['run_data'] or self.sample_json['analysis']['settings']['pool_dropout'] != True):
			print "%s already has the needed metrics. Skipping getRunInfo.sh"%run
		else:
			# QC_getRunInfo.sh gets the following metrics: % amps covered at the beg and end, Ts/Tv ratio,	# Total variants,	# HET variants, 	# HOM variants
			# It also gets the metrics from the report.pdf if it is available.
			# I had to put it all on one line because python kept complaining about formatting issues.
			qcgetruninfo="bash %s/QC_getRunInfo.sh "%self.__QCDirectory + \
					"--run_dir %s "%run_json['run_folder'] + \
					"--out_dir %s/Analysis_Files/temp_files "%run_json['run_folder'] + \
					"--amp_cov_cutoff %s "%self.sample_json['analysis']['settings']['min_amplicon_coverage'] + \
					"--depth_cutoff %s "%self.sample_json['analysis']['settings']['%smin_base_coverage'%pref] + \
					"--wt_hom_cutoff %s %s "%(self.sample_json['analysis']['settings']['%swt_cutoff'%pref], self.sample_json['analysis']['settings']['%shom_cutoff'%pref])+ \
					"--beg_bed  %s "%self.sample_json['analysis']['settings']['beg_bed'] + \
					"--end_bed %s "%self.sample_json['analysis']['settings']['end_bed'] + \
					"--project_bed %s "%str(self.sample_json['analysis']['settings']['project_bed']) + \
					"--ptrim_json %s/PTRIM.bam "%run_json['run_folder'] + \
					"--software_dir %s "%self.__QCDirectory
			#if [ "$CDS_BED" != "" ]; then
			#	qcgetruninfo="$qcgetruninfo --cds_bed $CDS_BED "
			# QC_getRunInfo's will run the pool dropout script if it hasn't already been calculated
			if self.sample_json['analysis']['settings']['pool_dropout'] == True and ('run_data' not in run_json or 'pools_total' not in run_json['run_data']):
				qcgetruninfo += "--pool_dropout "
			# cleanup will be done at the end of this script
			#run the qcgetruninfo command
			status = runCommandLine(qcgetruninfo)
			if status != 0:
				sys.stderr.write("%s QC_getRunInfo.sh had an error!!\n"%run)
				self.no_errors = False
		# if the median read length was not gathered from the report PDF, or if this is a merged bam file, then calculate the median read length
		# also mean_read_length, aligned_bases, and aq20_bases
		# Align_Stats is imported from alignment_stats.py
		#bam_file = "%s/%s"%(jsonData['run_folder'], jsonData['analysis']['files'][0])
		#output_folder = "%s/Analysis_Files"%jsonData['run_folder']
		#Align_Stats(bam_file, output_folder, run)
		Align_Stats(run)
		new_run_json = json.load(open(run))
		if 'run_data' in new_run_json and 'cutoffs' in self.sample_json['analysis']['settings'] and 'exp_median_read_length' in self.sample_json['analysis']['settings']['cutoffs']:
			# Now set the perc_exp_median_read_length
			new_run_json['run_data']['perc_exp_median_read_length'] = float(new_run_json['run_data']['median_read_length']) / float(self.sample_json['analysis']['settings']['cutoffs']['exp_median_read_length'])
			#write new json file
			write_json(run, new_run_json)


	# QC two runs with each other
	# For Tumor / Normal pairs, Run1 should be the normal run, and Run2 should be the tumor run.
	# Output will be put into a dir like: sample1/QC/Run1vsRun2
	def QC_2Runs(self, sample_json, run1, run2, pref1, pref2, merged=''):
		run1_json = json.load(open(run1))
		run2_json = json.load(open(run2))
		
		# set the paths
		if 'results_qc_json' in sample_json and 'qc_folder' in sample_json:
			qc_json = sample_json['results_qc_json']
			qc_folder = sample_json['qc_folder']
		else:
			qc_json = "%s/QC/results_QC.json"%sample_json['output_folder']
			sample_json['results_qc_json'] = qc_json
			sample_json['qc_folder'] = "%s/QC"%sample_json['output_folder']
			qc_folder = sample_json['qc_folder']
	
		# load the output QC json. will be used to check if this combination has already been made.
		qc_json_data = {}
		if os.path.isfile(qc_json):
			qc_json_data = json.load(open(qc_json))

		# run1 vs run2:
		run1vsrun2 = '%svs%s'%(run1_json['run_name'], run2_json['run_name'])

		if 'chromosomes_to_analyze'+merged not in sample_json['analysis']['settings']:
			if sample_json['project'] == 'PNET':
				sample_json['analysis']['settings']['chromosomes_to_analyze'+merged] = ['all']
			else:
				sample_json['analysis']['settings']['chromosomes_to_analyze'+merged] = sample_json['analysis']['settings']['chromosomes_to_analyze']

		# only analyze the entire exome for merged runs.
		if run1_json['json_type'] != 'merged' or run2_json['json_type'] != 'merged':
			merged = ''
		
		# IDEA: If the 'all' comparison has already been made, then pull the chr combination out of it.
		# Only the runs that pass the single run QC metrics will be QC'd together.
		# QC these two runs for every chr type that is listed in chromosomes to analyze.
		for chromosome in sample_json['analysis']['settings']['chromosomes_to_analyze'+merged]:
			# now set the output_dir
			output_dir = "%s/%s%svs%s"%(qc_folder, chromosome, run1_json['run_name'], run2_json['run_name'])
			comp_type = run1_json['run_type'] + "_" + run2_json['run_type']
	
			# If this QC comparison has already been made, skip it.
			if 'QC_comparisons' in qc_json_data and chromosome in qc_json_data['QC_comparisons'] and \
					comp_type in qc_json_data['QC_comparisons'][chromosome] and run1vsrun2 in qc_json_data['QC_comparisons'][chromosome][comp_type]:
				print "%s has already run QC_2runs. Skipping."%output_dir
			else:
				# QC these two runs. QC_2Runs.sh takes the two run dirs and finds a .bam, .vcf, and .cov.xls file in the same dir as the .bam file
				qc2runs = "bash %s/QC_2Runs.sh "%self.__QCDirectory + \
				"--run_dirs %s %s "%(run1_json['run_folder'], run2_json['run_folder']) + \
				"--json_out %s "%qc_json + \
				"--output_dir %s "%output_dir + \
				"--project_bed %s "%sample_json['analysis']['settings']['project_bed'] + \
				"-a %s "%sample_json['analysis']['settings']['min_amplicon_coverage'] + \
				"-jp %s %s "%(sample_json['analysis']['settings']['%stvc_json'%pref1], sample_json['analysis']['settings']['%stvc_json'%pref2]) + \
				"-d %s %s "%(sample_json['analysis']['settings']['%smin_base_coverage'%pref1], sample_json['analysis']['settings']['%smin_base_coverage'%pref2]) + \
				"-gt %s %s %s %s "%(sample_json['analysis']['settings']['%swt_cutoff'%pref1], sample_json['analysis']['settings']['%shom_cutoff'%pref1], sample_json['analysis']['settings']['%swt_cutoff'%pref2], sample_json['analysis']['settings']['%shom_cutoff'%pref2]) + \
				"--software_dir %s "%self.__QCDirectory
				#"--cleanup " # The main cleanup will be done at the end of this script because the PTRIM.bam is needed for QC_getRunInfo.sh, and the chr_subset is needed for each run comparison.

				# subset this specified chromosome
				if chromosome == "718":
					qc2runs += "--subset_bed %s "%sample_json['analysis']['settings']['subset_bed']
					qc2runs += "--all_vcfs_dir %s/all%svs%s "%(qc_folder, run1_json['run_name'], run2_json['run_name'])
				elif chromosome != "all":
					qc2runs += "--subset_chr %s "%chromosome
	
				# if the recalc option is specified, we might be able to pass in the total eligible and possible bases because those shouldn't change (unless the amplicon cutoff is changed from 30x).
				if self.recalc_3x3_tables and 'old_GTs' in qc_json_data and chromosome in qc_json_data['old_GTs'] and comp_type in qc_json_data['old_GTs'][chromosome] and \
					run1vsrun2 in qc_json_data['old_GTs'][chromosome][comp_type]:
					# get the bases and add them to the command
					# recalculate the total_possible_bases
					total_possible_bases = int(1 / (qc_json_data['old_GTs'][chromosome][comp_type][run1vsrun2]['perc_avail_bases'] / qc_json_data['old_GTs'][chromosome][comp_type][run1vsrun2]['total_eligible_bases']))
					qc2runs += "--bases %s %s "%(qc_json_data['old_GTs'][chromosome][comp_type][run1vsrun2]['total_eligible_bases'], total_possible_bases)
		
				#run the qc2runs command
				sys.stdout.write("Running: %s  at: %s\n"%(qc2runs, getTimestamp()))
				status = runCommandLine(qc2runs)
				if status != 0:
					sys.stderr.write("%s vs %s QC_2Runs.sh had an error!!\n"%(run1, run2))
					self.no_errors = False
		return sample_json, qc_json

	
	# @param run_json Update the pass/fail status of this run according to the cutoffs specified here.
	def update_run_status(self, run, num_runs):
		run_json = json.load(open(run))
		status = 'pass'
		# Update differently for runs vs merged bam files
		if 'json_type' not in run_json:
			print "ERROR: %s/%s unable to check the pass/fail cutoffs because the 'json_type' of 'run'/'merged' is not found in the json file"%(run_json['sample'], run_json['run_name'])
			return
		if 'cutoffs' not in self.sample_json['analysis']['settings']:
			print "Unable to update run status based on cutoffs because cutoffs are unavailable"
		elif run_json['json_type'] == 'run':
			# first check the % expected median read length. If this is an ffpe sample, then we don't fail the read length cutoff.
			if 'perc_exp_median_read_length' in self.sample_json['analysis']['settings']['cutoffs'] \
					and run_json['run_data']['perc_exp_median_read_length'] < self.sample_json['analysis']['settings']['cutoffs']['perc_exp_median_read_length'] \
					and not ('ffpe' in self.sample_json['analysis']['settings'] and self.sample_json['analysis']['settings']['ffpe'] == True):
				print "%s/%s fialed the 'perc_exp_median_read_length' flag. %.2f < %.2f"%(run_json['sample'], run_json['run_name'], run_json['run_data']['perc_exp_median_read_length'], self.sample_json['analysis']['settings']['cutoffs']['perc_exp_median_read_length'])
				status = 'FAIL'
			# check the coverage at the +10 and -10 positions
			if 'begin_end_amp_cov' in self.sample_json['analysis']['settings']['cutoffs'] \
					and run_json['run_data']['begin_amp_cov'] < self.sample_json['analysis']['settings']['cutoffs']['begin_end_amp_cov']:
				print "%s/%s fialed the 'begin_amp_cov' flag. %.2f < %.2f"%(run_json['sample'], run_json['run_name'], run_json['run_data']['begin_amp_cov'], self.sample_json['analysis']['settings']['cutoffs']['begin_end_amp_cov'])
				status = 'FAIL'
			if 'begin_end_amp_cov' in self.sample_json['analysis']['settings']['cutoffs'] \
					and run_json['run_data']['end_amp_cov'] < self.sample_json['analysis']['settings']['cutoffs']['begin_end_amp_cov']:
				print "%s/%s fialed the 'end_amp_cov' flag. %.2f < %.2f"%(run_json['sample'], run_json['run_name'], run_json['run_data']['end_amp_cov'], self.sample_json['analysis']['settings']['cutoffs']['begin_end_amp_cov'])
				status = 'FAIL'
			# checck the coverage
			if 'run_amp_cov' in self.sample_json['analysis']['settings']['cutoffs'] \
					and run_json['run_data']['amp_cov'] < self.sample_json['analysis']['settings']['cutoffs']['run_amp_cov']:
				print "%s/%s fialed the 'run_amp_cov' flag. %.2f < %.2f"%(run_json['sample'], run_json['run_name'], run_json['run_data']['amp_cov'], self.sample_json['analysis']['settings']['cutoffs']['run_amp_cov'])
				status = 'FAIL'
			# check the pool coverage
			#if self.sample_json['analysis']['settings']['pool_dropout'] and 'pools_between_10_and_50' in run_json['run_data'] and 'pools_less_than_10' in run_json['run_data']:
			#	# currently the only options are 10, 50, and 75 so just stick with this.
			#	if run_json['run_data']['pools_between_10_and_50'] > 0 or run_json['run_data']['pools_less_than_10'] > 0:
			#		print "%s/%s has pools with < 50%% covergae"%(run_json['sample'], run_json['run_name'])
			#		status = 'FAIL'
		# Update the merged bam status differently 
		elif run_json['json_type'] == 'merged':
			if 'merged_amp_cov' in self.sample_json['analysis']['settings']['cutoffs'] and run_json['run_data']['amp_cov'] < self.sample_json['analysis']['settings']['cutoffs']['merged_amp_cov']:
				print "%s/%s fialed the 'merged_amp_cov' flag. %.2f < %.2f"%(run_json['sample'], run_json['run_name'], run_json['run_data']['amp_cov'], self.sample_json['analysis']['settings']['cutoffs']['merged_amp_cov'])
				status = 'REQUEUE'
		# set the pass fail status for this run
		run_json['pass_fail_status'] = status
		# write this run's updated status to the json file
		write_json(run, run_json)

	# Update run status based on the 3x3 table error rates.
	#- If any of the 3x3 table error rates fail, then the runs for this sample won't be merged until the 3x3 tables are manually reviewed to determine the faulty run, or if the sample needs to be resequenced.
	def update_3x3_runs_status(self, sample_json, run1, run2, qc_json):
		perc_avail_bases = 0
		qc_data = json.load(open(qc_json))
		# set the default sample_status to be 'pending_merge' if it has not yet been set.
		# Otherwise it will be set to 'pending_3x3_review' if any of the 3x3 table comparisons failed the 'error_rate' cutoff
		if sample_json['sample_status'] != 'pending_3x3_review':
			sample_json['sample_status'] = 'pending_merge'
		# load the jsons
		run1_json = json.load(open(run1))
		run2_json = json.load(open(run2))
		# set the default status to pass for each run if it has not yet been set so as to not overwrite a run that is already set to 'fail' to 'pass' in the case of tumor/normal comparisons
		if 'pass_fail_3x3_status' not in run1_json or run1_json['pass_fail_3x3_status'] == 'pending':
			run1_json['pass_fail_3x3_status'] = "pass"
		if 'pass_fail_3x3_status' not in run2_json or run2_json['pass_fail_3x3_status'] == 'pending':
			run2_json['pass_fail_3x3_status'] = "pass"
		if 'cutoffs' in sample_json['analysis']['settings'] and 'error_rate' in sample_json['analysis']['settings']['cutoffs']:
			# run1vsrun2:
			comp_type = run1_json['run_type'] + "_" + run2_json['run_type']
			run1vsrun2 = run1_json['run_name'] + "vs" + run2_json['run_name']
			for chr in qc_data['QC_comparisons']:
				# Check to see if this comparison passes the 'error_rate' cutoff
				#if 'error_rate' not in sample_json['analysis']['settings']['cutoffs']:
					# Use this a default error_rate cutoff of 3.0e-05
					#sample_json['analysis']['settings']['cutoffs']['error_rate'] = 0.00003
				if comp_type in qc_data['QC_comparisons'][chr] and run1vsrun2 in qc_data['QC_comparisons'][chr][comp_type] \
						and qc_data['QC_comparisons'][chr][comp_type][run1vsrun2]['error_rate'] > sample_json['analysis']['settings']['cutoffs']['error_rate']:
					print "%s%svs%s Failed with an error rate too high of %.3e"%(chr, run1_json['run_name'], run2_json['run_name'], qc_data['QC_comparisons'][chr][comp_type][run1vsrun2]['error_rate'])
					run1_json['pass_fail_3x3_status'] = 'FAIL'
					run2_json['pass_fail_3x3_status'] = 'FAIL'
					sample_json['sample_status'] = 'pending_3x3_review'

				# add this variable for the final merged normal/tumor comparison
				if comp_type in qc_data['QC_comparisons'][chr] and run1vsrun2 in qc_data['QC_comparisons'][chr][comp_type]:
					perc_avail_bases = qc_data['QC_comparisons'][chr][comp_type][run1vsrun2]['perc_avail_bases']
		else:
			# if teh cutoffs are not available, then don't merge the runs yet.
			sample_json['sample_status'] = 'pending_3x3_review'
		# write the runs updated status to the json file
		write_json(run1, run1_json)
		write_json(run2, run2_json)

		#return the perc_avail_bases for the merged normal/tumor comparison
		return sample_json, perc_avail_bases

