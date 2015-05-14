#! /usr/bin/env python

# Goal:  Merge and run coverage analysis on the two Samples generated.
# Output:  A mered bam file, and coverage analysis on the merged bam file.

from optparse import OptionParser
import os
import os.path
import sys
import re
import datetime
import json
from QC_Run import QC_Run
from tools import *

class Merger:
	# @param bams_to_merge a list of the bam files to merge together
	# @param merged_dir the directory in which to place the merged bam file
	# @param sample_name the name of the sample. Used for the SM tag
	# @param cleanup Flag to delete the temporary files or not. Default: false
	def __init__(self, sample_json, recalc_3x3_tables):
		if sample_json:
			self.sample_json = json.load(open(sample_json))
		self.merge_dir = ''
		self.bams_to_merge = []
		self.runs_to_merge = []
		self.QC_Run = QC_Run(self.sample_json, recalc_3x3_tables)

	
	# merge the following runs
	def merge(self):
		# this could be just a temporary fix
		if os.path.isfile(self.path_to_merged_bam):
			print "%s already exists. Not making it again."%self.path_to_merged_bam
		else:
			print "Sample %s is merging the following runs:  %s"%(self.sample_name, self.bams_to_merge)
		
			merge_command = "java -jar /opt/picard/picard-tools-current/MergeSamFiles.jar "
		
			# Add each run's bam file to mergeJob.sh
			for bam in self.bams_to_merge:
				if not os.path.isfile(bam) or bam[-4:] != ".bam":
					print "ERROR: the bam file '%s' does not exist!"%bam
					sys.exit(4)
				merge_command += "INPUT=%s "%bam

			# make sure the merged_dir exists, or make it.
			runCommandLine("mkdir -p %s"%self.merged_dir)
			#if not os.path.isdir(merged_dir):
				#print "ERROR: the output dir '%s' does not exist!"%bam
				#sys.exit(4)
				
			# Now set the output file, and then run the merge command
			merge_command +=  "	OUTPUT=%s/merged_badHeader.bam "%self.merged_dir
		
			if runCommandLine(merge_command) != 0:
				print "ERROR: %s something went wrong with merging!"%self.sample_name
				sys.exit(1)
		
			#echo "fixing header for %s/merged_badHeader.bam"
			correct_header_command = "samtools view -H %s/merged_badHeader.bam > %s/merged.header.sam "%(self.merged_dir, self.merged_dir)
			if runCommandLine(correct_header_command) != 0:
				print "ERROR: samtools view -H failed!"
				sys.exit(1)
	
			# A better way would be to check to see if the SM tags already match. Then we would be able to use ionstats and such.	
			SM_check_command = "grep -Eo 'SM:[a-zA-Z0-9_&/-]*"
			# NEED TO TEST THIS COMMAND. Is there anything that comes before the next : that is important?
			# Change the SM: tag so that it matches for every run merged. (There should be one SM tag for each run merged)
			# This was the old command. We will keep using this, and then if there are problems, we can manually correct them.
			sed_command = 'sed "s/SM:[a-zA-Z0-9_&/-]*/SM:%s/" %s/merged.header.sam > %s/merged.headerCorrected.sam'%(self.sample_name, self.merged_dir, self.merged_dir)
			# this updated command will change the SM tag to match everything up to the next : after the SM tag.
			# this command deletes the KS: tag!! not good! I don't know why but some headers are tab delimited, and some are not it seems.
			#sed_command = 'sed -E "s/SM:[^:]*:/SM:%s:/" %s/merged.header.sam > %s/merged.headerCorrected.sam'%(self.sample_name, self.merged_dir, self.merged_dir)
			if runCommandLine(sed_command) != 0:
				print "ERROR: sed command failed!"
				sys.exit(1)
		
			# write the new header to merged.bam
			reheader_command = "samtools reheader %s/merged.headerCorrected.sam %s/merged_badHeader.bam > %s "%(self.merged_dir, self.merged_dir, self.path_to_merged_bam)
			if runCommandLine(reheader_command) != 0:
				print "ERROR: sed command failed!"
				sys.exit(1)

			# set some extra variables for the JSON file.
			self.merged_json = "%s/merged.json"%self.merged_dir

			# if there is already an index file from a previous merge try, delete it.
			if os.path.isfile(self.path_to_merged_bam + ".bai"):
				os.remove(self.path_to_merged_bam + ".bai")

			# IF specified, cleanup the temporary files
			#if self.cleanup:
			# Need to cleanup here inorder for TVC to work. there can only be one bam file in the merged dir.
			os.remove("%s/merged_badHeader.bam"%self.merged_dir)
			os.remove("%s/merged.headerCorrected.sam"%self.merged_dir)
			os.remove("%s/merged.header.sam"%self.merged_dir)
		
			print "%s finished merging "%self.merged_dir


	# Update the final merged run status
	def update_merged_run_status(self, run, merged_perc_avail_bases=0):
		pass_fail_merged_status = 'pass'
		run_json = json.load(open(run))
		if run_json['run_type'] == 'germline':
			merged_perc_aval_bases = run_json['run_data']['amp_cov']
		#print merged_perc_avail_bases, self.sample_json['analysis']['settings']['cutoffs']['merged_amp_cov']
		# check to see if >90% of the bases are shared between the tumor normal comparison 
		if 'merged_amp_cov' in self.sample_json['analysis']['settings']['cutoffs'] and merged_perc_avail_bases != '':
			if merged_perc_avail_bases < self.sample_json['analysis']['settings']['cutoffs']['merged_amp_cov']:
				pass_fail_merged_status = 'REQUEUE'
			# write the final statuses here
			run_json['pass_fail_merged_status'] = pass_fail_merged_status
			run_json['merged_perc_avail_bases'] = merged_perc_avail_bases
			write_json(run, run_json)


	# @param runs the runs of a sample 
	# @param run_name either '', 'Normal/' or 'Tumor/'
	# @param pref the prefix of this type of merge. either 'normal_' 'tumor_' or ''
	# @returns a list of the passing bam files to merge, and the path to the merged dir.
	def check_merge(self, runs, run_name='', pref=''):
		# vars to return
		merge = False
		self.bams_to_merge = []
		self.runs_to_merge = []
		# Use this count so that we won't have to write over past merges if there are multiple merges.
		if 'merged_%scount'%pref not in self.sample_json:
			self.sample_json['merged_%scount'%pref] = 0
		# first check to see if all of the runs pass. 
		# Get all of the passing bam files for this sample.
		pending_runs, passing_runs = self.QC_Run.get_runs_status(runs)
		if len(pending_runs) != 0:
			print "Not merging. After QC_runs, runs should either be 'pass' or 'fail', not 'pending'. Pending runs: ", pending_runs
		elif len(passing_runs) < 1:
			# if none of the runs are passing, then don't do anything.
			pass
		elif self.sample_json['sample_status'] != "pending_merge" and self.sample_json['sample_status'] != "merged":
			# If any runs of the sample are not ready to be merged either because of 3x3 table error rate questions or other reasons, don't merge this sample.
			print "%s the 'sample_status' is '%s'. Needs to be 'pending_merge' to merge the runs."%(self.sample_json['sample_name'], self.sample_json['sample_status'])
		elif self.sample_json['sample_status'] == 'pending_merge':
			# Merge these runs.
			# First get the passing bams from the passing runs.
			for run in passing_runs:
				run_json = json.load(open(run))
				self.bams_to_merge.append("%s/%s"%(run_json['run_folder'], run_json['analysis']['files'][0]))
				self.runs_to_merge.append(run_json['run_name'])
				# sort the run names
				self.runs_to_merge.sort()
	
			# If this sample has already been merged:  If the runs to generate the merged bam don't match the current list: 
				# then delete the last created bam file and merge these runs
			# else don't remerge these files
			if len(self.bams_to_merge) == 1:
				# There is only one run, so don't merge it. Set the "final_%sjson"%pref flag to show what the final run is
				self.sample_json["final_%sjson"%pref] = run
			# use the 'merged_json' flag rather than the 'final_json' flag because 'final_json' can be set by a single non-merged run.
			elif 'merged_%sjson'%pref in self.sample_json and os.path.isfile(self.sample_json['merged_%sjson'%pref]):
				merged_json_data = json.load(open(self.sample_json['merged_%sjson'%pref]))
				# If the runs used to generate the current merged.bam file dont match the current bams_to_merge, then merge them. Otherwise don't
				if merged_json_data['json_type'] == 'merged' and set(self.bams_to_merge) != set(merged_json_data['bams_used_to_merge']):
					# in order to manage space, delete the last merged folder that was created.
					if self.sample_json['analysis']['settings']['cleanup'] == True: 
						# IDEA delete the entire folder? Or just the bam file?
						merged_bam = "%s/%s"%(merged_json_data['run_folder'], merged_json_data['analysis']['files'][0])
						print "	Deleting the old merged bam file: %s"%merged_bam
						os.remove(merged_bam)
					# Add one to the merged_count
					self.sample_json['merged_%scount'%pref] += 1
					# set new path to the merged_json
					self.merged_dir = "%s/%sMerged_%d"%(self.sample_json['sample_folder'], run_name, self.sample_json['merged_%scount'%pref])
					merge = True
				else:
					# Don't merge these runs because they've already been merged.
					print "%s the runs: '%s' have already been merged"%(self.sample_json['sample_name'], self.bams_to_merge)
			else:
				# Merge these runs
				self.merged_dir = "%s/%sMerged"%(self.sample_json['sample_folder'], run_name)
				# Add one to the merged_count
				self.sample_json['merged_%scount'%pref] += 1
				merge = True
		return merge


	# merge the runs of a sample
	# @param runs the bam files to merge
	# @param merged_dir the ouptut_dir in which to place the merged bam file
	# @param pref the prefix (either '', 'normal_', or 'tumor')
	# @param run_type either germline, normal, or tumor.
	# @param run_name either Merged, Normal_Merged or Tumor_Merged. Used for the titles of the 3x3 tables.
	def merge_runs(self, run_type, run_name='', pref=''):
		# if the file already exists, then merging must have finished, and don't merge again.
		self.merged_json = "%s/merged.json"%self.merged_dir
		if os.path.isfile(self.merged_json):
			print "%s already exists so not merging the bam files again"%self.merged_json
		else:
			self.sample_name = self.sample_json['sample_name']

			# get today's date to format the mreged.bam file name
			curr_date = datetime.date.today()
			# the name follows this format: A_227_Tumor_Merged_02152015
			run_name = "%s_%sMerged_%02d%02d%s"%(self.sample_name, run_name, curr_date.month, curr_date.day, curr_date.year)
			merged_bam = "%s.bam"%(run_name)
			self.path_to_merged_bam = "%s/%s"%(self.merged_dir, merged_bam)

			self.merge()
	
			# now set the json files
			# create the merged_bam's json file here so that the merger.py script can run on its own if necessary.
			merged_json = {
					'analysis': {
						'files': [merged_bam] 
						},
					'bams_used_to_merge':self.bams_to_merge, 
					'sample_name': self.sample_name, 
					'merged_bam': self.path_to_merged_bam,
					'json_file': self.merged_json,
					"json_type": "merged",
					"pass_fail_status": "pending", 
					"project": self.sample_json['project'], 
					"run_folder": self.merged_dir, 
					"run_name": run_name,
					"run_num": self.sample_json['merged_%scount'%pref],
					"run_type": run_type, 
					"runs_used_to_merge": ', '.join(self.runs_to_merge),
					"sample": self.sample_json['sample_name'], 
					"sample_folder": self.sample_json['sample_folder'],
					"sample_json": self.sample_json['json_file']
				}
			#write new json file
			write_json(self.merged_json, merged_json)
				
		# QC the merged run.
		self.QC_Run.runTVC_COV(self.merged_json, pref)
		self.QC_Run.getRunInfo(self.merged_json, pref)
		
		# Update the merge pass/fail status based on the metrics gathered by QC_getRunInfo.sh
		self.QC_Run.update_run_status(self.merged_json, 1)

		# Also store the path to this merged bam file in the sample's json file. Not really necessary, but it seems like a good idea.
		#if 'merged' not in self.sample_json['analysis']['files']:
		#	self.sample_json['analysis']['files']['merged'] = {}
		#self.sample_json['analysis']['files']['merged']['%sbam'%pref] = merger.path_to_merged_bam	
		# store the path to this merged bam folder in the sample's json file.
		#self.sample_json['merged_%sjson'%pref] = merged_dir

		# If the merge_json passes the cutoffs, set it as the final_json
		merge_json = json.load(open(self.merged_json))
		# add the path to this merge even if it doesn't pass
		self.sample_json["merged_%sjson"%pref] = self.merged_json
		if merge_json['pass_fail_status'] == 'pass':
			# Add a path to the final merged_json
			self.sample_json["final_%sjson"%pref] = self.merged_json

		# write the modified sample_json file
		write_json(self.sample_json['json_file'], self.sample_json)


# If we need this script to run on its own, update it when it is needed
#if __name__ == '__main__':
#
#	# set up the option parser
#	parser = OptionParser()
#	
#	# add the options to parse
#	parser.add_option('-j', '--json', dest='json', help='The samples json file. Will be used to get the passing bams.')
#	parser.add_option('-o', '--merged_dir', dest='output', help='The output file. If no output file is specified, output will be written to the screen')
#	parser.add_option('-s', '--sample_name', dest='sample', help='The name of the sample. Will be used to fix the SM tag of the merged BAM file')
#	parser.add_option('-b', '--bams', dest='bams', action='append', help='Use a -b for for each bam to include in merging')
#	parser.add_option('-c', '--cleanup', dest='cleanup', action='store_true', help='option to cleanup the temporary files used in merging and such.')
#
#	(options, args) = parser.parse_args()
#
#	if options.json and (not options.output and not options.sample and not options.bams):
#		Merger(options.json)
#	# I don't have time to implement these other options yet... 
#	#elif not options.json and (options.output and options.sample and options.bams):
##		merger = Merger()
##		merger.merge()
##		Merger(options.bams, options.output, options.sample)
#	else:
#		print "USAGE_ERROR: -j or (-o, -s and -b) are required. If the json file is provided, do not provide the other options. If the other options are provided, do not provide a json file."
#		print "only -j is implemented so far..."
#		parser.print_help()
#		sys.exit(1)
#
