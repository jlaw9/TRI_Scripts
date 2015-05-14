#! /usr/bin/env python2.7

# Goal: push a sample through the QC process.

import traceback
from optparse import OptionParser
import os
import sys
import re
import json
import glob
import shutil
from QC_Run import QC_Run
from merger import Merger
from cleanup import Cleanup
from alignment_stats import Align_Stats
# write_json and runCommandLine are imported from tools
from tools import *

class QC_Sample:
	def __init__(self, options):
		self.__QCDirectory = "/rawdata/legos/scripts/QC"
		self.no_errors = True
		self.cleanup_sample = Cleanup()
		self.options = options
		self.sample_json = json.load(open(options.json))
		self.qc_run = QC_Run(self.sample_json, options.recalc_3x3_tables)
	
	# will find all of the runs in a sample and QC them with each other
	def QC_merge_runs(self):
		# if this is a germline sample, QC all of the normal runs with each other.
		if self.sample_json['sample_type'] == 'germline':
			# Use the sample_status here to not re-run the QC and to not overwrite run status. The 'sample_status' should be reset to 'pushed' when new runs are pushed..
			#if self.sample_json['sample_status'] != 'pending_merge' and self.sample_json['sample_status'] != 'pending_3x3_review' and self.sample_json['sample_status'] != 'merged':
			# if the user specified the '--pass_fail' option, then run this part still
			if self.sample_json['sample_status'] == 'pushed' or self.options.pass_fail or self.options.qc_all:
				# QC the normal runs with each other
				self.QC_runs(self.sample_json['runs'])
				# write the sample json file
				write_json(self.sample_json['json_file'], self.sample_json)
	
			# what if there is only one run that passes all of the metrics? It should be marked as the 'final_json' and have the 'pass_fail_merged' flag marked as pass.
			# make the merger
			merger = Merger(self.sample_json['json_file'], self.options.recalc_3x3_tables)
			# Check to see if the normal runs are ready to be merged.
			merge = merger.check_merge(self.sample_json['runs'])
			if merge == True:
				# merge the normal and/or tumor runs. Will only merge the passing runs with each other.
				merger.merge_runs('germline')

				# load the sample json file because merger edited it.
				self.sample_json = json.load(open(self.sample_json['json_file']))

				# update the merged run status
				merger.update_merged_run_status(self.sample_json['merged_json'])
	
				if json.load(open(self.sample_json['merged_json']))['pass_fail_merged_status'] == 'pass':
					# Set the sample_status
					self.sample_json['sample_status'] = 'merged_pass'
					# cleanup the individual run bam files
					self.cleanup_sample.cleanup_runs(self.sample_json['runs'], self.sample_json['analysis']['settings']['cleanup'], self.no_errors)
					# Cleanup the merged dir 
					self.cleanup_sample.cleanup_runs([self.sample_json['merged_json']], self.sample_json['analysis']['settings']['cleanup'], self.no_errors)
				else:
					self.sample_json['sample_status'] = 'awaiting_more_sequencing'
	

		# if this is a tumor_normal sample, find the normal and tumor runs, and then QC them with each other.
		elif self.sample_json['sample_type'] == 'tumor_normal':
			# Separate the runs into tumor and normal lists
			normal_runs, tumor_runs = self.getTumor_Normal()
	
			if self.sample_json['analysis']['settings']['type'] == 'all_tumor_normal':
				# Use the sample_status here to not re-run the QC and to not overwrite run status. The 'sample_status' should be reset to 'pushed' when new runs are pushed..
				#if self.sample_json['sample_status'] != 'pending_merge' and self.sample_json['sample_status'] != 'pending_3x3_review' and self.sample_json['sample_status'] != 'merged':
				# if the user specified the '--pass_fail' option, then run this part still
				if self.sample_json['sample_status'] == 'pushed' or self.options.pass_fail or self.options.qc_all:
					# QC the normal or tumor runs with each other
					self.QC_runs(normal_runs, 'normal_')
					self.QC_runs(tumor_runs, 'tumor_')
					# now QC the tumor and normal runs together.
					self.QC_normal_tumor_runs(normal_runs, tumor_runs)
					# make the excel spreadsheet containing the data and copy it back to the proton
					#self._make_xlsx()
					# write the sample json file
					write_json(self.sample_json['json_file'], self.sample_json)
	
				# make the merger
				merger = Merger(self.sample_json['json_file'], self.options.recalc_3x3_tables)
				# Check to see if the normal runs are ready to be merged.
				merge_normal = merger.check_merge(normal_runs, 'Normal/', 'normal_')
				if merge_normal == True:
					# merge the normal and/or tumor runs. Will only merge the passing runs with each other.
					merger.merge_runs('normal', 'Normal_', 'normal_')
	
				# Check to see if the tumor runs are ready to be merged.
				merge_tumor = merger.check_merge(tumor_runs, 'Tumor/', 'tumor_')
				if merge_tumor == True:
					merger.merge_runs('tumor', 'Tumor_', 'tumor_')

				# load the sample json file because merger edited it.
				self.sample_json = json.load(open(self.sample_json['json_file']))
	
				# If any runs were merged, QC them. If there are only 1 normal and tumor run, they won't be QCd again. 
				#if normal_merge_dir != '' or tumor_merge_dir != '' or (len(normal_passing_bams) == 1 and len(tumor_passing_bams) == 1):	
				# now QC the tumor and normal merged bams together if both normal and tumor runs are ready.
					# To only QC all for the actual merged runs (PNET), change the 'final' part to 'merged'.
				if (merge_normal or merge_tumor) and ('final_normal_json' in self.sample_json and 'final_tumor_json' in self.sample_json):
					self.sample_json, qc_json = self.qc_run.QC_2Runs(self.sample_json, self.sample_json['final_normal_json'], self.sample_json['final_tumor_json'], 'normal_', 'tumor_', '_merged')
					self.sample_json, merged_perc_avail_bases = self.qc_run.update_3x3_runs_status(self.sample_json, self.sample_json['final_normal_json'], self.sample_json['final_tumor_json'], qc_json)
					# update the merged run status 
					merger.update_merged_run_status(self.sample_json['final_normal_json'], merged_perc_avail_bases)
					merger.update_merged_run_status(self.sample_json['final_tumor_json'], merged_perc_avail_bases)

					# cleanup the individual run bam files
					if merged_perc_avail_bases > .9:
						final_qc_dir = "%s/%svs%s"%(self.sample_json['qc_folder'], json.load(open(self.sample_json['final_normal_json']))['run_name'], json.load(open(self.sample_json['final_tumor_json']))['run_name'])
						# annotate the final somatic variants
						command = "~/TRI_Scripts/Somatic_Variants/somatic_variants.sh %s %s"%(final_qc_dir, self.sample_json['sample_name'])
						if runCommandLine(command) != 0:
							sys.stderr.write("ERROR: somatic annotation failed!\n")
	
						# Cleanup the PTRIM.bam and chr bam files after all of the QC is done.
						# are there any other files to clean up?
						self.cleanup_sample.cleanup_runs(self.sample_json['runs'], self.sample_json['analysis']['settings']['cleanup'], self.no_errors)
						#self.cleanup_sample.delete_runs(runs, self.sample_json['analysis']['settings']['cleanup'], self.no_errors)
	
						# Cleanup after the merging QC is done.
						self.cleanup_sample.cleanup_runs([self.sample_json['final_normal_json'], self.sample_json['final_tumor_json']], self.sample_json['analysis']['settings']['cleanup'], self.no_errors)
	
						# Set the sample_status
						self.sample_json['sample_status'] = 'merged_pass'
					else:
						self.sample_json['sample_status'] = 'awaiting_more_sequencing'

		# print the final status
		if self.no_errors == False or self.qc_run.no_errors == False:
			sys.stderr.write("%s finished with errors. See %s/sge.log for more details"%(self.sample_json['sample_name'], self.sample_json['output_folder']))
			self.sample_json['sample_status'] == 'failed'
			write_json(self.sample_json['json_file'], self.sample_json)
			sys.exit(1)
		else:
			print "%s finished with no errors!"%(self.sample_json['sample_name'])

		# write the sample json file
		write_json(self.sample_json['json_file'], self.sample_json)

		# make the excel spreadsheet containing the data and copy it back to the proton
		self._make_xlsx()
				

	# Separate the runs into tumor and normal lists
	def getTumor_Normal(self):
		normal_runs = []  
		tumor_runs = []
		for run in self.sample_json['runs']:
			run_json = json.load(open(run))
			# temp fix for runs that have old JSON files (i.e. SEGA)
			if 'run_type' not in run_json or 'run_num' not in run_json:
				if re.search('N-', run):
					run_json['run_type'] = 'normal'
				else:
					run_json['run_type'] = 'tumor'
				run_json['pass_fail_status'] = 'pending'
				run_json['json_type'] = 'run'
				run_json['json_file'] = run
				run_json['run_name'] = run_json['name']
				run_json['run_num'] = run_json['run_name'][-1]
				run_json['sample_name'] = run_json['sample']
				if re.search('-', run):
					run_json['run_folder'] = '/'.join(run.split('/')[:-1])
					run_json['sample_folder'] = os.path.abspath('/'.join(run.split('/')[:-1]) + "/../..")
				write_json(run, run_json)
				# temp fix over
			if 'analysis' not in run_json:
				bam = glob.glob("%s/*.bam"%run_json['run_folder'])[0].split('/')[-1]
				run_json['analysis'] = {'files': [bam]}
				write_json(run, run_json)
			if run_json['run_type'] == 'normal':
				normal_runs.append(run)
			elif run_json['run_type'] == 'tumor':
				tumor_runs.append(run)
			else:
				print "ERROR run type is not normal or tumor."
		return normal_runs, tumor_runs


	# QC the normal runs with each other 
	def QC_runs(self, runs, pref=''):
		# first run TVC_CV and get the Run info to prepare for QC2Runs
		for run in runs:
			run_json = json.load(open(run))
			# only run these if this run has a status of pending.
			# This way the pass_fail_status can be manually overwritten.
			if run_json['pass_fail_status'] == "pending" or self.options.pass_fail:
				self.qc_run.runTVC_COV(run, pref)
				self.qc_run.getRunInfo(run, pref)
				# Update the run status based on the metrics gathered by QC_getRunInfo.sh
				self.qc_run.update_run_status(run, len(runs))
		# if there is only one run for this sample, then set the status to 'pending_merge' so that the only run will be set as the 'final_json'
		passing_runs = self.qc_run.get_runs_status(runs)[0]
		if len(passing_runs) == 1:
			self.sample_json['sample_stats'] == 'pending_merge'
		else:
			for run1 in runs:
				run1_json = json.load(open(run1))
				for run2 in runs:
					run2_json = json.load(open(run2))
					# check to see if these two runs should be QC'd together. Only QC the runs that pass the single run QC metrics.
					if int(run1_json['run_num']) < int(run2_json['run_num']) and ((run1_json['pass_fail_status'] == 'pass' and run2_json['pass_fail_status'] == 'pass') or self.options.qc_all): 
						self.sample_json, qc_json = self.qc_run.QC_2Runs(self.sample_json, run1, run2, pref, pref)
						self.sample_json, perc_avail_bases = self.qc_run.update_3x3_runs_status(self.sample_json, run1, run2, qc_json)


	# now QC the tumor and normal runs together.
	def QC_normal_tumor_runs(self, normal_runs, tumor_runs):
		for normal_run in normal_runs:
			for tumor_run in tumor_runs:
				normal_json = json.load(open(normal_run))
				tumor_json = json.load(open(tumor_run))
				# Only QC the runs that pass the single run QC metrics.
				if (normal_json['pass_fail_status'] == 'pass' and tumor_json['pass_fail_status'] == 'pass') or self.options.qc_all: 
					self.sample_json, qc_json = self.qc_run.QC_2Runs(self.sample_json, normal_run, tumor_run, 'normal_', 'tumor_')
					self.sample_json, perc_avail_bases = self.qc_run.update_3x3_runs_status(self.sample_json, normal_run, tumor_run, qc_json)

	# make the xlsx file to be copied back to the proton
	def _make_xlsx(self):
		xlsx_file = '%s/%s_QC.xlsx'%(self.sample_json['qc_folder'], self.sample_json['sample_name'])
		
		make_xlsx_command = "python2.7 %s/QC_generateSheets.py "%self.__QCDirectory + \
			"--sample_path %s "%self.sample_json['sample_folder'] + \
			"--sheet_per_sample " + \
			"--out %s "%xlsx_file + \
			"--ex_json %s "%(self.sample_json['json_file'])
	
		status  = runCommandLine(make_xlsx_command)
		if status != 0:
			print "unable to generate the excel file"
		else:
			print "Generated the QC spreadsheet successfully!"
			# t would be really really cool if I could send them an email with the xlsx file!!
			if self.options.email and 'emails' in self.sample_json:
				# TEMP add my email automatically
				if 'jlaw@childhooddiseases.org' not in self.sample_json['emails']:
					self.sample_json['emails'].append('jlaw@childhooddiseases.org')
				for email in self.sample_json['emails']:
					# this command will email the status of the sample, and attach the excel spreadsheet and somatic variants if it is found.
					# TODO add germline project variants as well.
					somatic_variants = "/home/ionadmin/jeff/Lung_Somatic/%s_somatic.vcf"%self.sample_json['sample_name']
					if os.path.isfile(somatic_variants):
						email_command = '\tprintf "%s finished with a status of %s. \\n`grep sample_status *.json`\\n" | (cat -; uuencode %s %s; uuencode %s %s_%s_somatic.tsv) | ssmtp -vvv %s >/dev/null 2>&1\n' % (self.sample_json['sample_name'], "pass", xlsx_file, xlsx_file.split('/')[-1], somatic_variants, self.sample_json['project'], self.sample_json['sample_name'], email)
					else:
						email_command = '\tprintf "%s finished with a status of %s. \\n`grep sample_status *.json`\\n" | (cat -; uuencode %s %s) | ssmtp -vvv %s >/dev/null 2>&1\n' % (self.sample_json['sample_name'], "pass", xlsx_file, xlsx_file.split('/')[-1], email)
					runCommandLine(email_command)
			# just send the email for now.
#			# I will copy the .xlsx file to every run of the sample
#			for run in self.sample_json['runs']:
#				run_json = json.load(open(run))
#				if 'server_ip' in run_json and 'orig_filepath_plugin_dir' in run_json:
#					copy_command = "scp %s ionadmin@%s:%s "%(xlsx_file, run_json['server_ip'], run_json['orig_filepath_plugin_dir'])
#					status = runCommandLine(copy_command)
#					if status == 0:
#						print "Copied the QC.xlsx file back to %s successfully!  %s"%(run_json['proton'], copy_command)
#					else:
#						print "Failed to copy the QC.xlsx file back to %s...  %s"%(run_json['proton'], copy_command)
#					# try to copy the log file back as well.
#					copy_command = "scp %s/sge.log ionadmin@%s:%s/QC.log "%(self.sample_json['sample_folder'], run_json['server_ip'], run_json['orig_filepath_plugin_dir'])
#					status = runCommandLine(copy_command)
#					# try to add the log file to the plugin's log file.
#					# this didn't work... 
#					#copy_command = "ssh ionadmin@%s:%s/QC.log 'cat %s/QC.log >> %s/drmaa_stdout.txt"%(run_json['server_ip'], run_json['orig_filepath_plugin_dir'], run_json['orig_filepath_plugin_dir'], run_json['orig_filepath_plugin_dir'])
#					#status = runCommandLine(copy_command)

	# if the update_json flag is specified, then update the cutoffs found in the normal json file.
	def update_cutoffs(self):
		# load the json file
		update_json = json.load(open(self.options.update_cutoffs))
		# set the cutoff settings to the example json's cutoff settings
		self.sample_json['analysis']['settings']['cutoffs'] = update_json['analysis']['settings']['cutoffs']
		# write the updated sample's json file.
		write_json(self.options.json, self.sample_json)

	# move the old 3x3 tables to the flag "old_GTs" 
	def recalc_3x3_tables(self):
		self.sample_json['sample_status'] = 'pushed'
		# load the output QC json. will be used to check if this combination has already been made.
		qc_json_data = {}
		if os.path.isfile(self.sample_json['results_qc_json']):
			qc_json_data = json.load(open(self.sample_json['results_qc_json']))
		# if the user specified to recalculate the 3x3 tables, do that here.
		if self.options.recalc_3x3_tables and 'QC_comparisons' in qc_json_data:
			# rearrange the old 3x3 tables to calculate the new 3x3 tables usingn the updated GT cutoffs
			qc_json_data['old_GTs'] = qc_json_data['QC_comparisons']
			del qc_json_data['QC_comparisons']
			write_json(self.sample_json['results_qc_json'], qc_json_data)


if __name__ == '__main__':
	
	# set up the option parser
	parser = OptionParser()
	
	# add the options to parse
	parser.add_option('-j', '--json', dest='json', help="A sample's json file which contains the necessary options and list of runs to QC with each other")
	parser.add_option('-e', '--email', dest='email', action='store_true', help="Option to send an email when job finishes")
	parser.add_option('-q', '--qc_all', dest='qc_all', action='store_true', help="Generate the 3x3 tables for all run comparisons, even if they fail.")
	parser.add_option('-p', '--pass_fail', dest='pass_fail', action='store_true', help="Overwrite the 'pass/fail' status of each run according to the cutoffs found in the json file. Normally this step is skipped if all runs have finished the QC process, but this option will overwrite the 'pass/fail' status found.")
	parser.add_option('-r', '--recalc_3x3_tables', dest='recalc_3x3_tables', action='store_true', help="recalculate the 3x3 tables (original use was for the new GT cutoffs")
	parser.add_option('-u', '--update_cutoffs', dest='update_cutoffs', help="Change the cutoffs found in the JSON file using an example json file with the corrected cutoffs. Will be done before anything else in the script.")
	(options, args) = parser.parse_args()

	# check to make sure the inputs are valid
	if not options.json:
		print "USAGE_ERROR: --json is required"
		parser.print_help()
		sys.exit(8)
	if not os.path.isfile(options.json):
		print "ERROR: the json file '%s' is not found"%options.json
		parser.print_help()
		sys.exit(4)
	if options.update_cutoffs and not os.path.isfile(options.update_cutoffs):
		print "ERROR: the update_cutoffs json file '%s' is not found"%options.update_cutoffs
		parser.print_help()
		sys.exit(4)

	try:
		qc_sample = QC_Sample(options)

		# if the update_json flag is specified, then update the cutoffs found in the normal json file.
		if options.update_cutoffs:
			qc_sample.update_cutoffs()

		# if the recalc_3x3_tables flag is specified, then rearrange the results_QC.json file so that the 3x3 tables will be recalculated.
		if options.recalc_3x3_tables:
			qc_sample.recalc_3x3_tables()

		# QC and merge all of the runs
		qc_sample.QC_merge_runs()

	except ValueError as e:
		#print sys.exc_traceback
		print traceback.format_exc().strip()
		print "Error: the JSON file is probably malformed"

