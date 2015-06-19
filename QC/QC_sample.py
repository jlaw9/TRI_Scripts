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
		self.no_errors = True
		self.cleanup_sample = Cleanup()
		self.options = options
		self.sample_json = json.load(open(options.json))
		# initialize the QC_Run object
		self.qc_run = QC_Run(self.sample_json, options.recalc_3x3_tables)
		self.__QCDirectory = "%s/QC"%self.sample_json['analysis']['software_directory']
	
	# will find all of the runs in a sample and QC them with each other
	def QC_merge_runs(self):
		# if this is a germline sample, QC all of the normal runs with each other.
		if self.sample_json['sample_type'] == 'germline':
			self.QC_germline()
	
		# if this is a tumor_normal sample, find the normal and tumor runs, and then QC them with each other.
		elif self.sample_json['sample_type'] == 'tumor_normal':
			self.QC_tumor_normal()
	
		# print the final status
		if self.no_errors == False or self.qc_run.no_errors == False:
			sys.stderr.write("%s finished with errors. See %s/sge.log for more details"%(self.sample_json['sample_name'], self.sample_json['output_folder']))
			self.sample_json['sample_status'] == 'failed'
			write_json(self.sample_json['json_file'], self.sample_json)
			sys.exit(1)
		else:
			print "%s finished with no errors"%(self.sample_json['sample_name'])
	
		# write the sample json file
		write_json(self.sample_json['json_file'], self.sample_json)
	
		# make the excel spreadsheet containing the data and copy it back to the proton
		self._make_xlsx()

	# if this is a germline sample, QC all of the normal runs with each other.
	def QC_germline(self):
		# Use the sample_status here to not re-run the QC and to not overwrite run status. The 'sample_status' should be reset to 'pushed' when new runs are pushed..
		#if self.sample_json['sample_status'] != 'pending_merge' and self.sample_json['sample_status'] != 'pending_3x3_review' and self.sample_json['sample_status'] != 'merged':
		# if the user specified the '--pass_fail' option, then run this part still
		if self.sample_json['sample_status'] == 'pushed' or self.options.pass_fail or self.options.qc_all:
			# QC the normal runs with each other
			self.QC_runs(self.sample_json['runs'])
	
		# what if there is only one run that passes all of the metrics? It should be marked as the 'final_json' and have the 'pass_fail_merged' flag marked as pass.
		# make the merger
		merger = Merger(self.sample_json, self.options.recalc_3x3_tables)
		# Check to see if the normal runs are ready to be merged.
		self.sample_json, merge = merger.check_merge(self.sample_json['runs'])
		if merge != True:
			if 'final_json' in self.sample_json:
				# update the final run status
				merger.update_merged_run_status(self.sample_json['final_json'])
		elif merge == True:
			# merge the normal and/or tumor runs. Will only merge the passing runs with each other.
			self.sample_json = merger.merge_runs('germline')
	
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
	
		# copy the final run's VCF file to the final_dir if it passes the "merged" coverage flag
		if 'final_json' in self.sample_json:
			final_json = json.load(open(self.sample_json['final_json']))
			if final_json['pass_fail_merged_status'] == 'pass':
				final_vcf = glob.glob("%s/*.vcf"%final_json['run_folder'])[0]
				final_project_dir = "/home/ionadmin/jeff/%s_Final_VCFs"%(self.sample_json['project'])
				print "copying %s to %s"%(final_vcf, final_project_dir)
				# check to make sure the final dir exists.
				if not os.path.isdir(final_project_dir):
					os.mkdir(final_project_dir)
				shutil.copy(final_vcf, "%s/%s.vcf"%(final_project_dir, self.sample_json['sample_name']))
				# now push the sample to s3 storage
				if self.sample_json['project'] == 'Einstein':
					print "pushing %s to amazon s3 storage"%self.sample_json['sample_name']
					self.push_sample_to_s3(final_json)

	# if this is a tumor_normal sample, find the normal and tumor runs, and then QC them with each other.
	def QC_tumor_normal(self):
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
	
			# make the merger
			merger = Merger(self.sample_json, self.options.recalc_3x3_tables)
			# Check to see if the normal runs are ready to be merged.
			self.sample_json, merge_normal = merger.check_merge(normal_runs, 'Normal/', 'normal_')
			if merge_normal == True:
				# merge the normal and/or tumor runs. Will only merge the passing runs with each other.
				self.sample_json = merger.merge_runs('normal', 'Normal_', 'normal_')
	
			# Check to see if the tumor runs are ready to be merged.
			self.sample_json, merge_tumor = merger.check_merge(tumor_runs, 'Tumor/', 'tumor_')
			if merge_tumor == True:
				self.sample_json = merger.merge_runs('tumor', 'Tumor_', 'tumor_')
	
			# If any runs were merged, QC them. If there are only 1 normal and tumor run, they won't be QCd again. 
			#if normal_merge_dir != '' or tumor_merge_dir != '' or (len(normal_passing_bams) == 1 and len(tumor_passing_bams) == 1):	
			# now QC the tumor and normal merged bams together if both normal and tumor runs are ready.
				# To only QC all for the actual merged runs (PNET), change the 'final' part to 'merged'.
			# The 'final_normal_json' and 'final_tumor_json' flags are set by merger.py in the function check_merge, line 157
			#if (merge_normal or merge_tumor) and ('merged_normal_json' in self.sample_json and 'merged_tumor_json' in self.sample_json):
			if 'final_normal_json' in self.sample_json and 'final_tumor_json' in self.sample_json:
				self.sample_json, qc_json = self.qc_run.QC_2Runs(self.sample_json, self.sample_json['final_normal_json'], self.sample_json['final_tumor_json'], 'normal_', 'tumor_', '_merged')
				self.sample_json, merged_perc_avail_bases = self.qc_run.update_3x3_runs_status(self.sample_json, self.sample_json['final_normal_json'], self.sample_json['final_tumor_json'], qc_json)
				# update the merged run status 
				merger.update_merged_run_status(self.sample_json['final_normal_json'], merged_perc_avail_bases)
				merger.update_merged_run_status(self.sample_json['final_tumor_json'], merged_perc_avail_bases)
	
				# cleanup the individual run bam files
				if merged_perc_avail_bases > .9:
					final_qc_dir = "%s/all%svs%s"%(self.sample_json['qc_folder'], json.load(open(self.sample_json['final_normal_json']))['run_name'], json.load(open(self.sample_json['final_tumor_json']))['run_name'])
					# annotate the final somatic variants
					command = "bash %s/Somatic_Variants/somatic_variants.sh %s %s"%(self.sample_json['analysis']['software_directory'], final_qc_dir, self.sample_json['sample_name'])
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
		pending_runs, passing_runs = self.qc_run.get_runs_status(runs)
		if len(passing_runs) == 1:
			self.sample_json['sample_status'] = 'pending_merge'
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
					somatic_variants = "%s/%s_somatic.xlsx"%(self.sample_json['qc_folder'], self.sample_json['sample_name'])
					if os.path.isfile(somatic_variants):
						email_command = '\tprintf "%s finished with a status of %s. \\n`grep sample_status *.json`\\n" | (cat -; uuencode %s %s; uuencode %s %s) | ssmtp -vvv %s >/dev/null 2>&1\n' % (self.sample_json['sample_name'], "pass", xlsx_file, xlsx_file.split('/')[-1], somatic_variants, somatic_variants, email)
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

	# send an email with the specified attachments
	# TODO finish this function
	def _send_email(self, emails, attachments, status):
		for email in emails:
			# this command will email the status of the sample, and attach the excel spreadsheet and somatic variants if it is found.
			email_command = '\tprintf "%s finished with a status of %s. \\n`grep sample_status *.json`\\n" | ssmtp -vvv %s >/dev/null 2>&1\n' % (self.sample_json['sample_name'], status, email)
			runCommandLine(email_command)

	# pushes the final run or merged files to amazon s3 storage.
	def push_sample_to_s3(self, final_json):
		# first get all of the files to push
		final_vcf = glob.glob("%s/*.vcf"%final_json['run_folder'])[0]
		target_vcf = "Einstein/%s/%s.vcf"%(self.sample_json['sample_name'], self.sample_json['sample_name'])
		final_cov = glob.glob("%s/*.amplicon.cov.xls"%(final_json['run_folder']))[0]
		target_cov = "Einstein/%s/%s"%(self.sample_json['sample_name'], final_cov.split('/')[-1])
		final_bam = "%s/%s"%(final_json['run_folder'], final_json['analysis']['files'][0])
		target_bam = "Einstein/%s/%s"%(self.sample_json['sample_name'], final_bam.split('/')[-1])
		final_bai = final_bam + ".bai"
		target_bai = target_bam + ".bai"
		final_json_file = final_json['json_file']
		target_json_file = "Einstein/%s/%s"%(self.sample_json['sample_name'], final_json_file.split('/')[-1])
	
		# call the push_files script to push each file to s3 storage
		status = os.system("bash /rawdata/scripts/TRI_Dev/push_files_s3.sh " + \
				"%s %s "%(final_vcf, target_vcf) + "%s %s "%(final_cov, target_cov) + \
				"%s %s "%(final_bam, target_bam) + "%s %s "%(final_bai, target_bai) + \
				"%s %s "%(final_json_file, target_json_file))
		if status != 0:
			self.no_errors = False
			print "ERROR: unable to push the sample to s3 storage"

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

	# get the alignment statistics for each run or merged bam file.
	def get_alignment_stats(self):
		# TEMP fix the runs.
		runs = []
		for run in glob.glob("%s/Normal/N-[0-9]/*.json"%self.sample_json['sample_folder']):
			runs.append(run)
		for run in glob.glob("%s/Tumor/T-[0-9]/*.json"%self.sample_json['sample_folder']):
			runs.append(run)
		self.sample_json['runs'] = runs
		write_json(self.sample_json['json_file'], self.sample_json)
		#
		## now get the alignment statistics
		for run in self.sample_json['runs']:
			print "getting alignment_status for: %s"%run
			Align_Stats(run)
#		if 'merged_normal_json' in self.sample_json:
#			normal_merged_path = json.load(open(self.sample_json['merged_normal_json']))['run_folder']
#			if not os.path.isfile('%s/Analysis_Files/ionstats_alignment.json'%json.load(open(self.sample_json['merged_normal_json']))['run_folder']):
#				# first fix the header of the merged.bam file
#				normal_merged_bam = "%s/%s"%(normal_merged_path, json.load(open(self.sample_json['merged_normal_json']))['analysis']['files'][0])
#				print "fixing header for %s"%normal_merged_bam
#				correct_header_command = "samtools view -H %s > %s/merged.header.sam "%(normal_merged_bam, normal_merged_path)
#				if runCommandLine(correct_header_command) != 0:
#					print "ERROR: samtools view -H failed!"
#					sys.exit(1)
#
#				# move the old bam
#				old_normal_merged_bam = "%s/bad_header.bam"%normal_merged_path
#				shutil.move(normal_merged_bam, old_normal_merged_bam)
#				# this command deletes the KS: tag!! not good! I don't know why but some headers are tab delimited, and some are not it seems.
#				sed_command = 'sed -E "s/SM:[^:]*:/SM:%s\tKS:/" %s/merged.header.sam > %s/merged.headerCorrected.sam'%(self.sample_json['sample_name'], normal_merged_path, normal_merged_path)
#				if runCommandLine(sed_command) != 0:
#					print "ERROR: sed command failed!"
#					sys.exit(1)
#				# write the new header to merged.bam
#				reheader_command = "samtools reheader %s/merged.headerCorrected.sam %s > %s "%(normal_merged_path, old_normal_merged_bam, normal_merged_bam)
#				if runCommandLine(reheader_command) != 0:
#					print "ERROR: sed command failed!"
#					sys.exit(1)
#				# make a new index file
#				runCommandLine("samtools index %s"%normal_merged_bam)
#				#remove the old bam
#				os.remove(old_normal_merged_bam)
#				os.remove("%s/merged.headerCorrected.sam"%normal_merged_path)
#				os.remove("%s/merged.header.sam"%normal_merged_path)
#				# then get the ionstats
#				Align_Stats(self.sample_json['merged_normal_json'])
#		if 'merged_tumor_json' in self.sample_json:
#			tumor_merged_path = json.load(open(self.sample_json['merged_tumor_json']))['run_folder']
#			if not os.path.isfile('%s/Analysis_Files/ionstats_alignment.json'%json.load(open(self.sample_json['merged_tumor_json']))['run_folder']):
#				# first fix the header of the merged.bam file
#				tumor_merged_bam = "%s/%s"%(tumor_merged_path, json.load(open(self.sample_json['merged_tumor_json']))['analysis']['files'][0])
#				print "fixing header for %s"%tumor_merged_bam
#				correct_header_command = "samtools view -H %s > %s/merged.header.sam "%(tumor_merged_bam, tumor_merged_path)
#				if runCommandLine(correct_header_command) != 0:
#					print "ERROR: samtools view -H failed!"
#					sys.exit(1)
#
#				# move the old bam
#				old_tumor_merged_bam = "%s/bad_header.bam"%tumor_merged_path
#				shutil.move(tumor_merged_bam, old_tumor_merged_bam)
#				# this command deletes the KS: tag!! not good! I don't know why but some headers are tab delimited, and some are not it seems.
#				sed_command = 'sed -E "s/SM:[^:]*:/SM:%s\tKS:/" %s/merged.header.sam > %s/merged.headerCorrected.sam'%(self.sample_json['sample_name'], tumor_merged_path, tumor_merged_path)
#				if runCommandLine(sed_command) != 0:
#					print "ERROR: sed command failed!"
#					sys.exit(1)
#				# write the new header to merged.bam
#				reheader_command = "samtools reheader %s/merged.headerCorrected.sam %s > %s "%(tumor_merged_path, old_tumor_merged_bam, tumor_merged_bam)
#				if runCommandLine(reheader_command) != 0:
#					print "ERROR: sed command failed!"
#					sys.exit(1)
#				# make a new index file
#				runCommandLine("samtools index %s"%tumor_merged_bam)
#				#remove the old bam
#				os.remove(old_tumor_merged_bam)
#				os.remove("%s/merged.headerCorrected.sam"%tumor_merged_path)
#				os.remove("%s/merged.header.sam"%tumor_merged_path)
#				# then get the ionstats
#				Align_Stats(self.sample_json['merged_tumor_json'])

		# copy the xlsx file here because it didn't get copied for a lot of samples
		#self._make_xlsx()

	# subset out the 718 gene set from the final merged PNET 3x3 tables
	def get_718_subset(self):
		# add the path to the 718 subset:
		self.sample_json['analysis']['settings']['subset_bed'] = '/rawdata/support_files/BED/PNET/AmpliSeqExome_PNET_subset.bed'
		self.sample_json['analysis']['settings']['chromosomes_to_analyze_merged'] = ['all', '718']
		if 'results_qc_json' not in self.sample_json and 'results_QC_json' in self.sample_json:
			self.sample_json['results_qc_json']	= self.sample_json['results_QC_json']
		self.sample_json['emails'] = ['jlaw@childhooddiseases.org']

		#self.sample_json['sample_status'] = 'pending_merge'
		write_json(self.sample_json['json_file'], self.sample_json)
		#Normal_Merged1vsTumor_Merged1
		#NMerge1vsTMerged1
		qc_comp_dir = ''
		if os.path.isdir("%s/allNMerged1vsTMerged1"%self.sample_json['qc_folder']):
			qc_comp_dir = "%s/allNMerged1vsTMerged1"%self.sample_json['qc_folder']
			qc_comp = "NMerged1vsTMerged1"
		elif os.path.isdir("%s/allNormal_Merged1vsTumor_Merged1"%self.sample_json['qc_folder']):
			qc_comp_dir = "%s/allNormal_Merged1vsTumor_Merged1"%self.sample_json['qc_folder']
			qc_comp = "Normal_Merged1vsTumor_Merged1"
		if qc_comp_dir != '':
			results_qc_json = json.load(open(self.sample_json['results_qc_json']))
			# fix the name of the folder and the name in the results_qc_json
			normal_merged_name = json.load(open(self.sample_json['merged_normal_json']))['run_name']
			tumor_merged_name = json.load(open(self.sample_json['merged_tumor_json']))['run_name']
			new_qc_comp = "%svs%s"%(normal_merged_name, tumor_merged_name)
			print "moving %s to %s"%(qc_comp_dir, "%s/all%s"%(self.sample_json['qc_folder'], new_qc_comp))
			shutil.move(qc_comp_dir, "%s/all%s"%(self.sample_json['qc_folder'], new_qc_comp))
	
			results_qc_json = json.load(open(self.sample_json['results_qc_json']))
			new_qc_comp_dict = results_qc_json['QC_comparisons']['all']['normal_tumor'][qc_comp]
			del results_qc_json['QC_comparisons']['all']['normal_tumor'][qc_comp]
	
			results_qc_json['QC_comparisons']['all']['normal_tumor'][new_qc_comp] = new_qc_comp_dict
			results_qc_json['sample_name'] = self.sample_json['sample_name']
			results_qc_json['sample'] = self.sample_json['sample_name']
			write_json(self.sample_json['results_qc_json'], results_qc_json)
			
		if 'merged_normal_json' in self.sample_json and 'merged_tumor_json' in self.sample_json:
			print "Running QC_2Runs"
			self.sample_json, qc_json = self.qc_run.QC_2Runs(self.sample_json, self.sample_json['merged_normal_json'], self.sample_json['merged_tumor_json'], 'normal_', 'tumor_', '_merged')
			print 'done'
			# done for now
			self._make_xlsx()

	def overlap(self):
		# add the merged_perc_avail_bases
			# fix the name of the folder and the name in the results_qc_json
			if 'merged_normal_json' in self.sample_json and 'merged_tumor_json' in self.sample_json:
				results_qc_json = json.load(open(self.sample_json['results_qc_json']))
				normal_merged_json = json.load(open(self.sample_json['merged_normal_json']))
				tumor_merged_json = json.load(open(self.sample_json['merged_tumor_json']))
				qc_comp = "%svs%s"%(normal_merged_json['run_name'], tumor_merged_json['run_name'])
				perc_avail_bases = results_qc_json['QC_comparisons']['all']['normal_tumor'][qc_comp]['perc_avail_bases']
				normal_merged_json['run_data']['merged_perc_avail_bases'] = perc_avail_bases
				tumor_merged_json['run_data']['merged_perc_avail_bases'] = perc_avail_bases
				write_json(self.sample_json['merged_normal_json'], normal_merged_json)
				write_json(self.sample_json['merged_tumor_json'], tumor_merged_json)

	def A_227(self):
		# compare the Normal merged file to all of the other tumor combinations
		# Separate the runs into tumor and normal lists
		#normal_runs, tumor_runs = self.getTumor_Normal()
		#for tumor_run in tumor_runs:
		#	# generate the 3x3 tables for only chr1. 
		#	self.sample_json, qc_json = self.qc_run.QC_2Runs(self.sample_json, self.sample_json['merged_normal_json'], tumor_run, 'normal_', 'tumor_', '_merged')
		self.sample_json, qc_json = self.qc_run.QC_2Runs(self.sample_json, self.sample_json['merged_normal_json'], self.sample_json['merged_tumor_json'], 'normal_', 'tumor_', '_merged')
		# merge tumor runs 1-3, 4-5, and 7-8.
        #tumor_1_2_3 = ["/mnt/Despina/projects/PNET/A_227/Tumor/T-1/A_227_T-1.json", "/mnt/Despina/projects/PNET/A_227/Tumor/T-2/A_227_T-2.json", "/mnt/Despina/projects/PNET/A_227/Tumor/T-3/A_227_T-3.json"]
        #tumor_7_8 = ["/mnt/Despina/projects/PNET/A_227/Tumor/T-7/A_227_T-7.json", "/mnt/Despina/projects/PNET/A_227/Tumor/T-8/A_227_T-8.json"]

		# done for now
		self._make_xlsx()
			
	def somatic_variants(self):
		# get the somatic variants from samples that pass the final overlapping coverage cutoff
		if 'final_normal_json' in self.sample_json and 'final_tumor_json' in self.sample_json:
			final_normal_json = json.load(open(self.sample_json['final_normal_json']))
			final_tumor_json = json.load(open(self.sample_json['final_tumor_json']))
			if "pass_fail_merged_status" in final_normal_json and final_normal_json["pass_fail_merged_status"] == 'pass':
				# get the path to the final QC comparison dir
				qc_comp_dir = "%s/QC/all%svs%s"%(self.sample_json['sample_folder'], final_normal_json['run_name'], final_tumor_json['run_name'])
				# get the somatic variants using the somatic_variants.sh script which utilizes Ozlem's scripts.
				#TODO "import" Ozlem's scripts into this pipeline
				command = "bash %s/Somatic_Variants/somatic_variants.sh %s %s"%(self.sample_json['analysis']['software_directory'], qc_comp_dir, self.sample_json['sample_name'])
				result = os.system(command)
				if result != 0:
					self.no_errors = False
				self._make_xlsx()

	def change_stringency(self):
		self.sample_json['analysis']['settings']['normal_tvc_json'] = "/rawdata/support_files/parameter_sets/Parameter_Tests/ampliseq_germline_lowstringency_pgm_parameters_jingwei_edits.json"
		self.sample_json['analysis']['settings']['tumor_tvc_json'] = "/rawdata/support_files/parameter_sets/Parameter_Tests/ch1_somatic_lowstringency_pgm_parameters_jingwei_edits.json"
		# initialize the QC_Run object
		qc_run_diff_settings = QC_Run(self.sample_json, recalc_3x3_tables=False)
		# change the stringency of the somatic analysis to see if we can find more somatic variants
		if 'final_normal_json' in self.sample_json and 'final_tumor_json' in self.sample_json:
			final_normal_json = json.load(open(self.sample_json['final_normal_json']))
			final_tumor_json = json.load(open(self.sample_json['final_tumor_json']))
			# check to see if the final runs passed the merged cutoff
			if "pass_fail_merged_status" in final_normal_json and final_normal_json["pass_fail_merged_status"] == 'pass':
				for final_json in [final_normal_json, final_tumor_json]:
					# re-run TVC on the individual bam files
					vcfs = glob.glob("%s/*.vcf"%final_json['run_folder'])
#					if len(vcfs) > 0: 
#						shutil.move(vcfs[0], "%s/Analysis_Files/4.2_TSVC_variants_High_String_Jingwei_edits.vcf"%final_json['run_folder']) 
					qc_run_diff_settings.runTVC_COV(final_json['json_file'], "%s_"%final_json['run_type'])
	
				# now fix the results_qc_json
				# move the results_qc_json analysis to a different "chromosome"
#				results_qc_json = json.load(open(self.sample_json['results_qc_json']))
#				qc_comp = "%svs%s"%(final_normal_json['run_name'], final_tumor_json['run_name'])
#				old_qc_comp_dict = results_qc_json['QC_comparisons']['all']['normal_tumor'][qc_comp]
#				del results_qc_json['QC_comparisons']['all']['normal_tumor'][qc_comp]
#		
#				results_qc_json['QC_comparisons']['High_String_Jingwei_Edits'] = {'normal_tumor': {qc_comp: old_qc_comp_dict}}
#				write_json(self.sample_json['results_qc_json'], results_qc_json)

				# get the path to the final QC comparison dir
				qc_comp_dir = "%s/QC/all%svs%s"%(self.sample_json['sample_folder'], final_normal_json['run_name'], final_tumor_json['run_name'])
#				high_string_qc_comp_dir = "%s/QC/High_String_Jingwei_Edits_all%svs%s"%(self.sample_json['sample_folder'], final_normal_json['run_name'], final_tumor_json['run_name'])
#				shutil.move(qc_comp_dir, high_string_qc_comp_dir)

				# QC the 2 final runs
				self.sample_json, qc_json = qc_run_diff_settings.QC_2Runs(self.sample_json, self.sample_json['final_normal_json'], self.sample_json['final_tumor_json'], 'normal_', 'tumor_', '_merged')

				# get the somatic variants using the somatic_variants.sh script which utilizes Ozlem's scripts.
				#TODO "import" Ozlem's scripts into this pipeline
				command = "bash %s/Somatic_Variants/somatic_variants.sh %s %s"%(self.sample_json['analysis']['software_directory'], qc_comp_dir, self.sample_json['sample_name'])
				result = os.system(command)
				if result != 0:
					self.no_errors = False
				#self._make_xlsx()

	def samtools_gatk(self):
		# 7 gene panel for PNET analysis
		Seven_Gene_Bed = "/rawdata/support_files/BED/PNET/7PNET_Genes_amplicons.bed"
		# run samtools and GATK on the SEGA dataset
		if 'final_normal_json' in self.sample_json and 'final_tumor_json' in self.sample_json:
			final_normal_json = json.load(open(self.sample_json['final_normal_json']))
			final_tumor_json = json.load(open(self.sample_json['final_tumor_json']))
			for final_json in [final_normal_json, final_tumor_json]:
				# if the analysis has not already been done, do it
				#if not os.path.isfile("%s/Analysis_Files/gatk_filtered_snps.vcf"%(final_json['run_folder'])) or not os.path.isfile("%s/Analysis_Files/7Genes_TSVC_variants.vcf"%(final_json['run_folder'])):
				bam = "%s/%s"%(final_json['run_folder'], final_json['analysis']['files'][0])

				# subset the 7 genes from the VCF file
				#vcf = glob.glob("%s/*.vcf"%(final_json['run_folder']))[0]
				#command = "bedtools intersect -header -a %s -b %s > %s/Analysis_Files/7Genes_TSVC_variants.vcf"%(vcf, Seven_Gene_Bed, final_json['run_folder'])
				#tvc_result = runCommandLine(command)
					
				# subset the 7 genes from the bam file to be uploaded to IR.
				#command = "samtools view -L /rawdata/support_files/BED/PNET/7PNET_Genes_amplicons.bed %s -b > %s/Analysis_Files/%s_%s_7genes_06032015.bam"%(bam, final_json['run_folder'], self.sample_json['sample_name'], final_json['run_type'])
				#subset_result = runCommandLine(command)
		
				# run samtools
				#command = "bash /home/ionadmin/TRI_Scripts/Variants/Samtools/run_samtools.sh %s %s %s/Analysis_Files"%(bam, Seven_Gene_Bed, final_json['run_folder'])
				#samtools_result = runCommandLine(command)
		
				# run GATK
				command = "bash /home/ionadmin/TRI_Scripts/Variants/GATK/run_gatk.sh %s %s %s/Analysis_Files"%(bam, Seven_Gene_Bed, final_json['run_folder'])
				gatk_result = runCommandLine(command)
		
				# I would push the bam file to IR now, but that would require a password...
				# Lastly, generate the Venn Diagrams. I need to copy the R code here and figure out how to copy the IR TVC results back to here as well.
			#	if samtools_result != 0 or gatk_result != 0 or tvc_result != 0 or subset_result != 0:
			#		sys.stderr.write("Error: samtools or gatk failed!\n")
			#		self.no_errors = False
			# done for now
			#self._send_email("jlaw@childhooddiseases.org")
			#self._make_xlsx()

	def variant_call(self):
		# re-run TVC to see which variants were filtered out
		#final_normal_dir = json.load(open(self.sample_json['final_normal_json']))['run_folder']
		final_tumor_dir = json.load(open(self.sample_json['final_tumor_json']))['run_folder']
		#runCommandLine("mv %s/4.2_TSVC_variants.vcf %s/Analysis_Folder"%(final_normal_dir, final_normal_dir))
		#runCommandLine("mv %s/4.2_TSVC_variants.vcf %s/Analysis_Folder"%(final_tumor_dir, final_tumor_dir))
		#self.qc_run.runTVC_COV(self.sample_json['final_normal_json'], "normal_")
		self.qc_run.runTVC_COV(self.sample_json['final_tumor_json'], "tumor_")

	def einstein_merged(self):
		# first find the merged dir
		merged_dir = "%s/Merged"%self.sample_json['sample_folder'] 
		if os.path.isdir(merged_dir):
			# if the merged file has >= 30x coverage, copy the VCF file
			merged_amp = glob.glob("%s/*.amplicon.cov.xls"%merged_dir)[0]
			command = "tail -n +2 %s | awk -v cutoff=30 '{ if ($10 >= cutoff) printf \".\"}' | wc -c"%merged_amp
			amp_cov = runCommandLine(command, get_output=True)
			command = "tail -n +2 %s | wc -l"%merged_amp
			num_amps = runCommandLine(command, get_output=True)
			if (float(amp_cov) / float(num_amps)) >= 0.9:
				# copy the VCF file
				var_file = glob.glob("%s/*.vcf"%merged_dir)[0]
				command = "cp %s /home/ionadmin/jeff/Einstein_passVCF/%s.vcf"%(var_file, self.sample_json['sample_name'])
				runCommandLine(command)
				# this sample should now be ready to push to amazons3 storage.


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
	parser.add_option('-s', '--get_alignment_stats', dest='get_alignment_stats', action='store_true', help="Gather the alignment stats for the runs and merged files")
	parser.add_option('-b', '--get_718_subset', dest='get_718_subset', action='store_true', help="Subset the 718 gene set for the 3x3 tables")
	parser.add_option('-o', '--overlap', dest='overlap', action='store_true', help="add teh % overlapping bases metric to the individual runs.")
	parser.add_option('-A', '--A_227', dest='A_227', action='store_true', help="compare Normal A_227 merged with Tumor runs")
	parser.add_option('-S', '--somatic_variants', dest='somatic_variants', action='store_true', help="annotate the somatic variants and store them in ~/jeff/Lung_Somatic/[sample]_somatic.vcf. Need to make this option more applicable to other projects")
	parser.add_option('-c', '--change_stringency', dest='change_stringency', action='store_true', help="Re-run the final somatic variant calling using a different stringency.")
	parser.add_option('-g', '--samtools_gatk', dest='samtools_gatk', action='store_true', help="Run Samtools and GATK")
	parser.add_option('-V', '--variant_call', dest='variant_call', action='store_true', help="Re-run the variant caller to see which variants were filtered out.")
	parser.add_option('-E', '--einstein_merged', dest='einstein_merged', action='store_true', help="get run metrics of the einstein merged files")
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

		# if the user specified, then only get_alignment_stats
		if options.get_alignment_stats:
			qc_sample.get_alignment_stats()
		elif options.get_718_subset:
			qc_sample.get_718_subset()
		elif options.overlap:
			qc_sample.overlap()
		elif options.A_227:
			qc_sample.A_227()
		elif options.change_stringency:
			qc_sample.change_stringency()
		elif options.somatic_variants:
			qc_sample.somatic_variants()
		elif options.samtools_gatk:
			qc_sample.samtools_gatk()
		elif options.variant_call:
			qc_sample.variant_call()
		elif options.einstein_merged:
			qc_sample.einstein_merged()
		else:
			# QC and merge all of the runs
			qc_sample.QC_merge_runs()

	except ValueError as e:
		#print sys.exc_traceback
		print traceback.format_exc().strip()
		print "Error: the JSON file is probably malformed"

