#! /usr/bin/env python2.7

# Goal: Create pretty 3x3 tables comparing mulitple runs, or Tumor/Normal pairs 

import sys
import os
import re
import subprocess
import json
import QC_stats
from optparse import OptionParser
try:
	import xlsxwriter
	from xlsxwriter.utility import xl_rowcol_to_cell
except ImportError:
	print "xlsxwriter is installed on python2.7 on triton. Use python2.7", sys.argv[0]
	sys.exit(1)


class XLSX_Writer():
	def __init__(self, workbook):
		#Dictionary containing all of the different self.formats
		# initialize self.formats with a blank string in case there are times where there is no format
		# this dictionary will hold all of the self.formats to be used by QC_genSheets
		self.formats = {}
		self.workbook = workbook
		self._setupWorkbook()


	# @param out_file_name the name of the outuput xlsx
	##### @returns a format of dictionary
	def _setupWorkbook(self):
		# The color for each alternating sample
		alt_sample_color = "#d5e8f8" # This color is considered azure by the color wheel. Go to this website to choose a new color

		# Bonnie requested that I change the default font size to 9.
		text_font_size = 9
		
		# set the formats for the run metrics sheet
		self.formats['header_format'] = self.workbook.add_format({'bold': True, 'font_size': 12, 'align': 'center', 'text_wrap': True})
		self.formats[''] = self.workbook.add_format({'font_size': text_font_size, 'align': 'center'})
		self.formats['center'] = self.workbook.add_format({'font_size': text_font_size, 'align': 'center'})
		self.formats['_azure'] = self.workbook.add_format({'font_size': text_font_size, 'align': 'center', 'bg_color': alt_sample_color})
		
		# set the number and perc formats
		self.formats['num_format'] = self.workbook.add_format({'num_format': '#,##0', 'font_size': text_font_size, 'align': 'center'})
		self.formats['perc_format'] = self.workbook.add_format({'num_format': '0.00%', 'font_size': text_font_size, 'align': 'center'})
		self.formats['dec3_format'] = self.workbook.add_format({'num_format': '0.000', 'font_size': text_font_size, 'align': 'center'})
		self.formats['num_format_azure'] = self.workbook.add_format({'num_format': '#,##0', 'font_size': text_font_size, 'align': 'center', 'bg_color': alt_sample_color})
		self.formats['perc_format_azure'] = self.workbook.add_format({'num_format': '0.00%', 'font_size': text_font_size, 'align': 'center', 'bg_color': alt_sample_color})
		self.formats['dec3_format_azure'] = self.workbook.add_format({'num_format': '0.000', 'font_size': text_font_size, 'align': 'center', 'bg_color': alt_sample_color})
		
		# formats for the 3x3 tables
		self.formats['header_format2'] = self.workbook.add_format({'bold': True, 'font_size': 14})
		self.formats['red'] = self.workbook.add_format({'font_size': text_font_size, 'align': 'center', 'bg_color': '#FF0000'})
		self.formats['top'] = self.workbook.add_format({'top': 2})
		self.formats['right'] = self.workbook.add_format({'right': 2})
		self.formats['bottom'] = self.workbook.add_format({'bottom': 2})
		self.formats['left'] = self.workbook.add_format({'left': 2})
		

	# Function to check if a value is greater than the maximum value (i.e. WT to HET is greater than 10 or something), then write it in red
	def _check_max_and_write(self, row, col, value, Max):
		if int(value) > int(Max):
			# write this cell in red
			self.MRsheet.write(row, col, value, self.formats['red'])
		else:
			self.MRsheet.write(row, col, value)
			# write this cell in normal white


	# Complex function to write cells. 
	# @param write_format is the format to write the cell in. 
	# @param Max will be the maximum threshold for writing in red, unless it's 0. If max if negative, it will be treated as a minimum threshold
	# MAx is not yet implemented.
	# @return returns 1 so there will be one less line of code (to incrament col). Maybe it could just incrament col, I would just rather not have global variables
	def _check_to_write(self, row, col, key, write_format, metrics, Sheet):
		if key in metrics:
			try:
				if re.search("=", write_format):
					if metrics[key] != "":
						cell1 = xl_rowcol_to_cell(row, col-2)
						cell2 = xl_rowcol_to_cell(row, col-1)
						Sheet.write_formula(row, col, "=%s-%s"%(cell1, cell2), self.formats[write_format[1:]])
				elif re.search("num_format", write_format):
					if not isinstance(metrics[key], int) and not isinstance(metrics[key], float):
						Sheet.write_number(row, col, int(metrics[key].replace(',','')), self.formats[write_format])
					else:
						Sheet.write_number(row, col, metrics[key], self.formats[write_format])
				elif re.search("perc_format", write_format):
					# begin and end amp_cov temporary fix
					if key == 'begin_amp_cov' or key == 'end_amp_cov':
						if float(metrics[key]) > 1:
							metrics[key] = metrics[key] / 2
					elif float(metrics[key]) > 2:
						metrics[key] = metrics[key] / 100
					Sheet.write_number(row, col, float(metrics[key]), self.formats[write_format])
				elif re.search("dec3_format", write_format):
					Sheet.write_number(row, col, float(metrics[key]), self.formats[write_format])
				# special case to write the formula for the +-10 bp col
				else:
					# if write_format is blank, then self.formats will also be blank
					Sheet.write(row, col, metrics[key], self.formats[write_format])
			except ValueError:
				Sheet.write(row, col, metrics[key], self.formats[write_format])
		return 1


	# little function to write a header cell.
	def _writeHeaderCell(self, col, text, width, Sheet):
		Sheet.write(0, col, text, self.formats['header_format'])
		Sheet.set_column(col,col,width, self.formats['center'])
		return 1


	# @param ex_json_data if it is None, then it will return 'false' in an if statement
	def _writeMergedMetricHeaders(self, ex_json_data=None):
		# RunSheet is where all of the metrics about each run will be written
		self.Merged = self.workbook.add_worksheet("Merged QC Metrics")
		self.Merged.freeze_panes(1,3)
		
		# First write the QC metrics for each run of each sample.
		# Write the header line. there could definitely be a better way of doing this, but this is what I figured out for now. Just comment and uncomment as needed.
		col = 0
		col += self._writeHeaderCell(col, "Sample #", None, self.Merged)
		# check to see if the N or T should be written
		if ex_json_data and 'sample_type' in ex_json_data and ex_json_data['sample_type'] == 'tumor_normal':
			col += self._writeHeaderCell(col, "Normal (N) or Tumor (T)", 7, self.Merged)
		col += self._writeHeaderCell(col, "Merged File Name", 30, self.Merged)
		col += self._writeHeaderCell(col, "Included Runs", 15, self.Merged)
		col += self._writeHeaderCell(col, "Total Aligned Basepairs", 12, self.Merged)
		col += self._writeHeaderCell(col, "Total AQ20 Basepairs", 12, self.Merged)
		col += self._writeHeaderCell(col, "Mean Read Length", None, self.Merged)
		col += self._writeHeaderCell(col, "Median Read Length", None, self.Merged)
		# see if we have the exp_read_length
		if ex_json_data and 'cutoffs' in ex_json_data['analysis']['settings'] and 'exp_median_read_length' in ex_json_data['analysis']['settings']['cutoffs']:
			col += self._writeHeaderCell(col, "%% expected read length (out of %d bp)"%ex_json_data['analysis']['settings']['cutoffs']['exp_median_read_length'], 12, self.Merged)
		else:
			col += self._writeHeaderCell(col, "% expected read length (out of XXX bp)", 12, self.Merged)
		col += self._writeHeaderCell(col, "Median Read Coverage Overall", None, self.Merged)
		col += self._writeHeaderCell(col, "% targeted bases >= 30x coverage", 12, self.Merged)
		col += self._writeHeaderCell(col, "% amplicons >= 30x coverage", 12, self.Merged)
		col += self._writeHeaderCell(col, "% amplicons >= 30x covered at bp +10 (considering fwd/rev read split)", 13, self.Merged)
		col += self._writeHeaderCell(col, "% amplicons >= 30x covered at bp n-10 (considering fwd/rev read split)", 13, self.Merged)
		col += self._writeHeaderCell(col, "% amplicons >= 30x covered at bp +10 (>=30x) that are also covered at n-10 bp (>=30x)", 13, self.Merged)
		col += self._writeHeaderCell(col, "Ts/Tv", None, self.Merged)
		col += self._writeHeaderCell(col, "# Total variants (single allele)", 13, self.Merged)
		col += self._writeHeaderCell(col, "# HET variants (single allele rates)", 13, self.Merged)
		col += self._writeHeaderCell(col, "# HOM variants (single allele rates)", 13, self.Merged)
		col += self._writeHeaderCell(col, "HET/HOM ratio (single allele rates)", 13, self.Merged)
		# check to see if the N or T should be written
		if ex_json_data and 'sample_type' in ex_json_data and ex_json_data['sample_type'] == 'tumor_normal':
			col += self._writeHeaderCell(col, "% Overlapping Bases", 8, self.Merged)
			col += self._writeHeaderCell(col, "Tumor/Normal Comparisons: No evidence of contamination", 10, self.Merged)
		try:
			# write this final header for tumor_normal comparisons
			if ex_json_data and 'merged_amp_cov' in ex_json_data['analysis']['settings']['cutoffs'] and ex_json_data['sample_type'] == "tumor_normal":
				col += self._writeHeaderCell(col, "Final Merged QC Status (PASS if the %% available bases is > %s in the final tumor/normal comparison"%ex_json_data['analysis']['settings']['cutoffs']['merged_amp_cov'], 15, self.Merged)
			elif ex_json_data and 'merged_amp_cov' in ex_json_data['analysis']['settings']['cutoffs'] and ex_json_data['sample_type'] == "germline":
				col += self._writeHeaderCell(col, "Final Merged QC Status (PASS if the %% amplicons covered is > %s)"%ex_json_data['analysis']['settings']['cutoffs']['merged_amp_cov'], 15, self.Merged)
		except (KeyError, TypeError):
			col += self._writeHeaderCell(col, "Final Merged QC Status", 10, self.Merged)
		self.Merged.set_column(col,col+20,12, self.formats['center'])
		
		self.Merged.set_row(0,100, self.formats['header_format'])


	# @param ex_json_data if it is None, then it will return 'false' in an if statement
	def _writeRunMetricHeaders(self, ex_json_data=None):
		# RunSheet is where all of the metrics about each run will be written
		self.RunSheet = self.workbook.add_worksheet("QC Metrics")
		self.RunSheet.freeze_panes(1,3)
		
		# First write the QC metrics for each run of each sample.
		# Write the header line. there could definitely be a better way of doing this, but this is what I figured out for now. Just comment and uncomment as needed.
		col = 0
		col += self._writeHeaderCell(col, "Sample #", None, self.RunSheet)
		# check to see if the N or T should be written
		if ex_json_data and 'sample_type' in ex_json_data and ex_json_data['sample_type'] == 'tumor_normal':
			col += self._writeHeaderCell(col, "Normal (N) or Tumor (T)", 7, self.RunSheet)
		col += self._writeHeaderCell(col, "Run #", 5, self.RunSheet)
		col += self._writeHeaderCell(col, "gDNA isolation date", 10, self.RunSheet)
		col += self._writeHeaderCell(col, "gDNA concentration (ng/ul)", 10, self.RunSheet)
		col += self._writeHeaderCell(col, "Library concentration (ng/ul)", None, self.RunSheet)
		col += self._writeHeaderCell(col, "Library prep date", 10, self.RunSheet)
		col += self._writeHeaderCell(col, "Run Date", 12, self.RunSheet)
		col += self._writeHeaderCell(col, "Run ID", None, self.RunSheet)
		col += self._writeHeaderCell(col, "Thermocycler Used", 12, self.RunSheet)
		col += self._writeHeaderCell(col, "Barcode used", 12, self.RunSheet)
		col += self._writeHeaderCell(col, "% Polyclonal", 12, self.RunSheet)
		col += self._writeHeaderCell(col, "Total Aligned Basepairs", 12, self.RunSheet)
		col += self._writeHeaderCell(col, "Total AQ20 Basepairs", 12, self.RunSheet)
		col += self._writeHeaderCell(col, "Mean Read Length", None, self.RunSheet)
		col += self._writeHeaderCell(col, "Median Read Length", None, self.RunSheet)
		# see if we have the exp_read_length
		if ex_json_data and 'cutoffs' in ex_json_data['analysis']['settings'] and 'exp_median_read_length' in ex_json_data['analysis']['settings']['cutoffs']:
			col += self._writeHeaderCell(col, "%% expected read length (out of %d bp)"%ex_json_data['analysis']['settings']['cutoffs']['exp_median_read_length'], 12, self.RunSheet)
		else:
			col += self._writeHeaderCell(col, "% expected read length (out of XXX bp)", 12, self.RunSheet)
		col += self._writeHeaderCell(col, "Median Read Coverage Overall", None, self.RunSheet)
		col += self._writeHeaderCell(col, "% targeted bases >= 30x coverage", 12, self.RunSheet)
		col += self._writeHeaderCell(col, "% amplicons >= 30x coverage", 12, self.RunSheet)
		col += self._writeHeaderCell(col, "% amplicons >= 30x covered at bp +10 (considering fwd/rev read split)", 13, self.RunSheet)
		col += self._writeHeaderCell(col, "% amplicons >= 30x covered at bp n-10 (considering fwd/rev read split)", 13, self.RunSheet)
		col += self._writeHeaderCell(col, "% amplicons >= 30x covered at bp +10 (>=30x) that are also covered at n-10 bp (>=30x)", 13, self.RunSheet)
		#col += self._writeHeaderCell(col, "total number of bases covered at 30x (the # of bases covered in the 'covered_bases region' region.)", 18, self.RunSheet)
		#col += self._writeHeaderCell(col, "% covered bases (n/83046)", 13, self.RunSheet)
		#col += self._writeHeaderCell(col, "% targeted bases (n/84447)", 13, self.RunSheet)
		# check to see if the 'pool_dropout' metrics are available for this spreadsheet
		if not ex_json_data or ('pool_dropout' in ex_json_data['analysis']['settings'] and ex_json_data['analysis']['settings']['pool_dropout'] == True):
			col += self._writeHeaderCell(col, "# of pools <10% median coverage", None, self.RunSheet)
			col += self._writeHeaderCell(col, "# of pools between 10%-50% median coverage ", None, self.RunSheet)
			col += self._writeHeaderCell(col, "# of pools between 50%-75% median coverage", None, self.RunSheet)
			col += self._writeHeaderCell(col, "# of pools passed", None, self.RunSheet)
			col += self._writeHeaderCell(col, "total # of pools", None, self.RunSheet)
		col += self._writeHeaderCell(col, "Ts/Tv", None, self.RunSheet)
		col += self._writeHeaderCell(col, "# Total variants (single allele)", 13, self.RunSheet)
		col += self._writeHeaderCell(col, "# HET variants (single allele rates)", 13, self.RunSheet)
		col += self._writeHeaderCell(col, "# HOM variants (single allele rates)", 13, self.RunSheet)
		col += self._writeHeaderCell(col, "HET/HOM ratio (single allele rates)", 13, self.RunSheet)
		#col += self._writeHeaderCell(col, "3x3 N-N pair (whole amplicon)", 13, self.RunSheet)
		#col += self._writeHeaderCell(col, "Total bases evaluated (>=30x in both runs) (whole amplicon)", 13, self.RunSheet)
		#col += self._writeHeaderCell(col, "% Available Bases (whole amplicon)", 13, self.RunSheet)
		#col += self._writeHeaderCell(col, "3x3 qc observed error counts  (whole amplicon)", 13, self.RunSheet)
		#col += self._writeHeaderCell(col, "3x3 qc error rate  (whole amplicon)", 13, self.RunSheet)
		#col += self._writeHeaderCell(col, "Z-Score error rate (whole amplicon)", 13, self.RunSheet)
		#col += self._writeHeaderCell(col, "Total bases evaluated (>=30x in both runs) (cds only)", 13, self.RunSheet)
		#col += self._writeHeaderCell(col, "% Available Bases (cds only)", 13, self.RunSheet)
		#col += self._writeHeaderCell(col, "3x3 qc observed error counts  (cds only)", 13, self.RunSheet)
		#col += self._writeHeaderCell(col, "3x3 qc error rate  (cds only)", 13, self.RunSheet)
		#col += self._writeHeaderCell(col, "Z-Score error rate (cds only)", 13, self.RunSheet)
		#if ex_json_data and 'cutoffs' in ex_json_data['analysis']['settings']:
		try:
			# print the actual cutoffs used in the script for the header line
			text = "Run Status FAIL if (i) %% expected median read length <%.2f%%, "%(ex_json_data['analysis']['settings']['cutoffs']['perc_exp_median_read_length']) + \
					"OR (ii) %% Amplicons covered at +10 and/or n-10th position <%.2f%%, "%(ex_json_data['analysis']['settings']['cutoffs']['begin_end_amp_cov']) + \
					"OR (iii) pools found <50% median coverage"
			col += self._writeHeaderCell(col, text, 20, self.RunSheet)
			text = "QC status (PASSED RUN status but FAILED with error rate > %.1e in QC table). "%(ex_json_data['analysis']['settings']['cutoffs']['error_rate']) + \
					"Look at only Normal-Normal or Tumor-Tumor comparison in the case of LOH candidates"
			col += self._writeHeaderCell(col, text, 20, self.RunSheet)
		except (KeyError, TypeError):
			col += self._writeHeaderCell(col, "run pass/fail status", None, self.RunSheet)
			col += self._writeHeaderCell(col, "3x3 table error rate pass/fail status", None, self.RunSheet)

		self.RunSheet.set_column(col,col+20,12, self.formats['center'])
		self.RunSheet.set_row(0,100, self.formats['header_format'])


	# @param run_metrics the dictionary containing all of the run_metrics
	def writeRunMetrics(self, run_metrics, ex_json_data=None):
		# first write the headers
		self._writeRunMetricHeaders(ex_json_data)
		self._writeMergedMetricHeaders(ex_json_data)
		merged_row = 0
		run_row = 0
		azure = '_azure'
		
		for sample in sorted(run_metrics):
			# for each sample, change the color
			if len(run_metrics[sample]['runs']) > 0:
				if azure == "":
					azure = "_azure"
				else:
					azure = ""
			for run, metrics in sorted(run_metrics[sample]['runs'].iteritems()):
				if 'json_type' in metrics and metrics['json_type'] == 'merged':
					Sheet = self.Merged
					merged_row += 1
					row = merged_row
				else:
					Sheet = self.RunSheet
					run_row += 1
					row = run_row
				col = 0
				#col += self._check_to_write(row, col, 'sample_num', "" + azure, metrics, Sheet)
				col += self._check_to_write(row, col, 'sample', "" + azure, metrics, Sheet)
				# write the run number with the N or T
				if ex_json_data and 'sample_type' in ex_json_data and ex_json_data['sample_type'] == 'tumor_normal':
					if 'run_type' in metrics and metrics['run_type'] == 'normal':
						Sheet.write(row, col, "N")
					else:
						Sheet.write(row, col, "T")
					col += 1
				# if this is a merged run, then put the merged name for the run number
				if 'json_type' in metrics and metrics['json_type'] == 'merged':
					col += self._check_to_write(row, col, 'run_name', "" + azure, metrics, Sheet)
					col += self._check_to_write(row, col, 'runs_used_to_merge', "" + azure, metrics, Sheet)
				else:
					col += self._check_to_write(row, col, 'run_num', "" + azure, metrics, Sheet)
					col += self._check_to_write(row, col, 'gDNA_isolation', "" + azure, metrics, Sheet)
					col += self._check_to_write(row, col, 'gDNA_conc', "" + azure, metrics, Sheet)
					col += self._check_to_write(row, col, 'lib_conc', "" + azure, metrics, Sheet)
					col += self._check_to_write(row, col, 'lib_prep_date', "" + azure, metrics, Sheet)
					col += self._check_to_write(row, col, 'run_date', "" + azure, metrics, Sheet)
					col += self._check_to_write(row, col, 'run_id', "" + azure, metrics, Sheet)
					col += self._check_to_write(row, col, 'thermocycler', "" + azure, metrics, Sheet)
					col += self._check_to_write(row, col, 'barcode', "" + azure, metrics, Sheet)
					col += self._check_to_write(row, col, 'polyclonality', "perc_format" + azure, metrics, Sheet)
				col += self._check_to_write(row, col, 'aligned_bases', "num_format" + azure, metrics, Sheet)
				col += self._check_to_write(row, col, 'aq20_bases', "num_format" + azure, metrics, Sheet)
				col += self._check_to_write(row, col, 'mean_read_length', "num_format" + azure, metrics, Sheet)
				col += self._check_to_write(row, col, 'median_read_length', "num_format" + azure, metrics, Sheet)
				col += self._check_to_write(row, col, 'perc_exp_median_read_length', 'perc_format' + azure, metrics, Sheet)
				col += self._check_to_write(row, col, 'median_coverage_overall', "num_format" + azure, metrics, Sheet)
				col += self._check_to_write(row, col, 'base_cov', "perc_format" + azure, metrics, Sheet)
				col += self._check_to_write(row, col, 'amp_cov', "perc_format" + azure, metrics, Sheet)
				col += self._check_to_write(row, col, 'begin_amp_cov', 'perc_format' + azure, metrics, Sheet)
				col += self._check_to_write(row, col, 'end_amp_cov', 'perc_format' + azure, metrics, Sheet)
				# give it the dummy 'end_amp_cov' key to write the function of +-10 bp difference. the = is for a function
				col += self._check_to_write(row, col, 'begin_end_amp_cov', 'perc_format' + azure, metrics, Sheet)
				if not ex_json_data or ('pool_dropout' in ex_json_data['analysis']['settings'] and ex_json_data['analysis']['settings']['pool_dropout'] == True) \
						and 'json_type' in metrics and metrics['json_type'] != 'merged':
					col += self._check_to_write(row, col, 'pools_less_than_10', "num_format" + azure, metrics, Sheet)
					col += self._check_to_write(row, col, 'pools_between_10_and_50', "num_format" + azure, metrics, Sheet)
					col += self._check_to_write(row, col, 'pools_between_50_and_75', "num_format" + azure, metrics, Sheet)
					col += self._check_to_write(row, col, 'pools_pass', "num_format" + azure, metrics, Sheet)
					col += self._check_to_write(row, col, 'pools_total', "num_format" + azure, metrics, Sheet)
				col += self._check_to_write(row, col, 'ts_tv', 'dec3_format' + azure, metrics, Sheet)
				col += self._check_to_write(row, col, 'total_vars', 'num_format' + azure, metrics, Sheet)
				col += self._check_to_write(row, col, 'num_het', 'num_format' + azure, metrics, Sheet)
				col += self._check_to_write(row, col, 'num_hom', 'num_format' + azure, metrics, Sheet)
				col += self._check_to_write(row, col, 'het_hom', 'dec3_format' + azure, metrics, Sheet)
				# Now write the pass/fail status
				# let the master spreadsheet handle conditional formatting
				if ex_json_data and 'cutoffs' in ex_json_data['analysis']['settings'] \
						and 'json_type' in metrics and metrics['json_type'] != 'merged':
					col += self._check_to_write(row, col, 'pass_fail_status', azure, metrics, Sheet)
					col += self._check_to_write(row, col, 'pass_fail_3x3_status', azure, metrics, Sheet)
				elif 'json_type' in metrics and metrics['json_type'] == 'merged':
					# if this is a tumor_normal comparison and the runs have been merged and QCd, write this metric
					if ex_json_data and 'sample_type' in ex_json_data and ex_json_data['sample_type'] == 'tumor_normal':
						col += self._check_to_write(row, col, 'merged_perc_avail_bases', 'perc_format' + azure, metrics, Sheet)
						col += self._check_to_write(row, col, 'pass_fail_3x3_status', azure, metrics, Sheet)
					col += self._check_to_write(row, col, 'pass_fail_merged_status', azure, metrics, Sheet)
				# if there are no cutoffs, then the script wont know if it passes or fails yet.
			
				# Set the color of this row according to the current color
				Sheet.set_row(row, None, self.formats[azure])
	
				#run_num = int(metrics['run_num'])
				#runs = int(metrics['runs'])

				# Old naming sheme:
				#col += self._check_to_write(row, col, 'polyclonal', "perc_format" + azure, metrics)
				#col += self._check_to_write(row, col, 'mean', "num_format" + azure, metrics)
				#col += self._check_to_write(row, col, 'median', "num_format" + azure, metrics)
				#col += self._check_to_write(row, col, 'end_amp_cov', '=perc_format' + azure, metrics)
	#			col += self._check_to_write(row, col, 'total_covered', 'num_format' + azure, metrics)
	#			col += self._check_to_write(row, col, 'perc_expected', 'perc_format' + azure, metrics)
	#			col += self._check_to_write(row, col, 'perc_targeted', 'perc_format' + azure, metrics)
				# Old naming scheme:
				#col += self._check_to_write(row, col, 'tstv', 'num_format' + azure, metrics)
		
				# This was never finished. The idea was to choose a representative 3x3 table and show its summary here, but it never really worked out.
#				# Now write the N-N pair and the T-T pairs
#				col += self._check_to_write(row, col, 'same_pair', "" + azure, metrics)
#				col += self._check_to_write(row, col, 'same_total_eligible_bases', 'num_format' + azure, metrics)
#				col += self._check_to_write(row, col, 'same_perc_available_bases', 'perc_format' + azure, metrics)
#				col += self._check_to_write(row, col, 'same_error_count', 'num_format' + azure, metrics)
#				col += self._check_to_write(row, col, 'same_error_rate', 'perc_format' + azure, metrics)
#				try:
#					if 'same_zscore' in metrics and float(metrics['same_zscore']) > 3:
#						self.RunSheet.write_number(row, col, float(metrics['same_zscore']), dec3_format_red)
#						metrics['same_status'] = 'Fail'
#						col += 1
#					else:
#						col += self._check_to_write(row, col, 'same_zscore', 'dec3_format' + azure, metrics)
#				except ValueError:
#					col += self._check_to_write(row, col, 'tn_zscore', '' + azure, metrics)
#				col += self._check_to_write(row, col, 'same_status', "" + azure, metrics)
#		
#			# Now write the T-N pairs
#				col += self._check_to_write(row, col, 'tn_pair', "" + azure, metrics)
#				col += self._check_to_write(row, col, 'tn_total_evaluated', 'num_format' + azure, metrics)
#				col += self._check_to_write(row, col, 'tn_perc_available_bases', 'perc_format' + azure, metrics)
#				col += self._check_to_write(row, col, 'tn_error_count', 'num_format' + azure, metrics)
#				col += self._check_to_write(row, col, 'tn_error_rate', 'perc_format' + azure, metrics)
#				try:
#					if 'tn_zscore' in metrics and float(metrics['tn_zscore']) > 3:
#						self.RunSheet.write_number(row, col, float(metrics['tn_zscore']), dec3_format_red)
#						metrics['tn_status'] = 'Fail'
#						col += 1
#					else:
#						col += self._check_to_write(row, col, 'tn_zscore', 'dec3_format' + azure, metrics)
#				except ValueError:
#					col += self._check_to_write(row, col, 'tn_zscore', '' + azure, metrics)
#				col += self._check_to_write(row, col, 'tn_status', "" + azure, metrics)


	# Write the multiple run Info 3x3 tables if specified
	# @param QC_comparisons the dictionary from all of the _QC.json files
	def write3x3Tables(self, QC_comparisons):
		sheet = ''
		if not self.sheet_per_sample:
			# self.MRsheet is the sheet where all of the 3x3 tables for each of the multiple runs of each sample will be written
			sheet = self.workbook.add_worksheet("3x3 tables")
	
		row = 1
		col = 1
		# QC_comparisons was made by QC_stats.py and has all of the needed 3x3 table stats.
		for sample in sorted(QC_comparisons):
			print "	Writing 3x3 tables for sample: " + sample
			# if specified, add a new sheet per sample
			if self.sheet_per_sample:
				# self.MRsheet is the sheet where all of the 3x3 tables for each of the multiple runs of each sample will be written
				sheet = self.workbook.add_worksheet(sample)
				row = 1
				col = 1
	
			#first make the header for sample.
			sheet.write(row, col, sample, self.formats['header_format2'])
			row += 2
	
			# iterate through the different data sizes (i.e. all, chr1, 718, etc.)
			for data_type in sorted(QC_comparisons[sample]):
				#print 'here with data_type: %s'%data_type
				# comp_type is is either normal_normal, tumor_normal, or tumor_tumor and contains the specified types of comparisons
				# row will be updated as the 3x3 tables are written so that row will always be in teh correct place.
				if "germline_germline" in QC_comparisons[sample][data_type]:
					#sheet.write(row, col, "%s: Normal vs Normal"%data_type, self.formats['header_format2'])
					row = self._write_comparisons(sheet, data_type, QC_comparisons[sample][data_type]['germline_germline'], row, col)
				if "normal_normal" in QC_comparisons[sample][data_type]:
					#sheet.write(row, col, "%s: Normal vs Normal"%data_type, self.formats['header_format2'])
					row = self._write_comparisons(sheet, data_type, QC_comparisons[sample][data_type]['normal_normal'], row, col)
				if "tumor_tumor" in QC_comparisons[sample][data_type]:
					#sheet.write(row, col, "%s: Tumor vs Tumor"%data_type, self.formats['header_format2'])
					row = self._write_comparisons(sheet, data_type, QC_comparisons[sample][data_type]['tumor_tumor'], row, col)
				if "normal_tumor" in QC_comparisons[sample][data_type]:
					#sheet.write(row, col, "%s: Normal vs Tumor"%data_type, self.formats['header_format2'])
					row = self._write_comparisons(sheet, data_type, QC_comparisons[sample][data_type]['normal_tumor'], row, col)

					
	# write the comparisons of a sample. 
	# Will also handle the row, col increases and such
	def _write_comparisons(self, sheet, data_type, qc_comparisons, start_row, start_col):
		row2 = start_row
		highest_row = 0
		col2 = start_col
		# loop through the different comparisons of the given data_type and comparison type.
		for runs_compared, table_values in sorted(qc_comparisons.iteritems()):
			# get the run number
			if 'run1_num' in table_values and 'run2_num' in table_values:
				run1_num = int(table_values['run1_num'])
				run2_num = int(table_values['run2_num'])
			else:
				run1_num = int(runs_compared.split('vs')[0][-1])
				run2_num = int(runs_compared.split('vs')[1][-1])
			# if the run types match, then traingulate by each row being a new run2_num, and each col being a new run1_num
			if table_values['run1_type'] == table_values['run2_type']:
				row2 = start_row + (run2_num-2)*11
				col2 = start_col + (run1_num-1)*6
			# if the run types don't match (i.e. tumor normal), then traingulate by each row being a new normal run(run1_num), and each col being a new tumor run (run2_num)
			else:
				row2 = start_row + (run1_num-1)*11
				col2 = start_col + (run2_num-1)*6
			if row2 > highest_row:
				highest_row = row2
			# now write the 3x3 tables on the given sheet for the given data_type and comparison type
			self._write3x3Table(sheet, data_type, runs_compared, table_values, row2, col2)
		return highest_row+11


	# function to write a 3x3 table to the given sheet
	def _write3x3Table(self, sheet, data_type, runs_compared, table_values, row, col):
		#first make the header for each table.
		header = "%s: %s    -  vs  -    %s"%(data_type, runs_compared.split('vs')[0], runs_compared.split('vs')[1]) # i.e. sample1  Run1  -  vs  - Run2
		sheet.write(row, col, header, self.formats['header_format2'])
	
		# Write the err rate, % available bases and the # of GTs reassigned to the spreadsheet
		header_error_rate = "%s:  %.3e"% ("error rate", float(table_values["error_rate"]))
		row += 1;
		sheet.write(row, col, header_error_rate)
		row += 1
		header_perc_bases = "%s:  %.2f"% ("% available bases", float(table_values["perc_avail_bases"])*100)
		sheet.write(row, col, header_perc_bases)
		row += 1
		if 'reassigned_GTs' in table_values:
			reassigned_GTs = "%s:  %d"% ("# of GTs reassigned", table_values["reassigned_GTs"])
			sheet.write(row, col, reassigned_GTs)
			row += 1;
		
		#now write the table headers
		sheet.write(row+1, col, "WT", self.formats['right'])
		sheet.write(row+2, col, "HET", self.formats['right'])
		sheet.write(row+3, col, "HOM", self.formats['right'])
		sheet.write(row+4, col, "Sum:")
		sheet.write(row, col+1, "WT", self.formats['bottom'])
		sheet.write(row, col+2, "HET", self.formats['bottom'])
		sheet.write(row, col+3, "HOM", self.formats['bottom'])
		sheet.write(row, col+4, "Sum:")
		# write the diagonal cell values:
		sheet.write(row+1, col+1, table_values["WT_WT"])
		sheet.write(row+2, col+2, table_values["HET_HET"])
		sheet.write(row+3, col+3, table_values["HOM_HOM"])
		# now write the off-diagonals
		sheet.write(row+2, col+1, table_values["HET_WT"])
		sheet.write(row+3, col+1, table_values["HOM_WT"])
		sheet.write(row+1, col+2, table_values["WT_HET"])
		sheet.write(row+3, col+2, table_values["HOM_HET"])
		sheet.write(row+1, col+3, table_values["WT_HOM"])
		sheet.write(row+2, col+3, table_values["HET_HOM"])
		# Add the totals for each row and column
		total = int(table_values["WT_WT"]) + int(table_values["WT_HET"]) + int(table_values["WT_HOM"]) + \
				int(table_values["HET_WT"]) + int(table_values["HET_HET"]) + int(table_values["HET_HOM"]) + \
				int(table_values["HOM_WT"]) + int(table_values["HOM_HET"]) + int(table_values["HOM_HOM"])
		# Write the totals for each row and column
		sheet.write(row+4, col+1, str(int(table_values["WT_WT"])+int(table_values["HET_WT"])+int(table_values["HOM_WT"])), self.formats['top'])
		sheet.write(row+4, col+2, str(int(table_values["WT_HET"])+int(table_values["HET_HET"])+int(table_values["HOM_HET"])), self.formats['top'])
		sheet.write(row+4, col+3, str(int(table_values["WT_HOM"])+int(table_values["HET_HOM"])+int(table_values["HOM_HOM"])), self.formats['top'])
		sheet.write(row+1, col+4, str(int(table_values["WT_WT"])+int(table_values["WT_HET"])+int(table_values["WT_HOM"])), self.formats['left'])
		sheet.write(row+2, col+4, str(int(table_values["HET_WT"])+int(table_values["HET_HET"])+int(table_values["HET_HOM"])), self.formats['left'])
		sheet.write(row+3, col+4, str(int(table_values["HOM_WT"])+int(table_values["HOM_HET"])+int(table_values["HOM_HOM"])), self.formats['left'])
		sheet.write(row+4, col+4, str(total))
	

# start here
if (__name__ == "__main__"):
	# parse the arguments
	parser = OptionParser()
	
	# All of the arguments are specified here.
	parser.add_option('-s', '--sample_path', dest='sample_paths', action='append', help='REQUIRED: /path/to/the/sample_dir. Can provide muliple paths')
	parser.add_option('-r', '--run_info_only', dest='run_info_only', action="store_true", default=False, help="Get only the individual run's info. Default is to get both run info and the 3x3 table qc_comparisons")
	parser.add_option('-q', '--qc_info_only', dest='qc_info_only', action="store_true", default=False, help="Get only the 3x3 table QC comparison info. Default is to get both run info and the 3x3 table qc_comparisons")
	parser.add_option('-S', '--sheet_per_sample', dest='sheet_per_sample', action="store_true", default=False, help="Will write a new sheet containing the 3x3 table comparisons per sample.")
	parser.add_option('-o', '--out', dest='out', default='QC.xlsx', help='Specify the output xlsx file [default: %default]')
	parser.add_option('-j', '--ex_json', dest='ex_json', help='An example sample json file containing analysis settings which we can use for this spreadsheet. (with the -r option)')
	
	# Gets all of the command line arguments specified and puts them into the dictionary args
	(options, args) = parser.parse_args()
	
	if not options.sample_paths:
		print "USAGE_ERR: --sample_path is required! Use -h for help"
		parser.print_help()
		sys.exit(8)
	
	# Ensure the user specifies a .xlsx ending if the -o/--output option is used
	if options.out[-5:] != ".xlsx":
		parser.error("-o (--out) output file must end in .xlsx")
	output_file_path = os.path.abspath("/".join(options.out.split("/")[:-1]))
	# if the dir to place the QC file doesn't exist then make it
	if not os.path.isdir(output_file_path):
		os.mkdir(output_file_path)
	
	# setup the xlsx workbook
	workbook = xlsxwriter.Workbook(options.out)
	# add the self.formats and such to the workbook
	xlsx_writer = XLSX_Writer(workbook)
	if options.sheet_per_sample:
		xlsx_writer.sheet_per_sample = True
	else:
		xlsx_writer.sheet_per_sample = False

	runs_json_data = {}
	QC_3x3_json_data = {}
	# Get the data available by finding the json files in the sample_path specified.
	for sample_path in options.sample_paths:
		sample_path = os.path.abspath(sample_path)
		# if the 3x3 tables are not the only thing you want, then get the QC run metrics
		if not options.qc_info_only:
			runs_json_data = dict(runs_json_data.items() + QC_stats.main_runs_only(sample_path).items())
		# if the QC run metrics are not the only thing you want, then get the 3x3 tables
		if not options.run_info_only:
			QC_3x3_json_data = dict(QC_3x3_json_data.items() + QC_stats.main_QC_only(sample_path).items())

	if not options.qc_info_only:
		ex_json_data = None
		if options.ex_json:
			if not os.path.isfile(options.ex_json):
				print "%s is not found. Unable to use it as an example"
			else:
				ex_json_data = json.load(open(options.ex_json))
		print "	Generating the main run QC run metrics Spreadsheet" 
		xlsx_writer.writeRunMetrics(runs_json_data, ex_json_data)
	if not options.run_info_only:
		print "	Generating the 3x3 QC table Spreadsheets" 
		xlsx_writer.write3x3Tables(QC_3x3_json_data)
	
	print "Finished generating the QC table: " + options.out
	xlsx_writer.workbook.close()

