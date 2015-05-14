#! /usr/bin/env python

# Goal: Remove duplicate variants, the variants that have Multi Allelic calls, and the variants where FRO + FAO < 30 (or a specified input) from a merged vcf file.


import sys
import re
import os
import subprocess
import json
from tools import *
from optparse import OptionParser

class QC_Filter:
	def __init__(self, options):
		self.options = options
		# input files for the program.
		try:
			# first load the files
			self.vcf = open(options.vcf, 'r') # input file1 vcf 
			self.filtered_vcf = open(options.filtered_vcf, 'w') # Output vcf that has the variants filtered.
			if options.vars_filtered:
				self.vars_filtered = open(options.vars_filtered, 'w')
		except IOError, ValueError:
			print "ERROR: Check that file's exist and that the json file is formatted correctly"
			sys.exit(1)
	
		# set the filtered variables 	
		self.last_chr_pos = "" # This variable will be set to the chr, and pos of the line immediately before the current line.
		self.duplicate_variants_removed = 0
		self.multiple_variants_removed = 0
		self.coverage_variants_removed = 0

		# go go go!
		self.main()
		self.vcf.close()
		self.filtered_vcf.close()
		if options.vars_filtered:
			self.vars_filtered.close()

	# filter the VCF file.
	def main(self):
		for line in self.vcf:
			write = False
			if line[0] == "#": # Write the header lines first
				self.filtered_vcf.write(line)
				if self.options.vars_filtered:
					self.vars_filtered.write(line)
			else:
				line_arr = line.strip().split('\t')
				# if the varian passes all of these filters, then keep it
				if self.duplicate_filter(line_arr) and self.multiple_alleles_filter(line_arr, self.options.multi_allele_max_freq) and self.coverage_filter(line_arr, self.options.min_base_cov):
					self.filtered_vcf.write(line)
				# otherwise throw it away
				elif self.options.vars_filtered:
					self.vars_filtered.write(line)
	
		# now store the metrics in the json file
		if self.options.json:
			json_data = json.load(open(self.options.json))
			if 'run_data' not in json_data:
				json_data['run_data'] = {}
			json_data['run_data']['vars_removed_dup'] = self.duplicate_variants_removed
			json_data['run_data']['vars_removed_multi'] = self.multiple_variants_removed
			json_data['run_data']['vars_removed_cov'] = self.coverage_variants_removed
			# write the json file
			write_json(self.options.json, json_data)
	
		# print the results
		print "# of duplicate variants removed: " + str(self.duplicate_variants_removed)
		print "# of multiple variants removed: " + str(self.multiple_variants_removed)
		print "# of variants removed that had low coverage: %s"%(self.coverage_variants_removed)

	# @param line_arr Check to see if this variant is a duplicate
	# @returns true or false
	def duplicate_filter(self, line_arr):
		chr_pos = line_arr[0] + "_" + line_arr[1]
		# If the current chr and pos are equal to the chr and pos of the line immediately before this line, this is a duplicate entry. This line will be skipped.
		if chr_pos != self.last_chr_pos:
			self.last_chr_pos = chr_pos
			return True
		else:
			self.duplicate_variants_removed += 1
			return False

	# @param line_arr Check to see if this variant has > 2 alleles
	# @param multi_allele_max_freq The maximum allele frequency for multiple alleles (i.e. if > 2 alleles have a frequency > specified cutoff, the variant will be removed)
	# @returns true or false
	def multiple_alleles_filter(self, line_arr, multi_allele_max_freq=.05):
		# Don't keep variants that have a Multi allelic call. We are diploid and should therefore only have two alleles	
		if not re.search(",", line_arr[4]): 
			return True 
		else:
			self.multiple_variants_removed += 1
			return False

	# @param line_arr Check to see if this variant has enough coverage
	# @param min_base_cov A list where each value is the minimum base coverage for each sample listed in the vcf.
	# @returns true or false
	def coverage_filter(self, line_arr, min_base_cov):
		passes = True
		# loop through the info columns. The coverages passed in by the min_base_cov option should correspond to the columns or samples in the vcf file.
		for i in range(0,len(line_arr)-9):
			# Creates a dictionary with the description as the key, and the actual value as the value.	
			vcfInfo = dict(zip(line_arr[8].split(":"), line_arr[9+i].split(":")))   
			# Returns the FRO + FAO values.
			depth = self.get_depth(vcfInfo)
			if depth != None and depth < min_base_cov[i]:
				passes = False
		if not passes:
			self.coverage_variants_removed += 1
		return passes
				
	# @param Takes in a dictionary with the key as the description of a variant, and the value as the value of that description.
	# @returns Returns the amount of variants that were filtered (depth was not > cutoff in either vcf file).
	def get_depth(self, vcfInfo):
		alt_depth = ''
		ref_depth = ''
		try: 
			# Ion Reporter says the following about using FAO and FRO rather than AO and RO:
			# FAO is usually equal to AO; however, due to complex alleles and/or downsampling*, FAO may differ from AO.
			# AF=FAO/(FAO+FRO) and not FAO/FDP. This is because FDP may include reads that don't fit the flow space profile of any hypothesis; in such cases, FDP>=FAO+FRO and this is not used in allele frequency calculation.
			# Exception: When flow correction is not performed and there are no F tags in the VCF file, then DP=AO+RO and AF=AO/DP.
			if 'FAO' in vcfInfo and vcfInfo['FAO'] != '.':   # Ozlem said that the frequency depth score is better than the regular depth score, so use that if we can.
				alt_depth = vcfInfo['FAO']
				ref_depth = vcfInfo['FRO']
			elif 'AO' in vcfInfo:
				alt_depth = vcfInfo['AO']
				ref_depth = vcfInfo['RO']
			# if this is a merged VCF, the variants that were found in one run but not the other will be labelled at a '.' We want to keep them for the hotspot refilling
			elif vcfInfo['GT'] == '.':
				return
			total_depth = int(alt_depth) + int(ref_depth)
			return total_depth
		except KeyError:
			print "KEY ERROR: in the vcf file. Leaving this variant."
			return

# ----------------------------------------------------
# ------------- PROGRAM STARTS HERE ------------------
# ----------------------------------------------------

if __name__ == "__main__":
	
	# set up the option parser
	parser = OptionParser()
	
	# add the options to parse
	parser.add_option('-v', '--vcf', dest='vcf', help='A VCF file to be filtered.')
	parser.add_option('-j', '--json', dest='json', help='The run json file used to store the information about how may variants were filtered')
	parser.add_option('-o', '--filtered_vcf', dest='filtered_vcf', help='Output filtered vcf file')
	parser.add_option('-O', '--variants_filtered_out', dest='vars_filtered', help='Optiont to store the variants that were filtered in another VCF file.')
	parser.add_option('-c', '--min_base_cov', dest='min_base_cov', action='append', type="int", help='The minimum coverage at each variant required. If this is a merged file, use this option twice (-c cutoff_for_run1 -c cutoff_for_run2)')
	parser.add_option('-m', '--multi_allele_max_freq', dest='multi_allele_max_freq', action='append', type="float", help='The maximum allele frequency for multiple alleles (i.e. if > 2 alleles have a frequency > specified cutoff, the variant will be removed)') 
	
	(options, args) = parser.parse_args()
	
	#check to make sure either ID or name was provided
	if(not options.vcf or not options.filtered_vcf or not options.min_base_cov): 
		sys.stderr.write("USAGE ERROR: -v, -o, and -c are all required\n")
		parser.print_help()
		#print "use -h for help"
		sys.exit(8)

	QC_Filter(options)
