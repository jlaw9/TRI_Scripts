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
				if self.duplicate_filter(line_arr) and self.coverage_allele_Filters(line_arr, self.options.min_base_cov, self.options.multi_allele_max_freq):
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
	# @param min_base_cov A list where each value is the minimum base coverage for each sample listed in the vcf.
	# @param multi_allele_max_freq The maximum allele frequency for multiple alleles (i.e. if > 2 alleles have a frequency > specified cutoff, the variant will be removed). I went with a more conservative approach just because.
	# @returns true or false
	def coverage_allele_Filters(self, line_arr, min_base_cov, multi_allele_max_freq):
		# loop through the info columns. The coverages passed in by the min_base_cov option should correspond to the columns or samples in the vcf file.
		for i in range(0,len(line_arr)-9):
			vcfInfo = dict(zip(line_arr[8].split(":"), line_arr[9+i].split(":")))   

			# if this is a merged VCF, the variants that were found in one run but not the other will be labelled at a '.' We want to keep them for the hotspot refilling
			if vcfInfo['GT'] == '.':
				return True

			try:
				# We are diploid and should therefore only have two alleles
				if 'FAO' not in vcfInfo:
					orig_alleles = vcfInfo['AO'].split(',') + [vcfInfo['RO']]
				else:
					orig_alleles = vcfInfo['FAO'].split(',') + [vcfInfo['FRO']]
				# remove the alleles for which we have no GT information. Those will be filled in with a hotspot later.
				alleles = filter(lambda x: x != '.', orig_alleles) 
				total_cov = reduce(lambda x,y: int(x)+int(y), alleles)
				valid_alleles = filter(lambda x : (float(x)/total_cov) > multi_allele_max_freq[i], alleles)
				total_valid_cov = reduce(lambda x,y: int(x)+int(y), valid_alleles)
				# need to implement a new filter so that the alleles that are representative of the person are kept
				# even if those are two alternate alleles. For now, just stick with what we had
				# Would have to update QC_Compare_VCFs. IF both alleles are alternates, should that cell be marked as HOM alt? If only one alternate allele is real, should the other alternate allele be removed?
				#if len(valid_alleles) > 2: 
				if len(orig_alleles) > 2: 
					self.multiple_variants_removed += 1
					return False
				elif total_valid_cov < min_base_cov[i]:
					self.coverage_variants_removed += 1
					return False
				else:
					return True
			except KeyError:
				print "KEY ERROR: in the vcf file. Leaving this variant."
				return True
			

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

	# set the default value if it's not specified
	if not options.multi_allele_max_freq:
		options.multi_allele_max_freq = [0.05]*len(options.min_base_cov)

	QC_Filter(options)
