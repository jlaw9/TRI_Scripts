#!/usr/bin/env python

from optparse import OptionParser
import os
import os.path
import sys
import re
import json

class Subset_Depths:
	def __init__(self, bed, depths, out):
		self.options = options
		self.depths = open(depths, 'r')
		self.bed = open(bed, 'r')
		self.out = open(out, 'w')
	
		self.get_depths()
	
		self.depths.close()
		self.bed.close()
		self.out.close()

	def get_depths(self):
		bed_arr = self.bed.readline().split('\t')
		# because the bed_arr is a subset of the depths file, 
		# we can always assume the bed file is behind the depths file
		for depth_line in self.depths:
			depth_arr = depth_line.split('\t')
			# match the chromosomes
			#while len(bed_arr) > 1 and depth_arr[0] != bed_arr[0] and int(depth_arr[1]) < int(bed_arr[1]):
			#while len(bed_arr) > 1 and depth_arr[0] != bed_arr[0]:
			#	print depth_arr, bed_arr[0], bed_arr[1]
			#	bed_arr = self.bed.readline().split('\t')
			while len(bed_arr) > 1 and depth_arr[0] == bed_arr[0] and int(depth_arr[1]) > int(bed_arr[2]):
				#print depth_arr, bed_arr[0], bed_arr[1]
				bed_arr = self.bed.readline().split('\t')
			# the position in the depth file has to be inside the positions of the bed file
			if len(bed_arr) > 1 and int(depth_arr[1]) > int(bed_arr[1]) and int(depth_arr[1]) <= int(bed_arr[2]):
				self.out.write(depth_line)
				
if __name__ == '__main__':

	# set up the option parser
	parser = OptionParser()
	
	# add the options to parse
	parser.add_option('-b', '--bed', dest='bed', help='The subset bed file')
	parser.add_option('-d', '--depth', dest='depth', help='The depth file created by samtools depth')
	parser.add_option('-o', '--output', dest='output', help='The output file. If no output file is specified, output will be written to the screen')

	(options, args) = parser.parse_args()

	if not options.bed or not options.depth or not options. output:
		print "USAGE_ERROR: all options requred"
		parser.print_help()
		sys.exit(1)

	Subset_Depths(options.bed, options.depth, options.output)

