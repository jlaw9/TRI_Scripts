#! /usr/bin/env python

from optparse import OptionParser
import os
import os.path
import sys
import re
import json

class Allele_Histo:
	def __init__(self, options):
		self.options = options
	
	def create_histogram(self, matched_var_file):
		with open(matched_var_file, 'r') as var_file:
			normal = []
			tumor = []
			for line in var_file:
				line = line.split('\t')
				normal.append(line[4])
				tumor.append(line[8])
		# now use matplotlib to generate a histogram


if __name__ == '__main__':

	# set up the option parser
	parser = OptionParser()
	
	# add the options to parse
	parser.add_option('-j', '--json', dest='json', help='The json file.')
	parser.add_option('-o', '--output', dest='output', help='The output file. If no output file is specified, output will be written to the screen')
	parser.add_option('-d', '--debug', dest='debug', action='store_true', help='Option to debug or include print statements in code')

	(options, args) = parser.parse_args()

	if not os.path.isfile(options.json):
		print "USAGE_ERROR: json file %s not found"%options.json
		parser.print_help()
		sys.exit(1)

	ex = Allele_Histo(options, reads)

