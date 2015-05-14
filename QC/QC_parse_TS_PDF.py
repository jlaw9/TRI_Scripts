#! /usr/bin/env python

# GOAL: Parse the backup TS Run logs to get some of the metrics per barcode for each run.

import sys
import os.path
import os
from subprocess import Popen, PIPE
import subprocess
import re
import json
from optparse import OptionParser

__author__ = 'jefflaw'

# Function to get the % polyclonality
# @param current line
# @return the % polyclonality
def getPolyclonality(line):
	if line[-1] == "%":
		polyclonality = float(line.split(' ')[-1].replace("%","")) # The polyclonality % should be the last item on the line.
	else:
		# Calculate the % polyclonal by adding up all of the reads, and diving the polyclonal reads by the total
		filtered_polyclonal = int(line.split(' ')[-1].replace(",",""))
		total_reads = filtered_polyclonal
		# add up the rest of the reads (there are 3 more sections. Filtered low quality, Filtered primer dimer, and Final Library ISPs)
		for i in range(0,3):
			txtFile.readline()
			line = txtFile.readline().strip()
			total_reads += int(line.split(' ')[-1].replace(",",""))
		polyclonality = (filtered_polyclonal / float(total_reads)) * 100
	return polyclonality

# Function to get the run Date
# @param current line
# @return the run date
def getDate(line):
	dateMatcher = re.compile('.+(\d\d\d\d) (\d\d) (\d\d).+')
	date = dateMatcher.match(line)
	year = date.group(1)
	month = date.group(2)
	day = date.group(3)
	date = "%s/%s/%s"%(month, day, year)
	return date

# Function to fix the bases in barcoded cases where teh sample has extra underscores
# @param line the current line
# @return the total bases for this barcode
def fix_bases(line):
	# example problematic line:
	#IonXpress 078 PNET BC373 12,511,218,3359,898,389,624 75,920,330 164 bp
	bases = ''
	bases_index = 3
	bases_gathered = False
	# the sample could have more than one '_' so keep looping until we pass # the sample
	while not bases_gathered:
		if bases_index >= len(line):
			bases = ''
			break
		# there shouldn't be cases where the number of bases is < 1,000
		if len(line[bases_index].split(',')) < 2:
			bases_index += 1
		# Giga bases should not have more than 4 commas. 
		elif len(line[bases_index].split(',')) > 4:
			problem_index = 0
			fixed = False
			while not fixed:
				if len(line[bases_index].split(',')[problem_index]) > 3:
					bases = ''.join(line[bases_index].split(',')[:problem_index])
					# only get the first three digits of the # problem index
					bases += line[bases_index].split(',')[problem_index][:3]
					fixed = True
					bases_gathered = True
				# heres an extra catch
				elif len(line[bases_index].split(',')[problem_index]) < 1:
					bases = ''
					bases_gathered = True
					break
				else:
					problem_index += 1
		# else this must be the bases. The numbers might not be jumbled together
		else:
			bases = line[bases_index]
			bases_gathered = True
	return bases


#set up the option parser
parser = OptionParser()

#add the options to parse
parser.add_option('-j', '--json', dest='json', help='The JSON file')
parser.add_option('-p', '--pdf', dest='pdf', help='The pdf file')
parser.add_option('-o', '--output_dir', dest='output_dir', help='The location to store the txt version of the pdf')
(options, args) = parser.parse_args()

if not options.json or not options.pdf or not options.output_dir:
	print "All options are required. Use -h for help"
	sys.exit(8)

if not os.path.isfile(options.json):
	# If the json file doesn't exist for some reason, then add these simple metrics to it
	json_path = os.path.abspath(options.json)
	sample = json_path.split("/")[-3]
	run_name = json_path.json.split("/")[-2]
	runData = {'sample': sample, 'name': run_name}
else:
	#read in the json file
	jsonData = open(options.json)
	runData = json.load(jsonData)

pdfFile = options.pdf
# Name of text file after converting pdf to text
pdfText = options.output_dir + "/backupPDF.txt"

# if pdf2txt has already been run for this rotID, then we don't need to rerun it.
if 'run_data' in runData and 'polyclonality' in runData['run_data']:
	print "info has already been gathered. exiting."
	sys.exit()
else:
	# first call pdf2txt.py to change the pdf into text
	command = "pdf2txt.py -L 0.03 -M 30 -F 0.95 %s > %s 2>/dev/null"%(pdfFile, pdfText) 
	returnCode = subprocess.call(command, shell=True)
	
	if returnCode == 0:
		pass
		#print '%s converted to %s successfully!'%(pdfFile, pdfText)
	else:
		print "Could not convert to PDF... check file path given."
		# remove the blank text file and quit.
		os.remove(pdfText)
		sys.exit(1)

# check if this is a barcoded run.
#if 'barcode' in runData:
#	barcoded = True
#else:
#	barcoded = False
command = "grep IonXpress %s >/dev/null 2>&1"%pdfText
returnCode = subprocess.call(command, shell=True)
if returnCode == 0:
	barcoded = True
else:
	barcoded = False


txtFile = open(pdfText, 'r')

date = ''
polyclonality = ''
bases = ''
mean = ''
median = ""
mean = ''
if barcoded:
	# Parse the pdf file as a barcoded run.
	runs = []
	# This is gonna be intense! I need to parse the PDF to get the metrics I want.
	line = txtFile.readline()
	while line != "":
		line = line.strip()
		#if line == "Filtered: Polyconal":
		#if re.search("Polyclonal", line):
		if re.search("Filtered: Polyclonal", line):
			polyclonality = getPolyclonality(line)
	
		
		# -----------------------------------
		# ------- READ BARCODE METRICS -------
		# -----------------------------------
	
		# The barcode information line should look like this:
		#	Barcode		 Sample  Bases		>=Q20	  Reads	  Mean Read Length
		#	IonXpress 001 NOSM 13,796,178 10,531,909 102,594 134 bp
		if re.search(r"^[IonXpress]", line) and re.search(r"[bp]$", line): 
			# Get all of the Barcode info on this page of the pdf
			while re.search("IonXpress", line):
				line = line.split(' ')
				barcode = line[0] + '_' + line[1]
				mean = line[-2]
				sample = line[2]
				bases = line[3]
				# it's not quite this simple... The sample could have an '_'
				# which would throw off the calculation. Also could have  numbers jumbled together... for example:
		#IonXpress 078 PNET BC373 12,511,218,3359,898,389,624 75,920,330 164 bp
				if len(line[3].split(',')) < 2 or len(line[3].split(',')) > 4:
					bases = fix_bases(line)
				runs.append('\t'.join([barcode, sample, bases, mean]))
				line = txtFile.readline()
				line = txtFile.readline().strip()
	
		# Now get the Run Date.
		if re.search(r"R \d\d\d\d \d\d \d\d", line): # There should only be one line with this pattern
			date = getDate(line)
			# Finished gathering all of the metrics needed from the file. break and write info
			break
	
		line = txtFile.readline()
	
else:
	# This is not a barcoded run. The pdf needs to be parsed differently
	txtFile.readline()
	txtFile.readline()
	txtFile.readline()
	txtFile.readline()
	line = txtFile.readline().strip()
	# get the Total Bases from the 5th line of the pdf.
	if re.search("G", line):
		bases = float(line.split("G")[0])
	elif re.search("M", line):
		bases = float(line.split("M")[0]) / 100
	else:
		print "total bases not found?"

	# Read through the rest of the pdf to get the Median and mean read lengths, % polyclonality, and the run date
	median = ""
	while line != "":
		line = line.strip()
		# Get the mean and Median read length
		if re.search("bp", line) and median == "":
			bps = line.replace(" ","").split("bp")
			if len(bps) > 3:
				mean = bps[0]
				median = bps[1]
			else:
				median = bps[0]

		# Get the % polyclonality
		if re.search("Filtered: Polyclonal", line):
			polyclonality = getPolyclonality(line)

		# Now get the Run Date.
		if re.search(r"R \d\d\d\d \d\d \d\d", line): # There should only be one line with this pattern
			date = getDate(line)
			# Finished gathering all of the metrics needed from the file. break and write info
			break

		line = txtFile.readline()


# add the 'run_data' dict to runData if it's not found
if 'run_data' not in runData:
	runData['run_data'] = {}

if barcoded:
	runData['run_data']['barcoded'] = True
	phantom_bases = 0
	total_bases = 0
	num_runs = 0
	# Finished gathering the metrics for this PDF. Now write the metrics to a json for each barcode.
	for run in runs:
		run = run.split("\t")
		try:
			# add up the bases for the run.
			if run[0][0:3] != "Ion":
				# This must be a phantom barocode or something.
				phantom_bases += int(run[2].replace(",",""))
			else:
				num_runs += 1
				total_bases += int(run[2].replace(",",""))
		except ValueError:
			print 'Unable to gather bases'
			phantom_bases = ''
			total_bases = ''
		# if the current barcode matches this runs barcode, then add the run metrics.
		if run[0] == runData['barcode']:
			# add the run metrics to the run_data dicitonary
			runData['run_data']['bases'] = run[2].replace(",","")
			runData['run_data']['mean_read_length'] = run[3]
	# these metrics are used to calculate the percentage of expected bases.
	runData['run_data']['total_chip_bases'] = total_bases - phantom_bases
	runData['run_data']['phantom_bases'] = phantom_bases
	runData['run_data']['num_runs_on_chip'] = num_runs
else:
	runData['run_data']['barcoded'] = False
	runData['run_data']['total_bases'] = bases
	runData['run_data']['median_read_length'] = median 
	# If the mean was available, then also store that.
	try:
		runData['run_data']['mean_read_length'] = mean
	except NameError:
		pass


runData['run_data']['polyclonality'] = polyclonality
runData['run_data']['run_date'] = date


txtFile.close()

#dump the json file
with open(options.json, 'w') as newJSONFile:
	json.dump(runData, newJSONFile, sort_keys=True, indent=4)

