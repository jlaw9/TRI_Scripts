#! usr/bin/env python

# goal: read csv file containing both normal and tumor depths 
# read also table annovar output 
# append normal and tumor depths to table output from annovar
# keep only canonical transcript 

import sys
import os
import numpy as np
import string
import csv

# required input files: 

inputCSV = open(sys.argv[1], 'r') # let's say somatic.csv 

inputTableAnnovar = open(sys.argv[2], 'r') # table_annovar.pl output using vcf file corresponding to the input CSV (for instance; somatic.vcf)

inputTranscriptFile = sys.argv[3] #canonical transcript List 

outputCSV = open(sys.argv[4], 'w')

# write the header into output CSV

outputCSV.write('\t'. join(['chr', 'pos', 'Ref', 'Alt', 'Normal GT', 'Normal AF', 'Normal Alt Depth', ' Normal Ref Depth', 'Tumor GT', 'Tumor AF', 'Tumor Alt depth', ' Tumor Ref depth', 'Variant function', 'Variant Exonic function', 'Variant Gene', ' AA/CDS change', 'Transcript Name', '1000g2012apr_all', 'snp129Flag', 'snp137Flag', 'cosmic38Flag', 'LJB23_SIFT_score', 'LJB23_SIFT_pred', 'LJB23_Polyphen2_HDIV_score', 'LJB23_Polyphen2_HDIV_pred, ''LJB23_Polyphen2_HVAR_score', 'LJB23_Polyphen2_HVAR_pred', 'LJB23_GERP++', 'LJB23_PhyloP', 'LJB23_SiPhy']) + '\n')



writer = csv.writer(outputCSV, delimiter='\t')

# read transcipt list into an array 

canonicalList = np.loadtxt(inputTranscriptFile, dtype='str')

# need to work on different formatting
#somaticVariants=np.loadtxt(inputCSV,dtype={'names':('chr', 'pos', 'ref', 'alt', 'normalGT', 'normalAF', 'normalAltDepth', 'normalRefdepth','tumorGT', 'tumorAF', 'tumorAltdepth', 'tumorRefdepth'), 'formats': ('S1', 'i4', 'S1', 'S1', 'S1', 'f4', 'i4', 'i4', 'S1', 'f4', 'i4', 'i4')})

# for now read all columns as text 

somaticVariants=np.loadtxt(inputCSV,dtype='str')

# first line is header 
line = inputTableAnnovar.readline()

line_count=0

while line != '':
    
    #print line_count
    line = inputTableAnnovar.readline()
    lineArr=line.split('\t')

    match_found = ''
    if len(line) == 0:
        break
    
    #print len(line)
    #print (lineArr[9])


    if lineArr[5] == "splicing"  :
        AAchanges=lineArr[7].split(',')
        
        for AA in AAchanges:
            AAinfo = AA.split(':')
            transcript = AAinfo[0]
            if transcript in canonicalList:
                match_found = AA
                break
        
        if match_found == '':
            match_found = AAchanges[-1]
                
        CDS_change = match_found.split(':')

                
        AA_report =  CDS_change[2][2:]
        transcript_name = CDS_change[0][0:]


    elif  lineArr[9] !="." and lineArr[9] !="UNKNOWN" :
		#split column 9 by ":"
        AAchanges=lineArr[9].split(',')
		#print (AAchanges)
		#AAno=len(AAchanges)
		#print (AAno)
        
        for AA in AAchanges:
            AAinfo = AA.split(':')
            transcript = AAinfo[1]
            if transcript in canonicalList:
                match_found = AA
                break

        if match_found == '':
            match_found = AAchanges[-1]

        protein_change = match_found.split(':')

        if len(protein_change) > 4:  # sometimes there is no protein change reported for nonframeshift indels in annovar output, need to place "." in thatcase
            AA_report =  protein_change[4][2:]
            transcript_name = protein_change[1][0:]
        else:
            AA_report = '.'
            transcript_name = '.'

    else:
        AA_report = '.'
        transcript_name = '.'

#print (AA_report + '\t')
#print (transcript_name)

#print chr name, chr pos, Ref, Alt etc from input CSV file

    outputCSV.write('\t'. join([somaticVariants[line_count,0], somaticVariants[line_count,1], somaticVariants[line_count,2], somaticVariants[line_count,3], somaticVariants[line_count,4], somaticVariants[line_count,5], somaticVariants[line_count,6], somaticVariants[line_count,7], somaticVariants[line_count,8], somaticVariants[line_count,9], somaticVariants[line_count,10], somaticVariants[line_count,11]]) + '\t')

# print variant function, variant exonic function, gene name from annotation output
    outputCSV.write('\t'. join([lineArr[5], lineArr[8], lineArr[6]]) + '\t')


# note that I am not cross-checking chr names and positions (read from csv file and table annovar output). these should match anyways!
    #print (AA_report + '\t' + transcript_name)

    outputCSV.write('\t'. join([AA_report, transcript_name]) + '\t')

    outputCSV.write('\t'. join([lineArr[10], lineArr[11], lineArr[12], lineArr[13], lineArr[14], lineArr[16], lineArr[17], lineArr[18], lineArr[19], lineArr[20], lineArr[38], lineArr[39], lineArr[40]]))

    line_count=line_count+1


#   break


#if (line_count != len(somaticVariants):
    #print "mismatch observed in variants; check CSV file and table annovar output"
    


