#! usr/bin/env python

# OY 10/28/14
# Goal: apply grep -E "WT.+HET" matched_variants.csv > somatic.csv  and provide the ouput from grep statement as first argument
# adjust GT settings in grep accrodingly
# provide as 2nd argument the tumor vcf file
# subset out the variants matched in tumor vcf

#example: python getVCF_v2.py path_to_somatic.csv path_toVCF2_Final.vcf  path_to_output.vcf

import sys
import os
import numpy as np
import string


inputFileCSV = open(sys.argv[1], 'r') # matched_variants.csv (either before or after GT update)
#inputFileVCF = open(sys.argv[2], 'r') # tumor vcf

pathToVCF=sys.argv[2]

print (pathToVCF)

#normal_GT=sys.argv[3]
#tumor_GT=sys.argv[4]


#outVCF = open(sys.argv[5], 'w')    # path to output vcf

file_handle = file(sys.argv[3], 'w')


tumorData = np.loadtxt(pathToVCF, dtype='str', skiprows=98) # NOTE: assuming that the number of header lines will always be 98, if this format changes will need to update!!!!


chrNo=np.array(tumorData[:,1],dtype='|S10') # need the chr positions
tumor_chrNo= chrNo.astype(np.int) # convert chr positions to integer

# print DataIn[[0,1,3], :] # printing certain rows to check the data loaded

# read first line. no header present when you use "grep" to generate input matched csv file
line = inputFileCSV.readline()

lineArr=line.split('\t')
normal_GT=lineArr[4]
tumor_GT=lineArr[8]

print (normal_GT)
print (tumor_GT)

# Loop through each line
while line != "":
    lineArr=line.split('\t')
    
    # if genotype for either normal or tumor is ".", skip this line
    # if genotypes assigned to normal and tumor do not match user input, skip this line
    
    if  lineArr[4] !="." and lineArr[8] !="." :
        
        
       # locate the line in tumor vcf for this chr position
       # already saved the chr positions, so will do a match on the position array
       
       matchResult=np.where(tumor_chrNo == int(lineArr[1]))

       if matchResult[0].size == 1 :
           indexVal=matchResult[0]
           #print (indexVal)
           #print (tumorData[indexVal,1])
           #print ( string.join(tumorData[indexVal,1])," ")
           #print (tumorData[indexVal,1])
           #print (tumorData[indexVal,7])
           
           #outVCF.write('\t'.join(''.join([tumorData[indexVal,0]), ''.join(tumorData[indexVal,1]), ''.join(tumorData[indexVal,2]), ''.join(tumorData[indexVal,3]), ''.join(tumorData[indexVal,4]), ''.join(tumorData[indexVal,5]), ''.join(tumorData[indexVal,6]), string.join((tumorData[indexVal,7]), " "), ':'.join(tumorData[indexVal,8], ':'.join(tumorData[indexVal,9]])))
                                           
                                           #outVCF.write('\t'.join([tumorData[indexVal,1].tolist(), tumorData[indexVal,7].tolist()]))
           np.savetxt(file_handle,tumorData[indexVal,:], fmt='%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' , delimiter='\t', newline='\n')


    line = inputFileCSV.readline()
    #break


print 'Finished'
inputFileCSV.close()
file_handle.close()
