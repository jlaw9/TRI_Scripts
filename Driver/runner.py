## @package Runner
#
# This class is responsible for running jobs on the command line

import subprocess
import os
import time

__author__ = 'mattdyer'

class Runner:
    ## The constructor
    # @param self The object pointer
    # @param mode The mode in which jobs will be run (CommandLine or SGE)
    def __init__(self, mode):
        self.__mode = mode

		# emails are sent using template.py and QC_sample.py
        # these are the default ones
        #self.__emails = []
        self.__queue = 'all.q'
        self.__priority = '0'

    ## Run a job on the command-line
    # @param self The object pointer
    # @param systemCall The system call to executed
    # @returns The exit status and the output
    def runCommandLine(self, systemCall):
        #run the call and return the status
        returnCode = subprocess.call(systemCall, shell=True)
        return(returnCode)

    ## Run a job on via SGE
    # @param self The object pointer
    # @param file The the shell script to run
    # @param fileData The dictionary of file data
    # @returns The SGE job number
    def submitToSGE(self, file, fileData):
        #see if there were other emails address provided
        #if 'emails' in fileData:
        #    self.__emails = self.__emails + fileData['emails']

        #see if another queue has been provided
        if 'queue' in fileData['analysis']['settings']:
            self.__queue = fileData['analysis']['settings']['queue']

        #see if a different priority has been give
        if 'priority' in fileData['analysis']['settings']:
            self.__priority = fileData['analysis']['settings']['priority']

        #build the command call
        systemCall = 'qsub -p %s -q %s -e %s/sge.log -o %s/sge.log %s' % (self.__priority, self.__queue, fileData['output_folder'], fileData['output_folder'], file)
        #systemCall = 'qsub -m e -M \'%s\' -p %s -q %s -e %s/sge.log -o %s/sge.log %s' % (','.join(self.__emails), self.__priority, self.__queue, fileData['output_folder'], fileData['output_folder'], file)

        #submit to SGE and grab the job id
        #print systemCall
        jobString = subprocess.check_output(systemCall, shell=True)
        tokens = jobString.split(' ')
        return(int(tokens[2]))
        #return(1)
