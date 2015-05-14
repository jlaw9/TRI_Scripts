#!/usr/bin/env python2.7
#this script is responsible for managing / driving the execution of jobs for Atlas

import glob
import json
import shutil
import os
import logging
import time
import datetime
import fnmatch
from runner import Runner
from template_writer import TemplateWriter
import sys
from optparse import OptionParser

__author__ = 'mattdyer'

class JobManager:
    ## The constructor
    # @param self The object pointer
    # @param outputDirectory The output directory path
    # @param sampleDirectory The sample directory path
    # @param jobFilters The type of jobs we want to launch
    def __init__(self, softwareDirectory, sampleDirectory, jobFilters):
        self.__softwareDirectory = softwareDirectory
        self.__sampleDirectories = sampleDirectory
        self.__jobFilters = jobFilters

    ## Manage the job by parsing the json and finding which analysis to kick off
    # @param self The object pointer
    # @param json The json analysis string pulled from the database
    # @param baseDir The sample base directory to use
    # @returns The SGE job id
    def manageJob(self, jobFile, baseDir):
        #we will want to capture the process exit code and SGE number
        jobNumberSGE = -1

        #load the job file
        jsonData = open(jobFile)
        fileData = json.load(jsonData)
        #logging.debug('%s - %s' % (getTimestamp(), fileData))

        #create the output folder
        if 'name' in fileData:
            outputFolder = '%s/%s' %(fileData['sample_folder'], fileData['name'])
        else:
            outputFolder = fileData['output_folder']
        logging.debug('%s - Creating output folder %s' % (getTimestamp(), outputFolder))
        fileData['output_folder'] = outputFolder

        #see if it exists first
        if not os.path.exists(outputFolder):
            os.makedirs(outputFolder)

        #write the job template
        fileData['json_file'] = jobFile
        templateWriter = TemplateWriter(outputFolder, self.__softwareDirectory)
        analysisFile = templateWriter.writeTemplate(fileData)

        #now we can pass the job to be executed over to the job runner
        runner = Runner("CommandLine")
        logging.info('%s - Starting %s' % (getTimestamp(), analysisFile))
        fileData['status'] = 'submitted'
        fileData['output_folder'] = outputFolder
		# if the --requeue option was specified, update the queue
        if options.requeue:
            fileData['analysis']['settings']['queue'] = options.requeue

        #update the json
        self.__updateJSON(jobFile, fileData)

        #submit the job to SGE
        sgeJobID = runner.submitToSGE('%s/job.sh' % (outputFolder), fileData)
        fileData['status'] = 'queued'
        fileData['sge_job_id'] = sgeJobID
        logging.info('%s - Submitted to SGE (%i)' % (getTimestamp(), sgeJobID))

        #update the json
        self.__updateJSON(jobFile, fileData)

    ## Update the json file
    # @param self The object pointer
    # @param jobFile The input json job file
    # @param job The json job object
    def __updateJSON(self, jobFile, job):
        #now overwrite the old file
        with open(jobFile, 'w') as newJobFile:
            json.dump(job, newJobFile, sort_keys=True, indent=4)

    ## Get the list of pending jobs
    # @param self The object pointer
    # @returns An array of pending job json files
    def getPendingJobs(self):
        #get all the json files that are in in the sample directories
        logging.debug("%s - Looking for JSON job files in %s" % (getTimestamp(), self.__sampleDirectories))
        files = {}
        jobsToProcess = {}
		# instantiate a runner object in case it's needed
        runner = Runner("CommandLine")

        #recurse through and find the json files
        for directory in self.__sampleDirectories:
			# set directory as the absolute path
            directory = os.path.abspath(directory)
            #see if this directory exists
            if os.path.isdir(directory):
                for root, dirnames, filenames in os.walk(directory):
                    for filename in fnmatch.filter(filenames, '*.json'):
                        #see if it is the right type
                        #logging.debug('%s - Looking at %s' % (getTimestamp(), os.path.join(root, filename)))
                        jsonData = open(os.path.join(root, filename))
                        try:
                            fileData = json.load(jsonData)
                        except ValueError:
                            print "\nERROR: The JSON file: %s could not be loaded by python as a JSON file.\n"%os.path.join(root, filename) + \
                                  "Please ensure that the JSON file is not currently open in vim and that it is formatted correctly. Otherwise, delete it."
                            sys.exit(1)

                    
                        #since other json files may be around, let's be sure they have the analysis type flag
                        #can use this to filter things too
                        if 'analysis' in fileData and 'type' in fileData['analysis'] and 'status' in fileData:
							# If the analysis type matches the given jobFilter, check to see if this job should be started
                            if fileData['analysis']['type'] in self.__jobFilters:

								# if the status is 'queued' and the user specified requeue, then delete the current pending or running job
                                if 'sge_job_id' in fileData and (fileData['status'] == 'queued' and options.requeue):
                                    logging.info('%s - Deleting Job ID %s' % (getTimestamp(), fileData['sge_job_id']))
									# delete this json file's job id
                                    command = "qdel %s"%fileData['sge_job_id']                           
                                    runner.runCommandLine(command) 
									# reset the fileData status to pending to requeue this job
                                    fileData['status'] = 'pending'

								# If the status is 'pending' or if the user specified a status type to rerun, then start this job
                                if fileData['status'] == 'pending' or fileData['status'] == options.rerun:
                                        #job was the right type so we can add to list
                                        files[os.path.join(root, filename)] = directory
                                    
                                    
        #process each of the json files
        for file in files:
            #rename the json file so it won't get picked up by the next thread
            logging.debug('%s - Found %s' % (getTimestamp(), file))
#           shutil.move(file, '%s' % file)
            #shutil.copy(file, '%s_read' % (file))
                                    
            #add the file to the array
            jobsToProcess['%s' % (file)] = files[file]

        #return the array
        return(jobsToProcess)

## simple method to get a time stamp
# @returns A timestamp
def getTimestamp():
    #get the time / timestamp
    currentTime = time.time()
    timeStamp = datetime.datetime.fromtimestamp(currentTime).strftime('%Y-%m-%d_%H-%M-%S')

    #return the stamp
    return(timeStamp)

#start here when the script is launched
if (__name__ == "__main__"):
	# setup the option parser
	parser = OptionParser()

	# add the arguments
	parser.add_option('-s', '--sample_dir', dest='sample_dir', action='append', help='Recurse through the given directory and look for *.json jobs. default=["/rawdata/projects", "/results/projects", "/mnt/Despina/projects"]')
	parser.add_option('-S', '--software_dir', dest='software_dir', default="/rawdata/legos", help='The software root directory. [default: %default]')
	parser.add_option('-j', '--job_filters', dest="job_filters", action="append", default=["qc_sample"], help="The type of job to start. [default: %default] Choices: 'qc_tvc', 'qc_compare', 'qc_sample'. Put a -j before each job type if running multiple jobs.")
	parser.add_option('-r', '--requeue', dest="requeue", help="Will qdel jobs that have a status of 'queued' and re-submit them with the new specified queue, as well as start 'pending' jobs with the specified queue")
	parser.add_option('-R', '--rerun', dest="rerun", help="Will 'rerun' jobs of the specified status type")
	
	#parse the arguments
	(options, args) = parser.parse_args()

	# put the default here because otherwise the tri_driver will look in all of these places all of the time.
	if not options.sample_dir:
		#options.sample_dir = ["/rawdata/projects", "/results/projects", "/Volumes/HD/mattdyer/Desktop/temp"]
		options.sample_dir = ["/rawdata/projects", "/results/projects", "/mnt/Despina/projects"]

	## Check the job_filters
	#if not options.job_filters:
	#	options.sample_dir = ["/rawdata/projects", "/results/projects", "/mnt/Despina/projects"]
	#	#print "USAGE ERROR: -j (--job_filters) is required. Available choices: 'qc_tvc', 'qc_compare', 'qc_sample'. Use -h for help."
	#	#parser.print_help()
	#	sys.exit(8)

	for job_filter in options.job_filters:
		if job_filter != "qc_tvc" and job_filter != "qc_compare" and job_filter != 'qc_sample':
			print "USAGE ERROR: %s is not an available job filter. Use -h for help"%job_filter
			parser.print_help()
			sys.exit(8)

	#set up the logging
	logging.basicConfig(level=logging.DEBUG)

	#create the job manager
	jobManager = JobManager(options.software_dir, options.sample_dir, options.job_filters)

	#find the pending jobs
	jobsToProcess = jobManager.getPendingJobs()
	logging.info('%s - Found %i analyses to process' % (getTimestamp(), len(jobsToProcess)))

	#process the jobs now
	for job in jobsToProcess:
		logging.info('%s - Working on %s' % (getTimestamp(), job))

		#process the job
		jobManager.manageJob(job, jobsToProcess[job])

		#send an email?

