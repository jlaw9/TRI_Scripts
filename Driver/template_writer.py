__author__ = 'mattdyer'

#This class is responsible for writing the different job wrappers use at TRI
class TemplateWriter:

	## The constructor
	# @param self The object pointer
	# @params outputDirectory The output directory
	# @param softwareDirectory The software root directory
	def __init__(self, outputDirectory, softwareDirectory, QCSettings, sendEmail):
		#save the output directories
		self.__outputDirectory = outputDirectory
		self.__softwareDirectory = softwareDirectory
		self.__QCSettings = QCSettings
		self.__sendEmail = sendEmail

	## Write the analysis template
	# @param self The object pointer
	# @param job The json job object
	# @returns A string to the shell script file to be run
	def writeTemplate(self, job):
		#create the output file
		outputFile = '%s/job.sh' % (self.__outputDirectory)
		fileHandle = open(outputFile, 'w')

		#now depending on the analysis time, add the content we need
		if job['analysis']['type'] == 'qc_tvc':
			#we want to merge, QC, then call variants
			self.__writeHeader(job, fileHandle)
			self.__writeStatusChange('running', job, fileHandle, False)
			self.__writeCovTVCTemplate(job, fileHandle)
			self.__writeStatusChange('finished', job, fileHandle, True)
		elif job['analysis']['type'] == 'qc_compare':
			self.__writeHeader(job, fileHandle)
			self.__writeStatusChange('running', job, fileHandle, False)
			self.__writeQCCompareTemplate(job, fileHandle)
			self.__writeStatusChange('finished', job, fileHandle, True)
		elif job['analysis']['type'] == 'qc_sample':
			self.__writeHeader(job, fileHandle)
			self.__writeStatusChange('running', job, fileHandle, False)
			self.__writeQCSampleTemplate(job, fileHandle)
			self.__writeStatusChange('finished', job, fileHandle, True)

		#close the file handle
		fileHandle.close()

		#return the file string
		return(outputFile)

	## Write the header stuff for a job
	# @param self The object pointer
	# @param job The json job object
	# @param file The file handle
	def __writeHeader(self, job, fileHandle):
		#this is just a dummy holder for now
		fileHandle.write('#! /bin/bash\n')
		fileHandle.write('#$ -wd %s\n' % self.__outputDirectory)
		if 'sample' in job:
			fileHandle.write('#$ -N %s.%s.%s\n' % (job['project'], job['sample'], job['name']))
		else:
			fileHandle.write('#$ -N %s.%s.%s\n' % (job['project'], job['sample_name'], job['analysis']['type']))
		fileHandle.write('#$ -V\n')
		fileHandle.write('#$ -S /bin/bash\n\n')

	## Write the code for merging BAM files
	# @param self The object pointer
	# @param job The json job object
	# @param file The file handle
	def __writeMergeTemplate(self, job, fileHandle):
		#create a seperate job for each file in the directory
		print ''

	## Write the code for updating the status in the JSON file
	# @param self The object pointer
	# @param status The json job object
	# @param job The json file to update the status in
	# @param fileHandle The file handle
	# @param wrap Boolean of whether or not to wrap the status
	def __writeStatusChange(self, status, job, fileHandle, wrap):
		#set the status and wrap if requested
		if not wrap:
			fileHandle.write('python %s/Driver/update_json.py -j %s -s %s\n' % (self.__softwareDirectory, job['json_file'], status))
			# I'm not sure how to send emails the normal way so I set up ssmtp on the linux server.
			#fileHandle.write('echo -e "%s beginning analysis" | ssmtp -vvv jlaw@childhooddiseases.org >/dev/null 2>&1\n' % (job['sample_name']))
		else:
			fileHandle.write('if [ $? -ne 0 ]; then\n')
			fileHandle.write('\tpython %s/Driver/update_json.py -j %s -s %s\n' % (self.__softwareDirectory, job['json_file'], "failed"))
			# send an error email automatically if something goes wrong
			if 'emails' in job:
			   for email in job['emails']:
				  fileHandle.write('\tprintf "%s finished with a status of %s. \\n`grep sample_status *.json`\\n" | ssmtp -vvv %s >/dev/null 2>&1\n' % (job['sample_name'], "failed", email))
			#fileHandle.write('\techo "%s finished with a status of %s. `grep sample_status *.json`" | ssmtp -vvv jlaw@childhooddiseases.org >/dev/null 2>&1\n' % (sample_name, "failed"))
			fileHandle.write('else\n')
			fileHandle.write('\tpython %s/Driver/update_json.py -j %s -s %s\n' % (self.__softwareDirectory, job['json_file'], status))
			# QC_sample is going to email the xlsx file when it finishes so don't send an email here.
			#fileHandle.write('\tprintf "%s finished with a status of %s. \\n`grep sample_status *.json`\\n" | ssmtp -vvv jlaw@childhooddiseases.org >/dev/null 2>&1\n' % (sample_name, status))
			#fileHandle.write('\techo "%s finished with a status of %s" | ssmtp -vvv jlaw@childhooddiseases.org >/dev/null 2>&1\n' % (sample_name, status))
			fileHandle.write('fi\n')


	## Write the code for running coverage analysis and tvc
	# @param self The object pointer
	# @param job The json job object
	# @param file The file handle
	def __writeCovTVCTemplate(self, job, fileHandle):
		#default is to not flag dups
		dupFlag = '--remove_dup_flags'

		#see if settings want to mark dups though
		if 'mark_dups' in job['analysis']['settings']:
			#if it is set to true, then we change the flag
			if job['analysis']['settings']['mark_dups'] == 'true':
				print 'here'
				dupFlag = '--flag_dups'

		#default is AmpliSeq for coverage analysis
		coverageAnalysisFlag = '--ampliseq'

		#see if the settings say targetseq
		if 'capture_type' in job['analysis']['settings']:
			#if it is set to true, then we change the flag
			if job['analysis']['settings']['capture_type'].lower() == 'targetseq' or job['analysis']['settings']['capture_type'].lower() == 'target_seq':
				coverageAnalysisFlag = '--targetseq'

		for file in job['analysis']['files']:
			fileHandle.write('bash %s/scripts/runTVC_COV.sh --ptrim PTRIM.bam --cleanup %s %s --cov %s %s --tvc %s %s --output_dir %s %s/%s\n' % (self.__softwareDirectory, dupFlag, coverageAnalysisFlag, job['analysis']['settings']['qc_merged_bed'], job['analysis']['settings']['qc_unmerged_bed'], job['analysis']['settings']['tvc_bed'], job['analysis']['settings']['tvc_parameter_json'], job['output_folder'], job['output_folder'], file))

	## Write the code for running qc comparisons
	# @param self The object pointer
	# @param job The json job object
	# @param file The file handle
	def __writeQCCompareTemplate(self, job, fileHandle):
		#default is to analyze all chromsomes
		chrFlag = '-chr chr1'
		# default is to cleanup
		cleanupFlag='-cl'

		# if analyze_all_chromosomes is set to true in the json file, then take off the chrFlag
		if job['analysis']['settings']['analyze_all_chromosomes'] == True:
			chrFlag = ''

		# If cleanup is specified to false in the json file, then remove the cleanupFlag
		if 'cleanup' in job['analysis']['settings'] and job['analysis']['settings']['cleanup'] == False:
			cleanupFlag = ''

		#let's check the type
		if job['analysis']['settings']['type'] == 'germline':
			#germline sample
			fileHandle.write('bash %s/scripts/QC/QC_sample.sh --beg_bed %s --end_bed %s -s %s -g %s %s %s %s -a %s -b %s %s %s\n' % (self.__softwareDirectory, job['analysis']['settings']['beg_bed'], job['analysis']['settings']['end_bed'], job['sample_folder'], job['analysis']['settings']['tvc_parameter_json'], job['analysis']['settings']['min_base_coverage'], job['analysis']['settings']['wt_cutoff'], job['analysis']['settings']['hom_cutoff'], job['analysis']['settings']['min_amplicon_coverage'], job['analysis']['settings']['project_bed'], chrFlag, cleanupFlag))
		elif job['analysis']['settings']['type'] == 'tumor_normal':
			#tumor_normal sample
			fileHandle.write('bash %s/scripts/QC/QC_sample.sh --beg_bed %s --end_bed %s -s %s -all %s %s %s %s %s %s %s %s -a %s -b %s %s %s\n' % (self.__softwareDirectory, job['analysis']['settings']['beg_bed'], job['analysis']['settings']['end_bed'], job['sample_folder'], job['analysis']['settings']['normal_tvc_json'], job['analysis']['settings']['normal_min_base_coverage'], job['analysis']['settings']['normal_wt_cutoff'], job['analysis']['settings']['normal_hom_cutoff'], job['analysis']['settings']['tumor_tvc_json'], job['analysis']['settings']['tumor_min_base_coverage'], job['analysis']['settings']['tumor_wt_cutoff'], job['analysis']['settings']['tumor_hom_cutoff'], job['analysis']['settings']['min_amplicon_coverage'], job['analysis']['settings']['project_bed'], chrFlag, cleanupFlag))
		#add other types later


	## Write the code for running qc comparisons
	# @param self The object pointer
	# @param job The json job object
	# @param file The file handle
	def __writeQCSampleTemplate(self, job, fileHandle):
		softwareDirectory = "%s/QC"%(self.__softwareDirectory)

		sendEmail = ''
		if self.__sendEmail:
			sendEmail = '-e'

		# the settings should all be contained within the json file.
		fileHandle.write('python2.7 %s/QC_sample.py --json %s %s %s\n'%(softwareDirectory, job['json_file'], self.__QCSettings, sendEmail))
		# use this for if the samples you're running don't have the cutoffs set or need to have teh pass/fail status set
		#fileHandle.write('python %s/QC_sample.py --json %s --update_cutoffs /home/ionadmin/jeff/PNET.json --pass_fail \n'%(softwareDirectory, job['json_file']))
		# use this for creating all 3x3 tables regardless of pass/fail status
		#fileHandle.write('python2.7 %s/QC_sample.py --json %s --qc_all\n'%(softwareDirectory, job['json_file']))
		# use this to recalculate the 3x3 tables (specificaly to update the GT cutoffs.
		#fileHandle.write('python2.7 %s/QC_sample.py --json %s --recalc_3x3_tables -e\n'%(softwareDirectory, job['json_file']))
		# use this to add the aligned stats (specificaly to update the GT cutoffs.
		#fileHandle.write('python2.7 %s/QC_sample.py --json %s --get_alignment_stats\n'%(softwareDirectory, job['json_file']))
		# use this to correct the merged bam name
		#fileHandle.write('python2.7 %s/QC_sample.py --json %s --fix_merged_names\n'%(softwareDirectory, job['json_file']))
		# use this to create the 718 gene list
		#fileHandle.write('python2.7 %s/QC_sample.py --json %s --get_718_subset\n'%(softwareDirectory, job['json_file']))
		# use this to get the overlapping metric
		#fileHandle.write('python2.7 %s/QC_sample.py --json %s --overlap\n'%(softwareDirectory, job['json_file']))
		# compare Normal A_227 merged with Tumor runs.
		#fileHandle.write('python2.7 %s/QC_sample.py --json %s --A_227 -e\n'%(softwareDirectory, job['json_file']))
		# annotate the somatic variants and store them in ~/jeff/Lung_Somatic/[sample]_somatic.vcf. Need to make this option more applicable to other projects
		#fileHandle.write('python2.7 %s/QC_sample.py --json %s --somatic_variants\n'%(softwareDirectory, job['json_file']))
		# run the coverage analysis from 0x - 30x
		#fileHandle.write('python2.7 %s/QC_sample.py --json %s --get_coverages\n'%(softwareDirectory, job['json_file']))
		# run the coverage analysis from 0x - 30x
		#fileHandle.write('python2.7 %s/QC_sample.py --json %s --samtools_gatk -e\n'%(softwareDirectory, job['json_file']))
		#add other types later


	## Write the code for running TVC
	# @param self The object pointer
	# @param job The json job object
	# @param file The file handle
	def __writeTVCTemplate(self, job, fileHandle):
		#this is just a dummy holder for now
		print ''
