
from optparse import OptionParser
import json
import logging

__author__ = 'mattdyer'

#start here when the script is launched
if (__name__ == '__main__'):
    #set up the option parser
    parser = OptionParser()

    #set up the logging
    logging.basicConfig(level=logging.DEBUG)

    #add the options to parse
    parser.add_option('-j', '--in', dest='json', help='The JSON file')
    parser.add_option('-s', '--status', dest='status', help='The status to update it to')
    (options, args) = parser.parse_args()

    #read in the json file
    jsonData = open(options.json)
    fileData = json.load(jsonData)

    #set the status
    fileData['status'] = options.status

    #dump the file
    with open(options.json, 'w') as newJSONFile:
            json.dump(fileData, newJSONFile, sort_keys=True, indent=4)