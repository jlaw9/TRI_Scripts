__author__ = 'mattdyer'

from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np
import json

# add labels to a plot
# @param points The plot object
# @param axis The axis object
def autolabel(points, axis, data):
    # attach some text labels
    for i, point in enumerate(points):
        height = float(point.get_height())

        axis.text(point.get_x()+point.get_width()/2., 1.05*height, '%.1f'%data[i],
                ha='center', va='bottom')

#start here when the script is launched
if (__name__ == '__main__'):

    #set up the option parser
    parser = OptionParser()

    #add the options to parse
    parser.add_option('-c', '--coverage', dest='cov', action='append', help='The coverage files')
    parser.add_option('-j', '--json', dest='json', help='The json file')
    parser.add_option('-o', '--output', dest='output', help='The output folder')
    parser.add_option('-b', '--bed', dest='bed', help='The BED file')
    parser.add_option('-p', '--position', dest='position', help='The position offset used in calculating coverage')
    (options, args) = parser.parse_args()

    #load bed file and set mappings to pools
    file = open(options.bed, 'r')
    startToPools = {}
    endToStart = {}
    maxPool = 0;

    for line in file:
        if not line.startswith('#'):
            line = line.rstrip('\n\r')
            tokens = line.split('\t')

            start = int(tokens[1]) + int(options.position)
            end = int(tokens[2]) - (int(options.position) - 1)
            key1 = '%s-%i' % (tokens[0], start)
            key2 = '%s-%i' % (tokens[0], end)

            #store teh end to start map so we have a way to map everything back to one position
            endToStart[key2] = key1

            #now parse the pools
            poolTokens = tokens[5].split('=')
            pools = poolTokens[len(poolTokens) -1].split(',')

            #store the pools

            if not key1 in startToPools:
                startToPools[key1] = []

            for pool in pools:
                startToPools[key1].append(pool)

                #keep track of max
                if int(pool) > maxPool:
                    maxPool = int(pool)

    #close the file
    file.close()
    coverage = {}

    #now read in the coverage files and sum per amplicon
    for fileName in options.cov:
        #print 'Working on %s' % fileName
        file = open(fileName, 'r')


        for line in file:
            line = line.rstrip('\n\r')
            tokens = line.split('\t')

            chromosome = tokens[0]
            position = int(tokens[1])
            count = int(tokens[2])
            key = '%s-%i' % (chromosome, position)

            if key in startToPools:
                #see we have already seen this
                if key in coverage:
                    #print 'adding %i to %s: %i' % (count, key, coverage[key])
                    coverage[key] += count
                    #print 'added to %s: %i' % (key, coverage[key])
                else:
                    coverage[key] = count
            elif key in endToStart:
                #map to start location and then store count
                key = endToStart[key]

                if key in coverage:
                    #print 'adding %i to %s: %i' % (count, key, coverage[key])
                    coverage[key] += count
                    #print 'added to %s: %i' % (key, coverage[key])
                else:
                    coverage[key] = count

        #close the file
        file.close()

    #now we can look through all the amplicons and generate our per stat pools
    pools = {}
    pools['All'] = []

    for start in startToPools:
        if start in coverage:
            count = coverage[start]

            #store the count
            pools['All'].append(count)

            #now we need to add it to the right pool
            for pool in startToPools[start]:
                if not pool in pools:
                    pools[pool] = []

                pools[pool].append(count)

        #else:
            #print 'Warning: %s not found in coverage dictionary' % (start)

    #now compute some stats
    medians = {}

    for pool in pools:
        sortedList = sorted(pools[pool])

        #grab the middle value
        middle = len(sortedList) / 2
        median = sortedList[middle]

        #if total length of items was odd then we need to average the two values in the middle
        if len(sortedList) % 2 == 1:
            median = (median + sortedList[middle + 1]) / 2

        #print '%s: %i' % (pool, median)
        medians[pool] = median

    #calculate how many pools fall into certain buckets
    poolsLessThan10 = 0
    poolsBetween10And50 = 0
    poolsBetween50And75 = 0
    poolsPass = 0
    allMedian = medians['All']

    graphMedians = []
    graphLabels = []
    graphColors = []

    #add the all group first
    graphMedians.append(allMedian)
    graphLabels.append('All')
    graphColors.append('black')

    #used for testing
    #medians['1'] = 10
    #medians['2'] = 50
    #medians['3'] = 150

    for i in range(maxPool):
        pool = '%i' % (i+1)
        if pool in medians:
            median = medians[pool]
            graphMedians.append(median)
            graphLabels.append(pool)

            if(median <= (allMedian * .10)):
                poolsLessThan10 += 1
                graphColors.append('red')
            elif(median > (allMedian * .10) and median <= (allMedian * .50)):
                poolsBetween10And50 += 1
                graphColors.append('orange')
            elif(median > (allMedian * .50) and median <= (allMedian * .75)):
                poolsBetween50And75 += 1
                graphColors.append('yellow')
            else:
                poolsPass += 1
                graphColors.append('green')

    #update the json, first read it in
    jsonData = open(options.json)
    fileData = json.load(jsonData)

    #add the key metrics
    if not 'run_data' in fileData:
        fileData['run_data'] = {}

    fileData['run_data']['pools_total'] = maxPool
    fileData['run_data']['pools_less_than_10'] = poolsLessThan10
    fileData['run_data']['pools_between_10_and_50'] = poolsBetween10And50
    fileData['run_data']['pools_between_50_and_75'] = poolsBetween50And75
    fileData['run_data']['pools_pass'] = poolsPass

    #write new json file
    with open(options.json, 'w') as newJobFile:
            json.dump(fileData, newJobFile, sort_keys=True, indent=4)

    #finally build the plots and spit out some numbers
    fig, ax1 = plt.subplots()

    #set up axes
    xPos = np.arange(len(graphLabels))
    plot1 = ax1.bar(xPos, graphMedians, align='center', alpha=0.4, color='b', label='Median Read Length')
    ax1.set_ylabel('Median Read Coverage')
    ax1.set_xlabel('Amplicon Pool')

    for i in range(len(plot1)):
        plot1[i].set_color(graphColors[i])


    #final additions and then plot
    autolabel(plot1, ax1, graphMedians)
    plt.xticks(xPos, graphLabels)
    plt.savefig('%s/median_read_length_by_pool.png' % (options.output))
