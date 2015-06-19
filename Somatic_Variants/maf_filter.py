#! usr/bin/env python

import sys

with open(sys.argv[1], 'r') as matched_variants:
	header = matched_variants.readline()
	for line in matched_variants:
		lineArr=line.split('\t')
		if lineArr[5] != '.' and lineArr[9] != '.':
			normal_AF = float(lineArr[5])
			tumor_AF = float(lineArr[9])

			if abs(normal_AF - tumor_AF) < 0.03:
				continue
			else:
				print line.strip()
