#! /usr/bin/env python2.7

# Goal: Create pretty 3x3 tables comparing mulitple runs, or Tumor/Normal pairs 

import sys
import os
import json

class Var_Writer():
	def __init__(self, xlsx_writer):
		self.xlsx_writer = xlsx_writer

	def writeVarHeaders(self, sheet):
		sheet.freeze_panes(1,4)
		col = 0
		# write the following headers
		col += self.xlsx_writer._writeHeaderCell(col, 'chr', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'pos', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Ref', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Alt', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Normal GT', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Normal AF', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Normal Alt Depth', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Normal Ref Depth', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Tumor GT', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Tumor AF', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Tumor Alt depth', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Tumor Ref depth', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Variant function', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Variant Exonic function', 15, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Variant Gene', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'AA/CDS change', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Transcript Name', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, '1000g2012apr_all', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'snp129Flag', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'snp137Flag', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'cosmic38Flag', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'SIFT_score', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'SIFT_pred', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Polyphen2_HDIV_score', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Polyphen2_HDIV_pred', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Polyphen2_HVAR_score', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'Polyphen2_HVAR_pred', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'GERP++', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'PhyloP', None, sheet)
		col += self.xlsx_writer._writeHeaderCell(col, 'SiPhy', None, sheet)

	def writeVarMetrics(self, CSV, sheet):
		# first write the Headers
		self.writeVarHeaders(sheet)
		# set variables
		row = 1
		col = 0
		azure = '_azure'

		with open(CSV, 'r') as csv:
			header = csv.readline()
			for line in csv:
				line = line.strip().split('\t')
				for item in line:
					sheet.write(row, col, item)
					col += 1
				row += 1
				col = 0
				# center the row
				sheet.set_row(row, None, self.xlsx_writer.formats[''])

