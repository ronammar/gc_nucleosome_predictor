#!/usr/bin/python

################################################################################
#AUTHOR: Ron Ammar
#DATE: February/8/2012
#MODIFIED: July/24/2012
#LAB: Bader/Giaever/Nislow Labs, Dept of Molecular Genetics
#USAGE: gc_nucelosome_predictor.py [-hg] <Genome FASTA File>
#EXAMPLE: gc_nucleosome_predictor.py -g 138 chr_all.fasta > output.txt
#DESCRIPTION: Uses GC content of a genome to predict nucleosome occupancy. 
#	Outputs the positions to a tab-delimited file with GC-smoothed values based 
#	on GC counts. This is done by binarizing the genome	to 1 where C/G and 0 
#	where A/T and then smoothing these binary values.
################################################################################

import math
import optparse
import re
import sys


##### Parse the command line options
usage= "usage: %prog [-hg] <Genome FASTA File>"
parser= optparse.OptionParser(usage)
parser.add_option("-g", help="Interval length for Gaussian smoother; in bp) [138]", \
	action="store", type="int", dest="gaussian_interval", default=138)
parser.add_option("--no-secondary-smoothing", \
	help="Turn off SMA (Simple Moving Average) to smooth afterwards for less jagged graph", \
	action="store_false", dest="secondary_smooth", default=True)
(options, args)= parser.parse_args()
if len(args) != 1:
	parser.error("incorrect number of arguments")


##### Data Values
previousExtremeCoverage= None


##### Functions ----------------------------------------------------------------

##### Find the Gaussian value with a certain period
#	Returns the value as a floating point number
#	NOTE: the AUC of a Gaussian = 1, so no need to scale it by 1/period.
def findGaussian(middle, coverageList, period, mean):
	middle= int(middle) # ensure that the centre coordinate is cast as an int
	smaStart= middle - (period / 2)
	smaStop= middle + (period / 2)
	
	# Edge cases (where interval falls off the list); Remove all indexes less
	# than 0 or greater/equal to length of the list
	if smaStart < 0:
		smaStart= 0
	if smaStop >= len(coverageList):
		smaStop= len(coverageList) - 1
	
	interval= range(smaStart, smaStop + 1)

	# Apply the Gaussian
	convolution= 0.0
	index= -(period/2)
	s= float(period/6) # standard deviation; division by 6 is key to encompass the Guassian in 3 std each way
	for i in interval:
		# ensure that all numbers are cast as floats, to avoid errors in python's
		# integer arithmetic
		coefficient= \
			(1.0 / (math.sqrt(2.0 * math.pi) * s)) * (math.e ** ((((index - mean)/(2.0 * s)) ** 2.0) * -1.0))
		convolution += coefficient * coverageList[i]

		index += 1 # should end at index == period/2
		
	return convolution


##### Find the simple moving average with a certain period
#	Returns the mean as a floating point number
def findSMA(middle, coverageList, period):
	middle= int(middle) # ensure that the centre coordinate is cast as an int
	smaStart= middle - (period / 2)
	smaStop= middle + (period / 2)
	
	# If the number is even, the middle coordinate cannot be in the middle, so
	# offset the smaStop value
	if period % 2 == 0:
		smaStop -= 1
	
	# Edge cases (where interval falls off the list); Remove all indexes less
	# than 0 or greater/equal to length of the list
	if smaStart < 0:
		smaStart= 0
	if smaStop >= len(coverageList):
		smaStop= len(coverageList) - 1
	
	interval= range(smaStart, smaStop + 1)
	
	# Calculate the mean for this middle coordinate
	sum= 0.0
	for i in interval:
		sum += coverageList[i]
	mean= sum/len(interval)
	
	return mean


##### Iterate along the list and smooth it using simple moving 
#	average smoothing. Store the smoothed data in a list.
#	Method can be either "SMA" or "Gaussian", defaults to Gaussian
#	PRECONDITION: The period (aka interval or bandwidth) must be odd
def smooth(coverageList, period, method="Gaussian"):
	smoothedCoverage= []		
	
	if method == "Gaussian":
		i= 0
		while i != len(coverageList):
			smoothedCoverage += [findGaussian(i, coverageList, period, 0.0)]
			i += 1
	elif method == "SMA":
		i= 0
		while i != len(coverageList):
			smoothedCoverage += [findSMA(i, coverageList, period)]
			i += 1
	
	return smoothedCoverage


##### Smooth binarized GC data and output for each chromosome
def smoothChromosome(chr, gcCalls):
	# Smooth the GC binary data (at any given position- either G/C or A/T) with
	# a Gaussian
	gcCallsSmoothed= smooth(gcCalls, options.gaussian_interval, "Gaussian")

	# Smooth out the little bumps so that extrema identification is not affected
	# if that is performed subsequently
	if options.secondary_smooth:
		gcCallsSmoothed= smooth(gcCallsSmoothed, 15, "SMA")

	# Output; Uses +1 because indeces in python are 0-based
	for index in range(0, len(gcCalls)):
		print "\t".join([chr, str(index + 1), '%.6f' % gcCallsSmoothed[index]])
#-------------------------------------------------------------------------------


##### Store the genome sequence
genomeFile= file(args[0], "r")
currentChromosome= ""
gcCalls= None
for line in genomeFile:
	line= line.strip()

	if re.search(r"^>", line): # in FASTA line
		# If we just completed storing a chromosome, smooth and output that data
		if gcCalls != None:
			smoothChromosome(currentChromosome, gcCalls)

		# Store the next chromosome
		currentChromosome= re.split(r"[>\s]", line)[1]
		gcCalls= []
	else:
		# calculate the GC
		currentGCString= []
		for nucleotide in line:
			if nucleotide == "C" or nucleotide == "c" or nucleotide == "G" or nucleotide == "g":
				currentGCString += [1]
			else:
				currentGCString += [0]
		
		gcCalls += currentGCString

genomeFile.close()

# Output the final chromosomal smoothed data
smoothChromosome(currentChromosome, gcCalls)
