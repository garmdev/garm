# ----------------------------------------------------------------------------
# v3.0 - 2011June27 -- Assumes files in grid*.csv, but for new version of
#	CDPOP v1.0. Just Dps.
# v2.0 - 2008Dec22 -- Assumes files are in the gridNNN format, but reads in 
#	*.csv, not grid*.csv
# v1.0 - March 2008 -- Assumes files are named grid*.csv
# GeneticDistance.py
# Author: Erin L Landguth
# Created: March 2008 
# Description: This program calculates a genetic distance matrix using:
#	Bray-Curtis Method:
#	formula 1 - 2W/(A+B), where W is the minimum value between the two comp-
#	arison's A and B.  The specific application is to calculate this distance
#	matrix for genotypes of a population with n individuals: Genetic Distance.
#	Proportion of Shared Alleles:
#	Nei's:
#	1 - sqrt(ithfreq*jthfreq)/loci
#	Proportion of Shared alleles: 
#	1 - proportion of shared alleles between individuals.	
# Program Input: directory of *.csv files
# Program Output: oldname+Gdmatrix.csv 
# Program Steps:
#	1. User input information.
#	2. fileList all of the *.csv in directory
#	3. Run genetic distance method
# ----------------------------------------------------------------------------

import glob							# The power of glob
from numpy import *					# General commands and functions
from numpy.random import *			# Random/statistic calculations
import time, datetime,pdb			# A time clock 

# Timing events: start
start_time = datetime.datetime.now()

# ------------------------------------------
# Step 1: Get user information
# ------------------------------------------

# Store directory path name
directory = 'C:/Users/brian/Documents/SchoolWork/PhDWork/PaperDrafts/BirdsProject/Data/St/500GensOnly/1318760945batchrun0mcrun0_inputvariables_50000invsqNoOverlapSt_tree_100_1_10.asc_St_XY.csv.cdmatrix.csv/'
loci = 12
noalleperlocus = 20
alleles = int(noalleperlocus)*ones(loci,int)
#alleles = array([6,10])
nogrids = 1655

# -------------------------------	
# Step 2: List files in directory
# -------------------------------
csvfileList = glob.glob(directory+'grid*.csv')
		
# Get length of files
nocsvfiles = len(csvfileList)

# -----------------------------------	
# Step 3: Run genetic distance method
# -----------------------------------
	
# ------------ Genetic Distance Matrix: Proportion of shared alleles -----------------------
# List all files with .csv extensions (but only the grid ones)

# Debugging break...
#pdb.set_trace()

print '\n'
print 'Creating the proportion of shared alleles genetic distance matrices...'
#pdb.set_trace()
# Get the first globbed file read in
for i in range(nocsvfiles):
						
	# Open file for reading
	inputfile = open(csvfileList[i],'r')
		
	# Read lines from the file
	lines = inputfile.readlines()
			
	#Close the file
	inputfile.close()
	
	# Create an empty matrix to append to
	x = []
	
	# Split up each line in file and append to empty matrix, x
	for l in lines:
		thisline = l.split(',')
		x.append(thisline)
			
	# Store genetic information: genes[bear], but need them as float values
	genes = []
	tempgenes = []
	for k in range(len(x)-1):
		# Get list from read in file
		tempgenes.append(x[k+1][7:int(7+sum(alleles))])
		# Create spot in genes
		genes.append([])
		for j in range(sum(alleles)):
			# Make each list spot an integer
			genes[k].append(float(tempgenes[k][j]))
	#pdb.set_trace()
	# Create a matrix of zeros to be filled
	gendmatrix = zeros((nogrids,nogrids),float)
	#pdb.set_trace()
	# Loop through each individual k
	for k in range(nogrids):
		# Compare individual k to every other inidividual j
		for j in range(nogrids):
			# Create a tempvariable to be written over for each comparison
			tempmin=[]
			# Loop through each allele value
			for alle in range(sum(alleles)):
				# Find the shared alleles between k and j checking the 4 conditions
				if genes[k][alle]==2.0:
					if genes[j][alle]==2.0:
						tempmin.append(2)
					elif genes[j][alle]==1.0:
						tempmin.append(1)
				elif genes[k][alle]==1.0:
					if genes[j][alle]==1.0:
						tempmin.append(1)
					elif genes[j][alle]==2.0:
						tempmin.append(1)
			# Write the Dps value to gendmatrix
			gendmatrix[k][j] = 1-float(sum(tempmin))/(2*loci)
	
	# Strip directory/filename of grid and add 'Gdmatrix.csv'
	gdpathname = csvfileList[i].replace('grid','Gdmatrix')
	
	# Create file to write matrix to
	outputfile = open(gdpathname,'w')
	
	# Sequence each row in the matrix
	for seqrow in gendmatrix:
	
		# Grab each element in each row and write element to outputfile
		for ele in range(len(seqrow)):
			outputfile.write(str(seqrow[ele]))
			# Add comma
			outputfile.write(',')
		
		# Return line
		outputfile.write('\n')
	
	# Close file
	outputfile.close()
	
	print '\n'
	print 'The genetic distance matrix '+gdpathname+' has been created.'
	
