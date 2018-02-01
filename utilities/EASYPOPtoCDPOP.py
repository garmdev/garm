# -------------------------------------------------------------------------------------------------
# EASYPOPtoCDPOP.py
# Author: Erin L Landguth
# Created: December 22 2008
# Description: This program will take the .dat genetic information files from EASYPOP
#	and convert them to CDPOP grid.csv format.
# Program Input: 
#	1. *.dat file or files in folder location
# Program Output: 
#	1. genes of individuals in grid.csv format - old name + .csv
# Program Steps:
#	1. Get User-defined information.
#	2. List .dat files in directory
#	3. Read in files, create new .csv file, and output format to this file
# ----------------------------------------------------------------------------

# Import packages
import glob						# The power of glob
from numpy import *				# General commands and functions
import time, datetime			# A time clock 

# Timing events: start
start_time = datetime.datetime.now()

# ----------------------------	
# Step 1: Get user information
# ----------------------------
directory = 'C:\\GISWork\\EasyPop_examplefiles\\'

# -------------------------------	
# Step 2: List files in directory
# -------------------------------
datfileList = glob.glob(directory+'*.dat')
		
# Get length of files
nodatfiles = len(datfileList)

# ---------------------------------------------------------------------------------	
# Step 3: Read in file, convert to CDPOP format, create new file, and write output
# ---------------------------------------------------------------------------------
# Loop through each .dat file in datfileList (order doesn't matter)
for i in range(nodatfiles):

	# Grab filename in List
	filename = datfileList[i]
			
	# Open file for reading
	inputfile = open(filename,'r')
				
	# Read lines from the file
	lines = inputfile.readlines()
				
	#Close the file
	inputfile.close()
				
	# Create an empty matrix to append to
	x = []
				
	# Split up each line in file by tab and append to empty matrix, x
	for l in lines:
		thisline = l.split('\t')
		x.append(thisline)
						
	# Store initial loci, allele, and number of individuals from .dat file
	noloci = int(x[0][1])
	noalleles =int(x[0][2])
	noindividuals = len(x)-noloci-1
	
	# --------------------------------------------------------------------
	# Convert EASYPOPgenes to CDPOPgenes = [individual][loci][alleles]
	# --------------------------------------------------------------------
	# Store genetic information: genes[individual][locus][allele]
	CDPOPgenes = []
	# Store EASYPOP gene output
	EASYPOPgenes = []
	
	# Loop through each individual
	for ithindividual in range(noindividuals):
		
		# Add gene individual spot
		CDPOPgenes.append([])
		
		# Get list from read in file
		EASYPOPgenes.append(x[ithindividual+noloci+1][1:(noloci+1)])
		
		# Loop through each locus
		for ithloci in range(noloci):
			
			# Loop through each allele
			for ithallele in range(noalleles):
				
				# Check first allele spot (case 1 at allele spot)
				if (int(EASYPOPgenes[ithindividual][ithloci][0:2])-1) == ithallele and (int(EASYPOPgenes[ithindividual][ithloci][2:4])-1) != ithallele:
					CDPOPgenes[ithindividual].append(1)
				
				# Check second allele spot (case 1 at allele spot)
				elif (int(EASYPOPgenes[ithindividual][ithloci][0:2])-1) != ithallele and (int(EASYPOPgenes[ithindividual][ithloci][2:4])-1) == ithallele:
					CDPOPgenes[ithindividual].append(1)			
				
				# Check if they are equal, then homogenous locus (case 2 at allele spot)
				elif (int(EASYPOPgenes[ithindividual][ithloci][0:2])-1) == ithallele and (int(EASYPOPgenes[ithindividual][ithloci][2:4])-1) == ithallele:
					CDPOPgenes[ithindividual].append(2)	
			
				# Else case 0 at allele spot
				elif (int(EASYPOPgenes[ithindividual][ithloci][0:2])-1) != ithallele and (int(EASYPOPgenes[ithindividual][ithloci][2:4])-1) != ithallele:
					CDPOPgenes[ithindividual].append(0)
					
	# --------------------------------------------------------------------
	# Create output table with titles of individuals + loci, allele titles
	# -------------------------------------------------------------------- 
	# Create file to write matrix to
	outputfile = open(filename+'.csv','w')
			
	# Write out the titles that match CDPOP grid output
	title = ['FID','id','sex','age','XCOORD','YCOORD']
		
	# Write out the title from xy points
	for j in range(len(title)):
		# Write out FID number
		outputfile.write(title[j])
		outputfile.write(',')
		
	# Write out the loci title info
	# Loop for loci length
	for j in range(noloci-1):
		# Loop for allele length
		for k in range(noalleles):
			outputfile.write('L'+str(j)+'A'+str(k))
			outputfile.write(',')
	# To get a return character on the end of the title
	for j in range(noalleles-1):
		outputfile.write('L'+str(noloci-1)+'A'+str(j))
		outputfile.write(',')
	outputfile.write('L'+str(noloci-1)+'A'+str(noalleles-1))
	outputfile.write('\n')
	
	# ------------------------------------------------------------
	# Write file information to the rest of the outputfile
	# ------------------------------------------------------------ 
	# Write out the xypoints + loci info		
	for j in range(noindividuals):
		outputfile.write(str(j))
		outputfile.write(',')
		outputfile.write('NA')
		outputfile.write(',')
		outputfile.write('NA')
		outputfile.write(',')
		outputfile.write('NA')
		outputfile.write(',')
		outputfile.write('NA')
		outputfile.write(',')
		outputfile.write('NA')
		outputfile.write(',')
		# Write genotype
		for allelespot in range(noloci*noalleles-1):
			outputfile.write(str(CDPOPgenes[j][allelespot]))
			outputfile.write(',')
		# To get return character on end
		outputfile.write(str(CDPOPgenes[j][noloci*noalleles-1]))
		outputfile.write('\n')
				
	print '\n'
	print filename+'.csv has been created.'
					
	# Close file
	outputfile.close()
	
		
print '\n'
print 'Total conversion time: ', datetime.datetime.now() - start_time