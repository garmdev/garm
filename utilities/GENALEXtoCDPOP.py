# -------------------------------------------------------------------------------------------------
# GENALEXtoCDPOP.py
# Author: Erin L Landguth
# Created: Feb 12 2009
# Description: This program will take the .csv genetic information files from GENALEX
#	and convert them to CDPOP grid.csv format.
# Program Input: 
#	1. *.csv file or files in folder location
# Program Output: 
#	1. genes of individuals in grid.csv format - old name + .csv
# Program Steps:
#	1. Get User-defined information.
#	2. List .csv files in directory
#	3. Read in files, create new .csv file, and output format to this file
# ----------------------------------------------------------------------------

# Import packages
import glob						# The power of glob
from numpy import *				# General commands and functions
import time, datetime			# A time clock 
import pdb						# iPython debugging

# Timing events: start
start_time = datetime.datetime.now()

# ----------------------------	
# Step 1: Get user information
# ----------------------------
directory = 'C:\\GISWork\\Ruth_MTBear\\Genalex_examplefiles\\'
studyarea = 'BMU103'
noloci = 6
loci1 = [187,189,191,194,195,197,199,201,203,205,207,209]
loci2 = [154,156,158,160,162,164,166,168]
loci3 = [237,239,241,243,245,249,251,252,253,254,255,256,257,259,261,263,265,269]
loci4 = [135,137,139,141,145,149,155,157,159,161,163,165,169,171]
loci5 = [231,233,235,237,239,241,243,245,247,249]
loci6 = [185,187,189,191,195,197,199,201,203,205,207]
noalleles = [len(loci1),len(loci2),len(loci3),len(loci4),len(loci5),len(loci6)]
noindividuals = 208
genes = [loci1,loci2,loci3,loci4,loci5,loci6]

# -------------------------------	
# Step 2: List files in directory
# -------------------------------
datfileList = glob.glob(directory+'*.csv')
		
# Get length of files
nodatfiles = len(datfileList)

# ------------------------------------------------------------------------------------------------------	
# Step 3: Read in file, convert to CDPOP format, create new file, and write output
# ------------------------------------------------------------------------------------------------------
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
		thisline = l.split(',')
		x.append(thisline)
	
	# --------------------------------------------------------------------
	# Convert to CDPOPgenes = [individual][loci][alleles]
	# --------------------------------------------------------------------
	# Store genetic information: genes[individual][locus][allele]
	CDPOPgenes = []
	# Store EASYPOP gene output
	Genalexgenes = []
	
	# Loop through each individual
	for ithindividual in range(noindividuals):
		
		# Add gene individual spot
		CDPOPgenes.append([])
		
		# Get list from read in file
		Genalexgenes.append(x[ithindividual+3][2:(noloci*2+2)])
		
		# Loop through each locus
		for ithloci in range(noloci):
			
			# Loop through each allele
			for ithallele in range(noalleles[ithloci]):
				
				# Check first allele spot (case 1 at allele spot)
				if int(Genalexgenes[ithindividual][2*ithloci]) == int(genes[ithloci][ithallele]) and int(Genalexgenes[ithindividual][2*ithloci+1]) != int(genes[ithloci][ithallele]):
					CDPOPgenes[ithindividual].append(1)
				
				# Check second allele spot (case 1 at allele spot)
				elif int(Genalexgenes[ithindividual][2*ithloci]) != int(genes[ithloci][ithallele]) and int(Genalexgenes[ithindividual][2*ithloci+1]) == int(genes[ithloci][ithallele]):
					CDPOPgenes[ithindividual].append(1)			
				
				# Check if they are equal, then homogenous locus (case 2 at allele spot)
				elif int(Genalexgenes[ithindividual][2*ithloci]) == int(genes[ithloci][ithallele]) and int(Genalexgenes[ithindividual][2*ithloci+1]) == int(genes[ithloci][ithallele]):
					CDPOPgenes[ithindividual].append(2)	
			
				# Else case 0 at allele spot
				else: 
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
	for i in range(noloci-1):
		# Loop for allele length
		for j in range(noalleles[i]):
			outputfile.write('L'+str(i)+'A'+str(j))
			outputfile.write(',')
	# To get a return character on the end of the title
	for i in range(noalleles[noloci-1]-1):
		outputfile.write('L'+str(noloci-1)+'A'+str(i))
		outputfile.write(',')
	outputfile.write('L'+str(noloci-1)+'A'+str(noalleles[noloci-1]-1))
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
		for allelespot in range(sum(noalleles)-1):
			outputfile.write(str(CDPOPgenes[j][allelespot]))
			outputfile.write(',')
		# To get return character on end
		outputfile.write(str(CDPOPgenes[j][sum(noalleles)-1]))
		outputfile.write('\n')
				
	print '\n'
	print filename+'.csv has been created.'
					
	# Close file
	outputfile.close()
	
		
print '\n'
print 'Total conversion time: ', datetime.datetime.now() - start_time