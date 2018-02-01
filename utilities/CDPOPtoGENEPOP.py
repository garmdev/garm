# -------------------------------------------------------------------------------------------------
# CDPOPtoGENEPOP.py
# Author: Erin L Landguth
# Created: Aug 21 2011
# Description: This program will take the grid{nthfile}.csv genetic information 
#	files from CDPOP and convert them to the GENEPOP.dat format.
# Program Input: 
#	1. folder location of the grid{nthfile}.csv files
# Program Output: 
#	1. genes of individuals in generalgrid{nthfile}.csv format 
# Program Steps:
#	1. Get User-defined information.
#	2. Preliminary work...
#	3. List .csv files in directory
#	4. Read in file, 
#	5. create output
#	6. create new format vector, 
#	7. and output format to new csv file

'''
Title line: "Grape populations in southern France"
ADH Locus 1
ADH #2
ADH three
ADH-4
ADH-5
mtDNA
Pop
Grange des Peres , 0201 003003 0102 0302 1011 01
Grange des Peres , 0202 003001 0102 0303 1111 01
Grange des Peres , 0102 004001 0202 0102 1010 01
Grange des Peres , 0103 002002 0101 0202 1011 01
Grange des Peres , 0203 002004 0101 0102 1010 01
'''

# --------------------------------------------------------------------------------------------------

# Import packages
import glob						# The power of glob
from numpy import *				# General commands and functions
import time, datetime			# A time clock 
import pdb						# iPython debugging

# Timing events: start
start_time = datetime.datetime.now()

# ------------------------------------	
# Step 1: Get user information
# ------------------------------------

# Directory location monte carlo folders - assume all folders to analyze
directory = 'C:/Landguth/CDFISH/GeneticVulnerability/barrier20110820/'

# Number of loci
noloci = 20

# Number of maximum possible alleles per loci
maxalleles = 20

# -------------------------------	
# Step 2: Prelimary work
# -------------------------------

# List folders in this dir
def listdirs(folder):
    return [d for d in (os.path.join(folder, d1) for d1 in os.listdir(folder)) if os.path.isdir(d)]
folderList = listdirs(directory)

# Create a genes vector, appending loci information with alleles to it
genes_genform = []
for iloci in xrange(noloci):
	locitemp = arange(1,maxalleles+1,1)
	genes_genform.append(list(locitemp))

# Create a storage vector for alleles
alleles = maxalleles*ones(noloci,int)
	
# Begin loop through folders if specificied
for ifolder in xrange(folderList):
	
	# --------------------------------------	
	# Step 3: List files in directory
	# --------------------------------------
	datfileList = glob.glob(folderList[ifolder]+'/'+'grid*.csv')
			
	# Get length of files
	nodatfiles = len(datfileList)
	
	# ---------------------------	
	# Step 4: Read in file 
	# ---------------------------
	# Loop through each grid
	for igrid in xrange(nodatfiles):

		# Grab filename in List
		filename = datfileList[igrid]
				
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
			
		# And grab the number of individuals in file here
		noindividuals = len(x)-1
		
		# And grab the rest of the information from the file
		sex_cdpop = []
		id_cdpop = []
		x_cdpop = []
		y_cdpop = []
		for ispot in xrange(noindividuals):
			sex_cdpop.append(x[ispot+1][2])
			id_cdpop.append(x[ispot+1][1])
			x_cdpop.append(float(x[ispot+1][4]))
			y_cdpop.append(float(x[ispot+1][5]))
					
		# Store genetic information: genes[individual][locus][allele]
		genes_cdpop = []
		for ispot in xrange(noindividuals):
			genes_cdpop.append([])
			for jspot in xrange(noloci):
				genes_cdpop[ispot].append(x[ispot+1][int(6+sum(alleles[0:jspot])):int(6+sum(alleles[0:jspot+1]))])
		
		# Delete x variable
		del(x)
		
		# --------------------------------------------------------	
		# Step 5: Create new output file for this grid
		# --------------------------------------------------------
		
		# Create file to write matrix to
		outputfilename = filename.split('\\')
		outputfile = open(outputfilename[0]+'/general'+outputfilename[1]+'.csv','w')
				
		# Write out the titles that match general grid format
		title = ['Population','ID','Sex','X','Y']
			
		# Write out the title
		for ititle in range(len(title)):
			outputfile.write(title[ititle])
			outputfile.write(',')
			
		# Write out the loci title 
		for i in range(noloci-1):
			outputfile.write('Locus'+str(i+1)+'a')
			outputfile.write(',')
			outputfile.write('Locus'+str(i+1)+'b')
			outputfile.write(',')
		# To get a return character on the end of the title
		outputfile.write('Locus'+str(noloci-1+1)+'a')
		outputfile.write(',')
		outputfile.write('Locus'+str(noloci-1+1)+'b')
		outputfile.write('\n')	

		# --------------------------------------------------------------------------	
		# Step 6: Create new format for each individual in grid file
		# --------------------------------------------------------------------------
			
		# Store general format gene output
		GenFormgenes = []
		
		# Loop through each individual
		for ithind in range(noindividuals):
			
			# Add gene individual spot 
			GenFormgenes.append([])
			
			# Loop through each locus
			for ithloci in xrange(noloci):
			
				# Add gene individual spot 
				GenFormgenes[ithind].append([])
				
				# Loop through each allele spot at that locu
				for ithallele in xrange(alleles[ithloci]):
					
					# Check if allele spot is 1
					if int(genes_cdpop[ithind][ithloci][ithallele]) == 1:
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
					
					# Check if allele spot is 2
					elif int(genes_cdpop[ithind][ithloci][ithallele]) == 2:
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
										
		# --------------------------------------------------------------------------	
		# Step 7: Write out genes and information to new outputfile
		# --------------------------------------------------------------------------
		
		# Loop through each individual
		for ithind in range(noindividuals):
			
			# Write out general information
			outputfile.write('NA')
			outputfile.write(',')
			outputfile.write(id_cdpop[ithind])
			outputfile.write(',')
			outputfile.write(str(int(sex_cdpop[ithind])+1))
			outputfile.write(',')
			outputfile.write(str(x_cdpop[ithind]))
			outputfile.write(',')
			outputfile.write(str(y_cdpop[ithind]))
			outputfile.write(',')			
			
			# Loop through each locus
			for ithloci in xrange(noloci-1):
			
				# Loop through each allele spot at that locus
				for ithallele in xrange(2):
				
					outputfile.write(str(GenFormgenes[ithind][ithloci][ithallele]))
					outputfile.write(',')
					
			# Return charater on end
			outputfile.write(str(GenFormgenes[ithind][noloci-1][0]))
			outputfile.write(',')
			outputfile.write(str(GenFormgenes[ithind][noloci-1][1]))
			outputfile.write('\n')	
			
		print '\n'
		print outputfilename[0],'/general',outputfilename[1],'has been created.'
						
		# Close file
		outputfile.close()		
			
	print '\n'
	print 'Total conversion time for foler',foldername[ifolder],': ',datetime.datetime.now() - start_time