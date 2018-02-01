# ----------------------------------------------------------------------------
# GdMatrix_Average.py
# Author: BKH
# Created: December 2011
# Description: Search through a directory tree to find all files with a certain
# keyword and then calculate the average on the matrices created from file.
# ----------------------------------------------------------------------------

import glob							# The power of glob
#from numpy import *					# General commands and functions
from numpy.random import *			# Random/statistic calculations
import time, datetime,pdb, os, shutil, re,sys	
import numpy as np

# Timing events: start
start_time = datetime.datetime.now()

# ------------------------------------------
# Step 1: Get user information
# ------------------------------------------

# Store directory path name
directory = ''
if len(sys.argv) >= 2:
	directory = sys.argv[1]

	
import os

file_path_list = []

start_dir = 'C:/Users/brian/Documents/SchoolWork/PhDWork/PaperDrafts/RMapGA/Data_12.22/CDPOP_MCMC_100runs/cdpop_bfid_270m_20111023_1_to_100'

#the keyword in the filename we want to key on to find specific GdMatrix files.
GD_file_key = 'Gdmatrix30.csv'

#size of population in Gdmatrix file
n = 50

#name of output file with new array
output_file_name = 'Gdmatrix_ave_100MC.csv'

for dirname, dirnames, filenames in os.walk(start_dir):
	for filename in filenames:
		filepath = os.path.join(dirname,filename)
		file_path_list.append(filepath)
		#print "Appending filepath to list: ", filepath


#file list name of matrix files to be read in
Gdmatrix_file_list = []
#main matrix we want to perform our function on (Ave)
data_matrix =  np.zeros((n,n))

for file_path in file_path_list:
	if file_path.count('Gdmatrix30.csv') > 0:
		
		Gdmatrix_file_list.append(file_path)
	


for file_path in Gdmatrix_file_list:

	# Open file for reading
	#inputfile = open(csvfileList[i],'r')
	infile = open(file_path,'r')
	# Read lines from the file
	lines = infile.readlines()
			
	#Close the file
	infile.close()
	
	# Create an empty matrix to append to
	x = []
	
	# Split up each line in file and append to empty matrix, x
	for ln in lines:
		
		thisline = ln.replace('\r','').replace('\n','').rstrip(',').strip().split(',')
		x.append(thisline)
	
	
	
	file_data = np.asarray(x,dtype = 'float').reshape(n,n)

	data_matrix += file_data
	
data_matrix = data_matrix/len(Gdmatrix_file_list)
	
		
# Create file to write matrix to
outputfile = open(output_file_name,'w')

# Sequence each row in the matrix
for seqrow in data_matrix:

	# Grab each element in each row and write element to outputfile
	for ele in range(len(seqrow)):
		outputfile.write(str(seqrow[ele]))
		if ele == len(seqrow)-1:
			
			continue
		else:
			
			# Add comma
			outputfile.write(',')
		
	# Return line
	outputfile.write('\n')

# Close file
outputfile.close()

			
			
			
		
