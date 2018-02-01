# ----------------------------------------------------------------------------
# CalculateEDmatrix_from file.py
# Author: Erin L Landguth
# Created: June 2009
# Description: This code calculates the ED from a set of points 
#	and outputs in a matrix.
# ----------------------------------------------------------------------------

# Import library
from scipy import *
import datetime
import time
import glob
import pdb
# --------------- Geographic Distance Matrix ---------------------
# Output the geographic distance matrix use xgrid and ygrid points

# Get some variables
directory = 'C:\\Users\\brian\\Documents\\SchoolWork\\PhDWork\\PaperDrafts\\BirdsProject\\Aus_birds_cdmatrices\\'
dir = ''
outfiletag= 'Edmatrix_'

fullfilenames = glob.glob(directory+'*XY.csv')
filenames = []

for f in fullfilenames:
	filenames.append(f.split('\\')[-1])

#pdb.set_trace()
for filename in filenames:
	
	print '\n'
	print 'Creating the Euclidean distance matrix for the xygrid...'	
	
	
	# Read in xypoints from file
	# Open file for reading
	inputfile = open(directory+filename,'r')

	# Read lines from the file
	lines = inputfile.readlines()

	#Close the file
	inputfile.close()

	# Create an empty matrix to append to
	x = []

	# Split up each line in file and append to empty matrix, x
	for i in lines:
		thisline = i.split(',')
		x.append(thisline)
			
	# Storage vectors
	xcoord = []
	ycoord = []
	nogrids = len(x)-1
			
	# Add x information to the created vectors
	for i in range(nogrids):
		xcoord.append(float(x[i+1][0]))
		ycoord.append(float(x[i+1][1]))	

	# Create a matrix of zeros to be filled
	geodmatrix = zeros((nogrids,nogrids),float)

	# Write to geodmatrix the euclidean distance values...
	for i in range(nogrids):
		for j in range(nogrids):
			geodmatrix[i][j] = float(sqrt((xcoord[i]-xcoord[j])**2+(ycoord[i]-ycoord[j])**2))

	# Transpose matrix to get vertical display (just being picky)
	geodmatrix = transpose(geodmatrix)

	# Create file to write matrix to
	outputfile = open(directory+outfiletag+'filename','w')

	# Sequence each row in the matrix
	for i in geodmatrix:

		# Grab each element in each row and write element to outputfile
		for j in range(len(i)):
			outputfile.write(str(i[j]))
			# Add comma
			outputfile.write(',')
		
		# Return line
		outputfile.write('\n')

	# Close file
	outputfile.close()

	# Free up memory...delete geodmatrix
	del(geodmatrix)