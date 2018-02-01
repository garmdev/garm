# ----------------------------------------------------------------------------
# v3.0 Mod - 2011 Nov 2nd BKH - Mod for Australian Bird Runs
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
directory = 'C:/Users/brian/Documents/SchoolWork/PhDWork/PaperDrafts/BirdsProject/Data/Tu1_BKH/Tu1_BKH/'
if len(sys.argv) >= 2:
	directory = sys.argv[1]
	#if os.name == 'nt':
	#	directory.replace('\\','/')
	#	directory = directory + '/'
		

run_types = ['rand60Overlap','randNoOverlap','invExp60Overlap','invExpNoOverlap','invsq60Overlap','invsqNoOverlap']
resist_surfs = ['We'] #'We','Tu','Sh','Mu','Ax','Cr','Un','Ha','Gl']
resist_types =['tree_100_1_2','tree_100_1_5','tree_100_1_10','tree_100_1_100','Edmatrix_1_1']
type_ext_one = '.asc'
type_ext_two = '_'
do_general_grid = True
test_type = 'smouse'
loci = 12
noalleperlocus = 20
alleles = int(noalleperlocus)*np.ones(loci,int)
#alleles = array([6,10])
nogrids = 1655
header_cols = 7
#dir_out = 'C:/Users/brian/Documents/SchoolWork/PhDWork/PaperDrafts/BirdsProject/Data/St/input/invExp60Overlap/batchrun0mcrun0_inputvariables_50000invExp60OverlapSt_tree_100_1_2_St_XY/'
input_dir = directory + 'input'
output_dir = directory + 'output'

#smouse_test = 'C:/Users/brian/Documents/SchoolWork/PhDWork/PaperDrafts/BirdsProject/Data/St_Modified_BKH_11.6.2011/output/invExp60Overlap/tree_100_1_10/Test/grid_mod500.csv'




# Get length of files
#nocsvfiles = len(csvfileList)

def calc_smouse(file):
	# Open file for reading
	#inputfile = open(csvfileList[i],'r')
	infile = open(file,'r')
	# Read lines from the file
	lines = infile.readlines()
	genes = []
	for ln in lines:
		genes.append(ln.replace('\n','').replace('\r','').strip().split(',')[7:])
			
	#Close the file
	infile.close()


	genes_ary = np.asarray(genes)
	genes_ary = np.asarray(genes_ary[1:],dtype='int16')
	psgd_vals = np.zeros((genes_ary.shape[0],genes_ary.shape[0]))
	
	for rows in range(psgd_vals.shape[0]):
		for rows_second in range(rows+1,psgd_vals.shape[0]):
			#pdb.set_trace()
			x = sum(((genes_ary[rows] - genes_ary[rows_second])**2)/2)
			psgd_vals[rows][rows_second] = x
			psgd_vals[rows_second][rows] = x
	
	#pdb.set_trace()
	
	pass

	
def calc_dps(nogrids,alleles,genes):
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
			
	return gendmatrix

def cleanDirectory(all_dirs,keyword):
	for type in run_types:
		for r_type in resist_types:	
			
			
			sub_in = input_dir + '/' + type + '/' +r_type
			sub_out = output_dir + '/' + type + '/' +r_type
			in_dir_list = glob.glob(sub_in + '/' + 'batchrun*')
			out_dir_list = glob.glob(sub_out + '/' + 'batchrun*')
			for d in in_dir_list:
				
				if checkCompleteness(d,'grid500*'):
					continue
				else:
					print ("Previous folder (%s) exists and is incomplete, removing..."%(d))
					shutil.rmtree(d)
			for d in out_dir_list:
				
				if checkCompleteness(d,'grid_mod500*'):
					continue
				else:
					print ("Previous folder (%s) exists and is incomplete, removing..."%(d))
					shutil.rmtree(d)



def checkCompleteness(dir,keyword):
	file_list = glob.glob(dir + '/' + keyword)
	if len(file_list) < 1:
		print("The folder: (%s) is incomplete, skipping..."%(dir))
		return False
	return True
	
def createHierachy():

	all_dirs = glob.glob(directory+'*batchrun*')
	print('Performing initial input file cleanup...')
	cleanDirectory(all_dirs,'*grid*500')
	#cleanDirectory(all_dirs,'*general_grid500*')
	#all_dirs = []

	print('End clean up')
	#pdb.set_trace()
	try:
		os.mkdir(input_dir + '/')
	except (IOError,OSError) as eMsg:
		pass
		#print("Error (%s) creating folder(%s)"%(eMsg,input_dir + '/'))
	try:
		os.mkdir(output_dir + '/')
	except (IOError,OSError) as eMsg:
		pass
		#print("Error (%s) creating folder(%s)"%(eMsg,output_dir + '/'))
	for dir in all_dirs:
		# new_file_dir =directory + '/batch' + dir.replace('.csv','').replace('.cdmatrix','').replace('.asc','').split('batch')[1]
		#all_dirs.append('batch' + dir.replace('.csv','').replace('.cdmatrix','').replace('.asc','').split('batch')[1])
		# pdb.set_trace()
		
		new_file_dir =dir.replace('.csv','').replace('.cdmatrix','')
		#pdb.set_trace()
		try:
			shutil.move(dir,new_file_dir)
		except:
			print('Folder name already exists... skipping')
			continue

	#Do some file clean-up with the directories
	all_dirs = glob.glob(directory+'*batchrun*')
	for type in run_types:
		try:
			os.mkdir(input_dir + '/' +type)
		except (IOError,OSError) as eMsg:
			pass
			#print("Error (%s) creating folder(%s)"%(eMsg,input_dir + '/' +type ))
		try:
			os.mkdir(output_dir + '/' +type)
		except (IOError,OSError) as eMsg:
			pass
			# print("Error (%s) creating folder(%s)"%(eMsg, output_dir + '/' +type))
		
		for r_type in resist_types:
			sub_in = input_dir + '/' + type + '/' +r_type
			sub_out = output_dir + '/' + type + '/' + r_type
			try:
				os.mkdir(sub_in)
			except (IOError,OSError) as eMsg:
				pass
				# print("Error (%s) creating folder(%s)"%(eMsg,sub_in))
			try:
				os.mkdir(sub_out)
			except (IOError,OSError) as eMsg:
				pass
				# print("Error (%s) creating folder(%s)"%(eMsg,sub_out))
				

	
	for dirs in all_dirs:
		#print dirs
		break_out = False
		for type in run_types:
			if not break_out:
			
				for r_type in resist_types:
					
					
					if dirs.count(type) > 0:
						if dirs.count(r_type +type_ext_one)> 0 or dirs.count(r_type + type_ext_two) > 0:
							bat_f = ''
							run_f = ''
							if dirs.lower().count('batch') > 0:
								bat_f = 'batch' + dirs.split('batch')[-1].split('_')[0]
								#bat_f = bat_f.lstrip('0123456789')

							else:
								print('Can\'t find batch number')
							
							if dirs.split('_')[-1].lower().count('run') > 0:
								run_f = dirs.split('_')[-1]
								run_f = run_f.replace('Run','')
							sub_in = input_dir + '/' + type + '/' +r_type
							sub_out = output_dir + '/' + type + '/' + r_type

							folder_in = dirs
							folder_out = sub_out + '/' + bat_f + run_f
						

									
							
							if checkCompleteness(folder_in,'*grid500*'):
								file_present = False
								try:
									os.mkdir(folder_out)
								except (IOError,OSError) as eMsg:
									print("Folder: (%s) exists already, skipping..."%(folder_out))
									break

									
								clean_files(folder_in,folder_out)
								print("Modified folder: (%s) to (%s)"%(folder_in, folder_out))
								new_folder = sub_in + '/' + bat_f + run_f
								shutil.move(dirs,new_folder)
								
							break_out = True
							break
								
							'''
							if not checkCompleteness(folder_out,'*generalgrid500*'):
								break_out = True
								break
							'''
								
						
								
						
						
					
						
					
			else:
				
				break
		if not break_out: print 'Could not find a match for:', dirs



def clean_files(folder_in,folder_out):
	#pdb.set_trace()
	
	input_files = glob.glob(folder_in + '/' + 'grid*')
	
	if do_general_grid:
		general_grid = glob.glob((folder_in + '/' + 'generalgrid*'))
	#pdb.set_trace()
	for file_num,file in enumerate(input_files):
		
		# Open file for reading
		#inputfile = open(csvfileList[i],'r')
		infile = open(file,'r')
		# Read lines from the file
		lines = infile.readlines()
				
		#Close the file
		infile.close()
		
		# Create an empty matrix to append to
		x = []
		
		# Split up each line in file and append to empty matrix, x
		for ln in lines:
			
			thisline = ln.strip(',').strip().replace('\r','').replace('\n','').split(',')
			
			x.append(thisline)
				

		new_grid = []
		for i,vals in enumerate(x):
			#vals[-1] = vals[-1].strip('\n')
			#pdb.set_trace()
			if len(vals) != sum(alleles) + header_cols:
				
				end_line = x.pop(i+1)
				#end_line[-1] = end_line[-1].strip('\n')
				new_line = vals + end_line
				
				new_grid.append(new_line)
			else:
				new_grid.append(vals)
		
		# Strip directory/filename of grid and add 'Gdmatrix.csv'

		
	
		filename = 'grid' + file.split('grid')[-1]
		gap = filename.split('grid')[-1].strip('.csv')
		if len(gap) == 3:
			gdpathname = filename.replace('grid','grid_mod')
		if len(gap) == 2:
			gdpathname = filename.replace('grid','grid_mod0')
		if len(gap) == 1:
			gdpathname = filename.replace('grid','grid_mod00')
			
			
			
		# Create file to write matrix to
		outputfile = open(folder_out + '/' +gdpathname,'w')
		
		# Sequence each row in the matrix
		for seqrow in new_grid:
		
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

	if do_general_grid:
		for file_num,file in enumerate(general_grid):
			# Open file for reading
			#inputfile = open(csvfileList[i],'r')
			inputfile = open(file,'r')	
			# Read lines from the file
			lines = inputfile.readlines()
					
			#Close the file
			inputfile.close()
			
			# Create an empty matrix to append to
			x = []
			
			# Split up each line in file and append to empty matrix, x
			for ln in lines:
				thisline = ln.strip(',').strip().replace('\r','').replace('\n','').split(',')
				x.append(thisline)
					

			new_grid = []
			for i,vals in enumerate(x):
				vals[-1] = vals[-1].strip('\n')
				
				if len(vals) != (2*loci + header_cols):
				
					end_line = x.pop(i+1)
					end_line[-1] = end_line[-1].strip('\n')
					new_line = vals + end_line
					
					new_grid.append(new_line)
				else:
					new_grid.append(vals)
			
			# Strip directory/filename of grid and add 'Gdmatrix.csv'

			
		
			filename = file.split('\\')[-1]
			gap = filename.split('grid')[-1].strip('.csv')
			if len(gap) == 3:
				gdpathname = filename.replace('generalgrid','generalgrid_mod')
			if len(gap) == 2:
				gdpathname = filename.replace('generalgrid','generalgrid_mod0')
			if len(gap) == 1:
				gdpathname = filename.replace('generalgrid','generalgrid_mod00')
			
			
			
			# Create file to write matrix to
			outputfile = open(folder_out + '/' +gdpathname,'w')
			
			# Sequence each row in the matrix
			for seqrow in new_grid:
			
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
	#for dir in all_dirs:
	
	


# -----------------------------------	
# Step 3: Run genetic distance method
# -----------------------------------
	
# ------------ Genetic Distance Matrix: Proportion of shared alleles -----------------------
# List all files with .csv extensions (but only the grid ones)

# Debugging break...
#pdb.set_trace()


#pdb.set_trace()
# Get the first globbed file read in
#csvfileList = ['C:/Users/brian/Documents/SchoolWork/PhDWork/PaperDrafts/BirdsProject/Data/St/input/invExp60Overlap/batchrun0mcrun0_inputvariables_50000invExp60OverlapSt_tree_100_1_2_St_XY/']
#nocsvfiles = 1

#calc_smouse(smouse_test)


createHierachy()


		
		
		
	


