# -----------------------------------------------------------------------------------------------------
# GARM.py
# Author(s): Brian Hand		     
# originally created  : March 2011
# GOAL: A GA to build the highest fitness resistance map using a mantel correlation test 
#		between cost distance and genetic distance.  Uses the UNICOR program to calculate
#		the cost distance matrix.
# Keywords: Genetic Algorithm, Resistance Map Creation 
#
#
#
#
# v1.2.1 7:25 AM 8/16/2013 - BKH Starting implementation on RDA into GARM as another statistical test. 
#	This will include the addition and need for frequency data. RDA uses genotype/haplotype frequency
#	data and a PCoA of the landscape distance matrix. I will need to include a new option in the config 
#	file in case one wants to run a genetic distance matrix. However, this might confuse things a bit
#	because the frequency data will be independent of a genetic distance matrix.
# v1.2	3:22 PM 7/23/2013 - BKH Working on implementing Circuitscape (v3.5.8) into GARM as an alternative
#	method to least-cost path modeling.  This will include a new option in the config file to choose
#	use of "resistance" modeling vs. "least_cost" modeling.
#	3:27 PM 7/25/2013 - BKH Have tested the circuitscape implmenetation up to outputs correctly 
#	returning from circuitscape and being tested correctly in GARM (Mantel, partial Mantel, and RV).
#
# v1.1  12:41 PM 3/23/2013 - BKH Staring to add the framework for a new CoIA test using the R package
#       ade4, and with dependencies on the rpy2 package which is an interface from python to R.  The CoIA
#       will be an additional fitness metric in the GA.  However, this implementation will only work with
#       linux to begin, at least for my set-up becuase of problems with installing rpy2.  This should be
#       doable on windows, but since the user loses parallel capabilities on windows in UNICOR, this
#		should not be a major issue.  I also plan to move GARM towards using fully the R interface to 
#		run tests for all fitness tests, including the Mantels.
#
# v1.0 12:14 AM 5/15/2013 - BKH Release version of GARM. Started the process of code clean-up and
#	documentation. 
#
# v0.9
#	10:32 AM 5/9/2013 - BKH I'm working to fix the hotstart functionality.  To do this I had to make 
#	sure the previous_pop_dict was being updated correctly.  Fitness is now case as a float
#	on hotstart. I've also turned off the track_visited switch as it has no functionality currently, and
#       should be removed from the start-up file.
# 
# v0.8 5:55 PM 4/14/2013 - BKH I worked on getting hotstart to function, did not take it through full 
#	full tests, however.
#	
#	9:24 PM 3/25/2013 - BKH Fixing some bugs with GARM to do with the last code update.  Some bugs to fix,
#	partial Mantels are not being run even with a euclidean distance matrix supplied, and some models
#	are showing none for all fitness metrics (suspect this is to do with a mutation bug), and also need to
#	update the code to change the number of permutations.  I'm also working towards moving UNICOR and GARM
#	towards unique run identifiers.  I've done this in other non-svn versions, but will now make this a 
#	live update.  #	
#       Additional Notes: Have also depreciated the tolerance and threshold values in the input file 
#		completely. Change mutation to now create new PopMember object to calculate the cdmatrix and 
#		fitness tests.  Also, mutation now happens after new children are created and their cdmatrix and 
#		fitness calculated so the previous list of population members is updated correctly before 
#		mutating weights.  This hopefully cleaned up the problems with mutation completely.  Noticed 
#		there was also a bug with mutation not updating the previous_pop_dict as was needed to keep 
#		duplicates from happening.
#
#
# v0.7
#       1:20 PM 3/9/2013 - BKH Working on changing top model output into the same R format as the current
#       member list and full member list.  Also looking into bug where top member was not being put back
#       into the full member list.  I don't know under what conditions it was being removed and not
#       being put back into the list.
#	UPDATED - There was no function to put mutated member back into all members list, this has been
#	changed and top member fitness now outputs in R format.
#	11:35 AM 3/11/2013 - BKH The copy.deepcopy caused a recursive error in copying a member.  Moved
#	to a copy.copy instead, the only catch is the wgt_dict gets overriden by mutation so, the wgts
#	are reset for the member_copy before being added to the all_member list
#
# v0.6
#	1:14 PM 1/5/2012 - BKH I've turned off the TOLERANCE and THRESHOLD abilities for the time being,
#	and the function haveSimilarWeights.  These have been replaced by the hasVisited function for now as 
#	this is a simplier, not as powerful, way to check and remove duplicates.  
#	
#	I also fixed a bug during child creation where the wgt_tup was not being updated, which then needs
#	to be checked against the previously visited individuals to avoid duplicates.  
#
# v0.5
#	5:40 PM 7/2/2012 - BKH Working on the mutation bug which is causing in rare cases inviduals which 
#	have already been considered to be recreated because of a mutation change.  This causes duplicates
#	in present generation lists.
#	
#	Also beginning work on implementing the new regression framework to provide an alternative method to
#	partial mantel tests.
#
# v0.4
#   7:26 PM 6/13/2012 - BKH Working on full member hotstart and the duplicates bug in the 
#	present generation list. 
#   Full member hotstart implemented and quickly tested.  
#  
#	9:04 PM 5/16/2012 - BKH Working on full member hotstart import function.  This is needed
#	to fill the checked member population so they are not visited again when hotstarting.  
# v0.3
#   10:20 AM 4/10/2012 - BKH Continuing on with R output format and reading in hotstart
#	files in new format.  Using an input R formatted file, GARM correctly calculated fitness
#	values for each of the input members.  
# 
#	12:05 4/4/2012 - BKH There is still a problem with the visited part of the code
#	Values are still showing up more than once in the all_members_list.
#	*Solution: Visited values are now added to the list as soon as they are created.  
#   
# 9:09 PM 4/3/2012 - BKH
#	Working on debugging hostart, and output format.  A change for hotstart input
#	to match R output format and keep I/O in one general form for now. 
# v0.2 12:53 PM 4/2/2012 BKH -
# 	1) Hotstart output - using the R output format, now available, not well tested
# 	2) R output file for easy data analysis - should be same format as Hostart I/O
# 	One difficulty is the output should have one with all generational members
# 	and one with just previous generation members. now available, not well tested
# v0.1 2:59 PM 3/30/2012 BKH - 
# 	1) previously visited solution tracking - 
# 	*Now imposing rounding on float weight values so they are now 
# 	rounding to the nearest tenth.  This is needed for visited 
# 	indexing to keep it simple
# 7:28 AM 3/28/2012 BKH - 
# 	Work for this version includes:
# 	1) previously visited solution tracking
# 	2) garm hotstart output
# 	3) an output for R which shows all members created during run time
# 	Notes: depreciated listCheck() as it is not being used anywhere in the 
# 	code.				
# 3.3.2012 - Beginning work on a visited solution dictionary to track answers which have 
# 	already been considered. This will be a new input parameter which the user can
# 	turn on and off.
	
# 2011.7.17 Upcoming Features
# - average crossover - child gets average value of parents
# - I/O for the weight_dictionary to be added
# - Rule Sets for weights, including greater and less than assignments
#   maximum weight combinations, Repair Operator
# - Int and float switch for weights
# ------------------------------------------------------------------------------------------------------
appRele = ""
appName = "GARM"
appVers = "1.2.1"
authorNames = "Brian Hand"


#File relative paths for importing functions
UTILITIES_PATH =  "../utilities/"
UNICOR_PATH = "../unicor/"
CIRCUITSCAPE_PATH = "../Circuitscape-3.5.8/"

# Import Modules with Except/Try statements
# For debugging
import pdb

# Platform and system functions
try:
	import os, sys, random, math, tempfile, shutil
	from operator import attrgetter                    
except ImportError as eMsg:
	print("ImportError: (%s) OS and SYS and random required!"%(eMsg))
	import sys
	sys.exit(-1)



#Import the package specific folders
CSV_folder = os.path.dirname(os.path.abspath(UTILITIES_PATH+"CSVParse"))
FILEIO_folder = os.path.dirname(os.path.abspath(UTILITIES_PATH+"FileIO"))
UNICOR_folder = os.path.dirname(os.path.abspath(UNICOR_PATH+"UNICOR"))
LOG_folder = os.path.dirname(os.path.abspath(UTILITIES_PATH+"Log"))
RANDOM_folder = os.path.dirname(os.path.abspath(UTILITIES_PATH+"RandomFun"))
CIRCUITSCAPE_folder = os.path.dirname(os.path.abspath(CIRCUITSCAPE_PATH+"Circuitscape-3.5.8"))


#Insert folders into the system path if they are not there already
if CSV_folder not in sys.path:
     sys.path.insert(0, CSV_folder)
if FILEIO_folder not in sys.path:
     sys.path.insert(0, FILEIO_folder)
if UNICOR_folder not in sys.path:
     sys.path.insert(0, UNICOR_folder)
if LOG_folder not in sys.path:
     sys.path.insert(0,LOG_folder)
if RANDOM_folder not in sys.path:
     sys.path.insert(0,RANDOM_folder)
if CIRCUITSCAPE_folder not in sys.path:
     sys.path.insert(0,CIRCUITSCAPE_folder)


# Numpy and SciPy Functions
try:
	import numpy as np
	import scipy.stats as spstats
except ImportError as eMsg:
	print("ImportError (%s) Numpy required."%(eMsg))
	sys.exit(-1)

# UNICOR functions
try:
	from UNICOR import main as RunUNICOR
except ImportError as eMsg:
	print("ImportError: (%s) UNICOR.py is required"%(eMsg))
	sys.exit(-1)

#Circuitscape functions
try:
	from cs_run import main as RunCircuitScape

except ImportError as eMsg:
	print("ImportError: (%s) cs_run.py is required"%(eMsg))
	sys.exit(-1)

# Utility functions
try: 
	import FileIO
	import CSVParse 
	import Log
	import RandomFun
	from RipMgr import mgrParseStanzaInputs    # Joe's Rip Manager for inputs
	from RipMgr import *
except ImportError as eMsg:
	print("ImportError: (%s) FileIO.py and CSVParse.py are required."%(eMsg))
	sys.exit(-1)

import copy, time, glob, re

#v1.1 3-23-2013 BKH
try:

	from rpy2 import robjects 
	from rpy2.robjects.packages import importr
	from rpy2.robjects.numpy2ri import numpy2ri
	robjects.conversion.py2ri = numpy2ri
	'''
	2018_01_01. See comments below, in def nparray2rmatrix.
	'''
	robjects.numpy2ri.activate()

except (IOError,OSError) as eMsg:
	print("Error (%s) rpy2 package is required to run Rv coefficient test")

	
def buildRMap(gen_member,maps_list):
	'''This function builds the resistance map from an input Pop member, and the 
	maps list.  Output is a resistance map in GIS ascii format ready for UNICOR.
	'''
	
	map_size = maps_list[0].map_data.shape
	res_map = np.zeros((map_size))
	########### DEVELOPER NOTE #################
	# This algorithm is assuming
	# points of NODATA are negative 
	# and should be set to zero in case
	# there is existing data from other maps for 
	# these points and therefore will result 
	# in a positive value when adding the layers together
	############################################
	
	for i, map in enumerate(maps_list):
		temp_map = np.zeros((map_size))
		NODATA_val =float(map.header_dict['NODATA_value'])
		for j, wgt_val in enumerate(gen_member.wgt_dict[map.name]):
			if map.breakpoints != None:
				#Use breakpoints to determine class cutoffs for 
				#each class associated with a map
				if j == 0:
					temp_map[map.map_data < map.breakpoints[j]] = wgt_val 
				elif j == int((map.num_of_classes) - 1):
					temp_map[map.map_data > map.breakpoints[-1]] = wgt_val 
				else:
					temp_map[map.map_data > map.breakpoints[j-1]] = wgt_val
					temp_map[map.map_data > map.breakpoints[j]] = 0.0
		
			if map.class_values != None:
				#Use straight class values to weight a mask of categorical 
				#data
				temp_map[map.map_data == map.class_values[j]] = wgt_val 
			
		temp_map[map.map_data == NODATA_val] = 0.
		res_map += temp_map
	#this adjustment assumes values of 0.0 are points of no data and should be ignored/set
	#to the default NODATA_value
	if 'NODATA_value' in map.header_dict:
		res_map[res_map == 0.0] = map.header_dict['NODATA_value']
	return res_map
	
'''
2018_01_19, Ted adds a pram one_dim_map, to be used if the model_type is "1D." See remarks above 
def load1dMap.
'''
def calcCDMatrix(model_type,gen_list, maps_list, con_model_user_file,grid_filename,log_handle,msg_verbose,log, one_dim_map ):
	'''Calculate the cdmatrix for a list of Pop members using UNICOR
	'''
	for i, member in enumerate(gen_list):
		if member.cd_matrix == None:
			FileIO.logMsg(log_handle,"Calculating cost-distance matrix for individual %s "%(i+1),msg_verbose)
			string = ''
                        for nums in member.wgt_tup:
                                string = string + '_' + str(nums)
                        grid_out_name = str(grid_filename.strip('.rsg')) + string+'.rsg'
			res_map = buildRMap(member, maps_list)
			
			FileIO.outputGrid(grid_out_name,res_map,header_dict = maps_list[0].header_dict,keys = keys)
			
			if model_type.lower() == 'least_cost': 
                        	header_dict, data = FileIO.loadFile(con_model_user_file, header_lines = 21)
                        	header_dict['Grid_Filename'] = grid_out_name
                        	FileIO.outputGrid(con_model_user_file,['Pregenerated UNICOR File'], \
						header_dict = header_dict)
				member.cd_matrix = RunUNICOR(con_model_user_file)
			elif model_type.lower() == 'resistance':
						
				log.profile('Start Circuitscape Run',verboseOverride =True)
				replace_file_line(con_model_user_file,'habitat_file =', 'habitat_file = ' + \
											grid_out_name + '\n') 
		
				#replace_file_line(con_model_user_file,'source_file =', 'source_file = ' + \
				#							grid_out_name + '\n') 
				replace_file_line(con_model_user_file,'output_file =', 'output_file = ' + \
								'./' +grid_out_name.strip('.rsg') + '\n')
				#a tuple is returned by Ciccuitscape, with the first item in the tuple
				#being a numpy array of resistance distances.
				raw_cd_matrix = RunCircuitScape(con_model_user_file)[0]
				
				#Circuitscape returns a numpy array with the members number in the 
				#first column and first row for identification.  This is not needed in
				#GARM, so remove this row and column.
				member.cd_matrix = raw_cd_matrix[1:,1:]
					
				log.profile('End Circuitscape Run',verboseOverride =True)
			
			elif model_type.lower() == "1d":
				'''
				2018_01_19. Ted adds new cd-matrix producer, details to be implemented,
				for now simply a cd matrix read in from a file.
				'''
				if  one_dim_map is not None:
					member.cd_matrix=one_dim_map
				else:
					print 'The model type is specified as \"1D\" but there is no 1D map.'
					sys.exit(-1)
				#end if 1D model exists, else does not
			else:
				print 'Unknown connectivity model specified, ' + str( model_type ) + '. Check GARM config file.'
				sys.exit(-1)				
			#This calls to the os to remove the created resistance map because it will have a unique
			#identifier to enable UNICOR to put out a unique cdmatrix filename.  
			os.remove(grid_out_name)

'''
2018_01_19 Ted added def load1dMap, to load a 1D map from file. This is for
cases in which the user has a pre-exisiting path (hence 1-dimensional map)
on the landscape, and so has no use for generating a path via UNICOR or the 
Circuitscape optimizations).  For the first implementation I'm just reading 
in from file an NXN (N=numpops) numpy matrix, as would have been generated, 
for example, by UNICOR.
'''
def load1dMap( filename ):
	'''
	In the current, simple implementation this def
	expects a file with floats separated by whitespace (space-key or tab).
	It also expects uniform number of floats per line.
	'''
	mymap=np.loadtxt( filename )

	return mymap
#end load1dMap

			
def calcFitness(gen_list, genetic_dist_ary,gen_counts_ary,ed_dist_ary,test_type,R_lib_home,test_nperms):
	
	#This fit constant is for Mantel and partial Mantel tests to get rid of negative correlation
	#models if one is using p-values to rank models
	poor_fit_constant = 0.001
	
	for member in gen_list:
	

		#v1.1 BKH adding new test type Rv coefficient
		#Going to make this test happen no matter what as well
		#may need to change this later on
		

		#create an r object to perform differnt functions
		r = robjects.r
		
		#import the correct R librarys
		importr('permute',lib_loc=R_lib_home) 
		ade4 = importr('ade4',lib_loc=R_lib_home)
		vegan = importr('vegan',lib_loc=R_lib_home)	
		
		
		#wrap the as.dist function to make it easier to call
		asdist = r['as.dist']	
		
		#This is needed in ade4 for input
		nf = member.cd_matrix.shape[0]

		gen_matrix = nparray2rmatrix(genetic_dist_ary)
		land_matrix = nparray2rmatrix(member.cd_matrix)
		gen_dist_matrix = asdist(gen_matrix)
		land_dist_matrix = asdist(land_matrix)
		
		#print("Untransformed distance matrix is Euclidean")
		#print(ade4.is_euclid(land_dist_matrix))	
		gen_dist_trans = ade4.cailliez(gen_dist_matrix)
		land_dist_trans = r.sqrt(land_dist_matrix)

		#print("Transformed distance matrix is Euclidean")
		#print(ade4.is_euclid(land_dist_trans))	
		#perform the Mantel calculations using vegan
		mantel_test = vegan.mantel(gen_matrix,land_matrix,permutations=test_nperms)
		member.mantel_R = mantel_test.rx2('statistic')[0]
		member.mantel_pvalue = mantel_test.rx2('signif')[0]

		#perform the partial Mantel test using vegan
		##### -- temp rem out and replace
		#if ed_dist_ary != None:	
		if ed_dist_ary.any() != None:	
		##### end rem out and replace
			ed_matrix = nparray2rmatrix(ed_dist_ary)
			ed_dist_matrix = asdist(ed_matrix)
			partial_test = vegan.mantel_partial(gen_matrix,land_matrix,ed_dist_matrix,permutations=test_nperms)
			member.partial_R = partial_test.rx2('statistic')[0]
			member.partial_pvalue = partial_test.rx2('signif')[0]
		
		#perform the calculations for the CoIA test, using ade4 functions
		land_pco = ade4.dudi_pco(land_dist_trans,scannf=False,nf=nf)
		gen_pco = ade4.dudi_pco(gen_dist_trans,scannf=False,nf=nf)


		rv_test = ade4.RV_rtest(land_pco.rx2('tab'),gen_pco.rx2('tab'),test_nperms)		
		
		member.coia_rv = rv_test.rx2('obs')[0]
		member.coia_pvalue = rv_test.rx2('pvalue')[0]
		
		'''
		2017_12_27.  Ted-debugging. There is an  error thrown on conditional 
		"if getn_counts_ary !=None": "The truth value of an array with more 
		than one element is ambiguous.  Use a.any() or a.all()"  I'm trying
		the "is not None" syntax, which seems to allow the test on an array.
		(Numpy apparently objects to simple boolean test on array).
		'''

#		if gen_counts_ary != None:
		if gen_counts_ary is not None:
			gen_counts_matrix = nparray2rmatrix(gen_counts_ary)
			gen_counts_matrix_chord = vegan.decostand(gen_counts_matrix,"norm")
			#use only the n-1 elements for RDA or all varaince will be explained
			#and there will be no unconstrained axis.
			land_len = len(land_pco.rx2('tab'))-1
			land_tab = robjects.r['data.frame'](land_pco.rx2('tab')[0:land_len])
			
			
			land_rda = vegan.rda(gen_counts_matrix_chord,land_tab)	
			radj = vegan.RsquareAdj(land_rda)
			pvalue = vegan.anova_cca(land_rda)
		
			member.rda_adj_r = radj.rx2('adj.r.squared')[0]
			member.rda_pvalue = pvalue.rx2('Pr(>F)')[0]	
	
		#v1.0 rearranged to have all the logic for selecting
		#which fitness to choose in one place
		if test_type == 'mantel_R':
			if member.mantel_R > 0.0:
				member.fitness = member.mantel_R
			else: member.fitness = poor_fit_constant
		elif test_type == 'mantel_pvalue': 
			member.fitness = 1. - member.mantel_pvalue
		elif test_type == 'partial_pvalue':
			if member.partial_R > 0.:
				member.fitness = 1. - pvalue
			else: member.fitness = poor_fit_constant
		elif test_type == 'partial_R':
			if member.partial_R > 0.:
				member.fitness = member.partial_R
			else: member.fitness = poor_fit_constant
		elif test_type == 'CoIA_RV':
			member.fitness = member.coia_rv
		elif test_type == 'CoIA_pvalue':
			member.fitness = 1.-member.coia_pvalue
		elif test_type == 'RDA_adj_R':
			if member.rda_adj_r > 0.:
				member.fitness = member.rda_adj_r
			else: member.fitness = poor_fit_constant
		elif test_type == 'RDA_pvalue':
			member.fitness = member.rda_pvalue
		else: 
			print('Unknown fitness test specified.')		
			sys.exit(-1)
		
			
def nparray2rmatrix(x):


	'''This function creates a numpy array appropriate for an rpy2 object
	'''
	nr, nc = x.shape

	''' 
	2017_12_27. Ted-debugging -- the call to matrix() is
	causing an error, related to rpy2's inability to
	take numpy ints for nrow, ncol -- some byte-alignment
	issue.  This problem I believe is due to changes in rpy2
	since the version originally used. The problem seems
	to have been solved by adding, after importing
	robjects.conversion.py2ri = numpy2ri, 
	the line robjects.numpy2ri.activate() (see
	https://stackoverflow.com/questions/2447454/
	converting-python-objects-for-rpy2)
	'''
	xvec = robjects.FloatVector(x.transpose().reshape((x.size)))
	xr = robjects.r.matrix(xvec, nrow=nr, ncol=nc)

	return xr



def crossover(child,weight_dict,cross_type = 'traditional'):

	'''
	random alleles - grab one weight at a time from each gene and combine
					into the child's weight dictionary
	traditional - takes one half gene split from each parent and 
				combines them together to form the child map
	random - takes two parents and randomly chooses genes from parent
			one or two and gives them to the child object
	randsplit - use one random split point to combine genes from
				the two parents
	average - take the average of the parents weights
	'''	
	
	#Randomize which parent is first to randomize combination of parent
	#genes otherwise the highest fitness parent will always be first
	ran_int = np.random.random_integers(0,1)
	if ran_int == 0:
		par_one_genes = copy.deepcopy(child.parents[0].wgt_dict)
		par_two_genes = copy.deepcopy(child.parents[1].wgt_dict)
	else:
		par_two_genes = copy.deepcopy(child.parents[0].wgt_dict)
		par_one_genes = copy.deepcopy(child.parents[1].wgt_dict)

	if cross_type == 'random_alleles':
		for map_var, wgt_vals in weight_dict.iteritems():
			genes_list = []			
			for i in range(int(wgt_vals['number_of_classes'])):
				ran_int = np.random.random_integers(0,1)
				if ran_int == 0:
					genes_list.append(par_one_genes[map_var][i])
				else:
					genes_list.append(par_two_genes[map_var][i])
			child.wgt_dict[map_var] = genes_list
		
	elif cross_type == 'random':
		for map_var, wgt_vals in weight_dict.iteritems():
			ran_int = np.random.random_integers(0,1)
			if ran_int == 0:
				child.wgt_dict[map_var] = par_one_genes[map_var]
			else:
				child.wgt_dict[map_var] = par_two_genes[map_var]
	elif cross_type == 'randsplit':
		ran_split = np.random.random_integers(0,len(weight_dict))
		for i, map_var in enumerate(weight_dict):
			if i < ran_split:
				child.wgt_dict[map_var] = par_one_genes[map_var]
			else:
				child.wgt_dict[map_var] = par_two_genes[map_var]
	elif cross_type == 'average':
		for key,wgts in par_one_genes.iteritems():
			child_wgts = []
			for i,par_one_wgt in enumerate(wgts):
				ave_wgt = (par_one_wgt + par_two_genes[key][i])/2
				child_wgts.append(ave_wgt)

			child.wgt_dict[key] = child_wgts
	elif cross_type == 'traditional': 
		for i, map_var in enumerate(weight_dict):
			if i < (len(weight_dict)/2):
				child.wgt_dict[map_var] = par_one_genes[map_var]
			else:
				child.wgt_dict[map_var] = par_two_genes[map_var]

	
'''
2018_01_19. Ted added param one_dim_map (for call to calcCDMatrix).
'''
def mutation(pres_gen_list,m_rate,weight_dict,previous_pop_dict,all_members_list,\
		maps_list,con_model_user_file,grid_filename,log_handle,msg_verbose,\
		genetic_dist_ary, gen_counts_ary, ed_dist_ary, test_type,R_lib_home,\
		test_nperms,model_type,log, one_dim_map):
	'''mutation takes a list of Pop members and does random mutation on the weights of all
	members.  This function is self-contained and does all recalculations for cdmatrix and
	fitness inside the function.
	'''
	
	for member in pres_gen_list:

		n = np.random.geometric(m_rate)
		total = 0
		cur_mem_wgt_dict = copy.deepcopy(member.wgt_dict)
		
			
		for map_var, wgt_val in cur_mem_wgt_dict.iteritems():
			for i in range(len(wgt_val)):
				if n > (total+1):
					total+=1
					continue
				if n == (total+1):
					low,high = weight_dict[map_var]['min_class_weight'][i], \
						 weight_dict[map_var]['max_class_weight'][i]

					#create a random integer between the min and max values of the
					#weight
					#if weight_dict[map_var]['dtype'] == 'float':
					#	weight = low + np.random.rand(1) * (high-low)	
					#elif weight_dict[map_var]['dtype'] == 'int':   
					weight = np.random.random_integers(low,high) 		
					
					n += np.random.geometric(m_rate)
					total +=1
					wgt_val[i] = np.asarray([weight])
					
					
		#Now check to see if mutation creates a weight set which has already been visited.		
		new_member = PopMember()
		new_member.wgt_dict = cur_mem_wgt_dict
		new_member.list_weights()
			
		
		if not has_visited(previous_pop_dict,new_member):
			#added in v0.7 BKH to add a copy of the previous member to the all members list			
			calcCDMatrix(model_type,[new_member], maps_list,con_model_user_file,grid_filename,\
			log_handle,msg_verbose,log, one_dim_map )
			calcFitness([new_member], genetic_dist_ary,gen_counts_ary,ed_dist_ary, test_type,R_lib_home,test_nperms)
			all_members_list.append(member)
			pres_gen_list.remove(member)	
			pres_gen_list.append(new_member)	
			add_visited(previous_pop_dict,new_member)
		
		
def propFitSelection(gen_list, num_elites, gen_pop_size,previous_pop_dict,all_members_list):
	'''
	Do proportional fitness calculations using a given
	fitnesss array and a maps_ary, the number of elite members
	to save and the generation population size.
	Returns an updated fitness_ary and maps_ary with only 
	the members of the next generation
	'''
	next_gen_list = []
	#remove all 0 fitness individuals from the list before 
	#selection
		

	sorted_gen_list = sorted(gen_list, key=attrgetter('fitness'))
	for i in range(num_elites):
		next_gen_list.append(sorted_gen_list.pop(-1))
	while len(next_gen_list) < (gen_pop_size):
		prob_map = np.zeros((len(sorted_gen_list)))
		tot_fitness = 0.
		
		for member in sorted_gen_list:
			tot_fitness += member.fitness
		for i,member in enumerate(sorted_gen_list):
			prob_map[i] = member.fitness/tot_fitness + prob_map[i-1]
		#This rand_index returns the first value in a list of lists,
		#the common output for a where search, finding all values greater
		#than the prob_map value and being closest
		rand_index = np.where(prob_map >  np.random.random_sample())[0][0]
	
		#Pop the member in the random index from the present generation
		#to be included in the next generation
		next_gen_list.append(sorted_gen_list.pop(rand_index))
	
	while len(sorted_gen_list) > 0:
		old_member = sorted_gen_list.pop()
		all_members_list.append(old_member)
	return next_gen_list, previous_pop_dict, all_members_list
		
		

def sumClasses(weight_dict, keyword = 'number_of_classes'):
	'''
	Small helper function to sum the total number of objects 
	for all keys in a dictionary given a value keyword for the
	sub dictionary, assuming it returns a numerical representation
	'''
	sum = 0
	for map_var in weight_dict.iterkeys():
		sum += int(weight_dict[map_var][keyword])
	return sum 
	
def trnySelect(members, t_size):
	'''Do tourney style selection using a list of members and a tourney size
	Output is the child, a new PopMember object, that knows its parents
	'''
	if t_size < 2:
		t_size == 2
		print("Setting tourney size to default size of 2.")
		
	trny_members = random.sample(members,t_size)
	sorted_trny_members = sorted(trny_members, key=attrgetter('fitness'))
	child = PopMember()
	child.parents = (sorted_trny_members[0],sorted_trny_members[1])
	return child

def createMemberList(num_members,member_list, weight_dict,previous_pop_dict):
	'''This function creates a list of randomly weighted PopMember objects
	'''
	#Create the first generation of resistance maps
	for i in range(num_members):
		create_new = True
		while create_new:
			new_member = PopMember(True, weight_dict)
			create_new = False
			for member in member_list:	
				if has_visited(previous_pop_dict,new_member):
					create_new = True
					break
			if not create_new:
				previous_pop_dict = add_visited(previous_pop_dict,new_member)
				member_list.append(new_member)	
	#return member_list
				

def createHotPop(hot_mem_list,hot_headers,weight_dict,maps_list,pres_gen_list):
	'''This function creates a hot start population from an input hotstart file.
	'''
	#function to create a member population from hot start input
	hot_headers = hot_headers[0].split(',')
	for member in hot_mem_list:
		new_member = PopMember(False, weight_dict)
		
		member = member[0].split(',')
		for i,header in enumerate(hot_headers):
			
			if member[i] == 'None':
				continue
			elif header == 'fitness':
				new_member.fitness = float(member[i])
			elif header == 'mantel_R':
				new_member.mantel_R = float(member[i])
			elif header == 'mantel_pvalue':
				new_member.mantel_pvalue = float(member[i])
			elif header == 'partial_R':
				new_member.partial_R = float(member[i])
			elif header == 'partial_pvalue':
				new_member.partial_pvalue = float(member[i])
			elif header == 'coia_rv':
				new_member.coia_rv = float(member[i])
			elif header == 'coia_pvalue':
				new_member.coia_pvalue = float(member[i])
			elif header == 'rda_adj_r':
				new_member.rda_adj_r = float(member[i])
			elif header == 'rda_pvalue':
				new_member.rda_pvalue = float(member[i])
		for map in maps_list:
			index = 1
			weight_list = []
			for i,header in enumerate(hot_headers):
				if header == map.name + '_' + str(index):
					index += 1
					weight_list.append(np.asarray([int(member[i])]))
							
			new_member.wgt_dict[map.name] = weight_list
		new_member.list_weights()		
		pres_gen_list.append(new_member)
		
def createWeightDict(feature_map_filenames):
	'''This function creates the weight dictionary assuming the .rip manager 
	stanza format as input, and calling files from the input GARM .rip file.
	'''
	weight_dict = {}
	for filename in feature_map_filenames:
		frip = mgrParseStanzaInputs(filename)
		
		feature_name = frip.kwdGetValue('name')
		
		weight_dict[feature_name] = {}
		
		
		stanza_label = 'DEFAULT_STANZA'
		# unpack ALL the keywords owned by this stanza
		kLis = frip.mgrStanzaListNames(stanza_label)
		
		for kwd in kLis:
			value = frip.kwdGetValue(kwd)
			weight_dict[feature_name][kwd] = value
	
	return weight_dict
	
def has_visited(previous_pop_dict,new_member):
	'''Small helper	function to check if a PopMember object has already
	been created previously and is in the previous_pop_dict
	'''
	if previous_pop_dict.has_key(new_member.wgt_tup):
		return True
	else:
		return False
def add_visited(previous_pop_dict,new_member):
	'''Small helper function to add a new member to the previous_pop_dict 
	to track visited solutions
	'''
	new_member.list_weights
	previous_pop_dict[new_member.wgt_tup] = new_member
	return previous_pop_dict	

def report(population,log_handle,msg_verbose):
	'''Used to report various info on a population of PopMembers to make
	end of run reporting quick and easy, or in debugging
	'''
	fitness_list = []
	count = 0
	for member in sorted(population, key=attrgetter('fitness')):
		count +=1
		FileIO.logMsg(log_handle,"Reporting on individual %s "%(count),msg_verbose)
		FileIO.logMsg(log_handle,"fitness %s "%(member.fitness),msg_verbose)
		FileIO.logMsg(log_handle,"mantel_R %s "%(member.mantel_R),msg_verbose)
		FileIO.logMsg(log_handle,"mantel_pvalue %s "%(member.mantel_pvalue),msg_verbose)
		FileIO.logMsg(log_handle,"coia_rv %s"%(member.coia_rv),msg_verbose)
		FileIO.logMsg(log_handle,"coia_pvalue %s"%(member.coia_pvalue),msg_verbose)
		if member.partial_R != None:
			FileIO.logMsg(log_handle,"partial_R %s "%(member.partial_R),msg_verbose)
			FileIO.logMsg(log_handle,"partial_pvalue %s "%(member.partial_pvalue),msg_verbose)
		FileIO.logMsg(log_handle,"Member weight dictionary \n%s "%(member.wgt_dict),msg_verbose)
		
		fitness_list.append(member.fitness)
	FileIO.logMsg(log_handle,"Average member fitness value %s "%(np.average(fitness_list)),msg_verbose)
	return np.average(fitness_list)

def write_fit_file(output_fit_info, test_type,top_member,num_immigrants,next_gen_list,log_handle,msg_verbose,out_fit_filename):
	'''Outputs a fitness file for a generation of models, giving top model fitness, average fitness,
	the number of immigrants
	'''
	output_headers= ['top_member_R','top_member_pvalue','top_member_fit','num_immigrants','average_fitness']
	
	if test_type == 'partial_R' or test_type == 'partial_pvalue':
		R = top_member.partial_R
		pvalue = top_member.partial_pvalue
	
	elif test_type == 'CoIA_RV' or test_type == 'CoIA_pvalue':
		R = top_member.coia_rv
		pvalue = top_member.coia_pvalue

	elif test_type == 'RDA_adj_R' or test_type == 'RDA_pvalue':
		R = top_member.rda_adj_r
		pvalue = top_member.rda_pvalue

	else:
		R = top_member.mantel_R
		pvalue = top_member.mantel_pvalue 
	
	dict_keys = top_member.wgt_dict.keys()
	out_keys = copy.deepcopy(keywords)

	info_list = [R,pvalue,top_member.fitness,num_immigrants,report(next_gen_list,log_handle,msg_verbose)]
	
	for keys in dict_keys:
		for num in range(len(top_member.wgt_dict[keys])):
			output_headers.append(keys+'_'+str(num+1))

	for key in dict_keys:
		for wgt_val in top_member.wgt_dict[key]:
			info_list.append(wgt_val[0])
	
	output_fit_info.append(info_list)
	
	FileIO.outputGrid(out_fit_filename,output_fit_info,headers=output_headers,delimiter=',')
	
def out_R_format(out_filename,keywords,member_list):
	'''Generic helper function to output a file in an R format, used for the full
	and current member list outputs
	'''
	out_list = []
	dict_keys = member_list[0].wgt_dict.keys()
	out_keys = copy.deepcopy(keywords)
	
	for member in member_list:
		line_list = []
		for key in keywords:
			
			line_list.append(getattr(member, str(key)))
		#add enviromental layers from the wgt_dictionary to the list of output values
		for key in dict_keys:
			for wgt_val in member.wgt_dict[key]:
				line_list.append(wgt_val[0])
		
		out_list.append(line_list)
	
	for keys in dict_keys:
		for num in range(len(member_list[0].wgt_dict[keys])):
			out_keys.append(keys+'_'+str(num+1))

	
	FileIO.outputGrid(out_filename,out_list,headers = out_keys,delimiter = ',')

def replace_file_line(file_path, old_line, new_line):
    #Create temp file
    fh, abs_path = tempfile.mkstemp()
    new_file = open(abs_path,'w')
    old_file = open(file_path)
    for line in old_file:
	if old_line in line:
        	new_file.write(new_line)
	else:
		new_file.write(line)
    #close temp file
    new_file.close()
    os.close(fh)
    old_file.close()
    #Remove original file
    os.remove(file_path)
    #Move new file
    shutil.move(abs_path, file_path)


class PopMember:
	'''PopMember is the model class which holds the information for each different model
	solution.
	'''
	def __init__(self, create_weights = False, weight_dict = None,KEY_DECIMALS = 1):
		self.wgt_dict = {}
		self.parents = []
		self.cd_matrix = None
		self.fitness = None
		self.pvalue = None
		self.mantel_R = None
		self.mantel_pvalue = None
		self.partial_R = None
		self.partial_pvalue = None
		self.id = None	
		self.wgt_tup = None
		self.coia_rv = None
		self.coia_pvalue = None
		self.rda_adj_r = None
		self.rda_pvalue = None
		if create_weights:
			self.randomizeWeights(weight_dict,KEY_DECIMALS)
			self.list_weights()
		
		
		
	def randomizeWeights(self,weight_dict,KEY_DECIMALS):
		'''This function is to start a PopMember object with random starting weights and
		is the default action for most runs.
		'''
		#create a list to append our random weights
		for map_variable, traits_dict in weight_dict.iteritems():
			weight_list = []
			for i in range(int(traits_dict['number_of_classes'])):
				#create a random integer between the min and max values of the
				#weight
				low,high = traits_dict['min_class_weight'][i] \
				,traits_dict['max_class_weight'][i]
				#Imposing rounding on float numbers to make indexing easier for visited solutions
				#3.30.2012 BKH
				
				#Taking out all floating weight values for now v1.2 BKH
				#if traits_dict['dtype'] == 'float':	
				#	weight = np.around(low + np.random.rand(1) * (high-low), decimals = KEY_DECIMALS)
				#if traits_dict['dtype'] == 'int':
				
					#extra list and array casting is necessary to match the float
					#convention above for the time being 2-19-2012 BKH
				weight = np.asarray([np.random.random_integers(low,high)])
				weight_list.append(weight)
			self.wgt_dict[map_variable] = weight_list
			
	def list_weights(self):
		'''This function is to create a tuple for each model for its set of weights, it is
		use to create a dictionary entry in the previous pop dict to track all visited solutions.
		'''
		weights = []
		#print self.wgt_dict
		for key,wgts in self.wgt_dict.iteritems():
			for wgt in wgts:
				weights.append(wgt[0]) 		
		self.wgt_tup = tuple(weights)
		
		
class Maps:
	'''Maps class is to create objects for map variables that are self-contained and 
	self-referential.  
	'''
	def __init__(self, class_name, traits_dict,header_lines):
		self.name = class_name
		self.filename = traits_dict['filename']
		#depreciated v1.2 BKH
		#self.class_type = traits_dict['class_type']
		self.num_of_classes = int(traits_dict['number_of_classes'])
		self.min_weight = traits_dict['min_class_weight']
		self.max_weight = traits_dict['max_class_weight']
		self.header_dict, map_data = FileIO.loadFile(self.filename,header_lines)
		self.breakpoints = None
		self.class_values = None
		#shut off for now, assuming int inputs only v1.2 BKH
		#self.dtype = traits_dict['dtype']	
		
		#Turn the map data list into a numpy array 
		self.map_data = np.asarray(map_data, dtype='float')
		if 'breakpoints' in traits_dict:
			self.breakpoints = traits_dict['breakpoints']
		if 'class_values' in traits_dict:
			self.class_values = traits_dict['class_values']
		if 'class_names' in traits_dict:
			self.class_names = traits_dict['class_names']
		
		#if 'dtype' in traits_dict:
		#	self.dtype = traits_dict['dtype']
		#else: self.dtype = 'float'
# ----------------------------------------------------------
# Global symbols, if any :))
#-----------------------------------------------------------
#These are the keys for an input GIS ascii file format.
keys = ['ncols','nrows','xllcorner','yllcorner',\
	'cellsize','NODATA_value']
#These are keywords for outputing current and full member lists in R format
keywords= ['fitness', 'mantel_R', 'mantel_pvalue','partial_R','partial_pvalue','coia_rv','coia_pvalue','rda_adj_r','rda_pvalue']

#------------------------------------------------------------
# Begin main file execution
#------------------------------------------------------------ 


def main(*arguements):
	
	# --------------------------------------------------------
	# Parse command line, to fetch session input file
	# --------------------------------------------------------
	if len(sys.argv) >= 2:
		ripFilePath = sys.argv[1]		
	
	# If user did not specify .rip file
	else:
		print "User must specify input file name (e.g., at command line type GARM.py user_input.rip)."
		sys.exit(-1)	
	
	# If .ip file does not exist
	if not os.path.exists(ripFilePath):
		print("Cannot find or open runtime inputs file(%s)"%(ripFilePath))
		sys.exit(-1)
	
	# --------------------------------------------------
	# Get user defined program input
	# --------------------------------------------------
			
	# create a RipMgr instance, via parser itself
	rip = mgrParseStanzaInputs(ripFilePath)   

	sessionLbl = rip.kwdGetValue('Session_label')
	sessionLbl = sessionLbl.strip(' ')
	logSessionPath = sessionLbl + ".log"
	try:
		log_handle =open(logSessionPath,'w')
	except (IOError,OSError) as eMsg:
		print("Error (%s) opening session logfile(%s)"%(eMsg,logSessionPath))
		sys.exit(-1)
	
	#Start up the log file and input starting information, verbose messages will print to screen
	msg_verbose = rip.kwdGetValue('verbose_messages')
	FileIO.logMsg(log_handle,"\n%s Release %s Version %s\n"%(appName,appRele,appVers),msg_verbose)
	FileIO.logMsg(log_handle,"Author(s): %s"%(authorNames)+'\n',msg_verbose)
	FileIO.logMsg(log_handle,"Session runtime inputs from: (%s)"%(ripFilePath)+'\n\n',msg_verbose)   
	FileIO.logMsg(log_handle,"Log output directed to     : (%s)"%(logSessionPath),msg_verbose)
	
	#This looks for a list of filenames of rip files which are feature map inputs.  A
	#dictionary is created to handle all the weight information for different map inputs.
	feature_map_filenames = rip.kwdGetValue('feature_map_filenames')
	weight_dict = createWeightDict(feature_map_filenames)
	
	# Read in the main input parameters from the input file
	con_model_user_file = rip.kwdGetValue('con_model_user_file')
	num_mutations = float(rip.kwdGetValue('num_mutations'))
	gen_pop_size = int(rip.kwdGetValue('gen_pop_size'))
	children_per_gen = int(rip.kwdGetValue('children_per_gen'))
	grid_filename = rip.kwdGetValue('grid_filename')
	out_top_filename = rip.kwdGetValue('top_map_filename')
	t_size = int(rip.kwdGetValue('t_size'))
	num_elites = int(rip.kwdGetValue('num_elites'))
	max_time = int(rip.kwdGetValue('max_time'))
	cross_type = rip.kwdGetValue('cross_type')
	max_gens = int(rip.kwdGetValue('max_generations'))	
	hot_start = rip.kwdGetValue('hot_start')
	test_type = rip.kwdGetValue('test_type')
	genetic_filename = rip.kwdGetValue('genetic_filename')
	edmatrix_filename = rip.kwdGetValue('ed_filename')
	test_nperms = int(rip.kwdGetValue('test_nperms'))	
	R_lib_home = rip.kwdGetValue('R_lib_home')
	model_type = rip.kwdGetValue('model_type')

	gen_counts_file = rip.kwdGetValue('gen_counts_filename')
	one_dim_map_filename = rip.kwdGetValue( 'one_dim_map_filename' )
		
	
	L = sumClasses(weight_dict)
	m_rate = num_mutations/float(children_per_gen*L)
	pres_gen_list = []
	maps_list = []
	header_lines = 6
	output_fit_info = []
	
	# This dict trackes previously visited solutions and fitness values are 
	# used as the key in the dict to easily compare to new members
	previous_pop_dict = {}
	all_members_list = []

	#file names for the output files for current and full member lists
	all_members_out = 'full_members.csv'
	gen_members_out = 'current_gen_members.csv'
	
	#Start the log file
	log = Log.Log()
 	log.open()
	out_fit_filename = 'fitness_vals_' + sessionLbl + '.out'
	
	#Load the genetic distance matrix and turn into numpy array
	header_dict,data = FileIO.loadFile(genetic_filename,delimiter=',')
	genetic_dist_ary = np.asarray(data,dtype='float')
	
	''' 
	2018_01_05.  Ted revision, noting that the code here seems to want
	to allow for absence of a gen_counts_file, in that it will then
	retain the data read in from the "genetic_filname", instead of replacing
	it with data from the gen_counts_file.  However, because the rip manager
	instance returns object None, not string 'None,' when there is
	no entry for a gen_counts_file, the absence then throws an error when
	the FileIO mod tries to open a file with None for a file name. I'm
	adding to the != 'None' formulation a second test,  "is not None"
	'''

	#if gen_counts_file != 'None':
	if gen_counts_file != 'None' and gen_counts_file is not None:
		try:	

			header_dict,data = FileIO.loadFile(gen_counts_file,header_lines=1,delimiter=',')
		
		except (IOError,OSError) as eMsg:
			print("Error (%s) opening session logfile(%s)"%(eMsg,ed_filename))
			print("Check that the genotype counts are available or if \
				not, please specify ed_filename None in GARM .rip file")
	
			sys.exit(-1)
		
		gen_counts_ary = np.asarray(data,dtype='int')
		
	else:
		gen_counts_ary = None	

	#v0.8 The logic of searching for a Euclidean distance matrix was changed to check
	#the input file to see if ed_filename is set to None for the run.  	

	'''
	2018_01_05.  See comment above the "if gene_counts_file != 'None'" statement --
	I suspect the same problem here -- should test for object None, not string 'None.'
	'''
	#if edmatrix_filename != 'None':
	if edmatrix_filename != 'None' and edmatrix_filename is not None:
		try:
			header_dict,ed_data = FileIO.loadFile(edmatrix_filename ,delimiter=',')
		except (IOError,OSError) as eMsg:
			print("Error (%s) opening session logfile(%s)"%(eMsg,ed_filename))
			print("Check that the euclidean distance matrix is available or if \
				not, please specify ed_filename None in GARM .rip file")
	
			sys.exit(-1)
	
		ed_dist_ary = np.asarray(ed_data,dtype='float')
	else: 
		ed_dist_ary = None
	
	'''
	Ted added 2018_01_19.
	'''
	if one_dim_map_filename != 'None' and one_dim_map_filename is not None:
		try:
			one_dim_map=load1dMap( one_dim_map_filename )
		except Exception as eMsg:
			print("Error (%s) opening session logfile(%s)"%(eMsg,one_dim_map_filename))
			print("Check that the 1-dimensional map is available or if \
					not, remove the one_dim_map_filename entry and use a model_type  \
					other than \"1D\"" )
	
			sys.exit(-1)

		#end try ... except	
	else:
		one_dim_map=None
	#end if 1d map file name exists, else No map

	#Write starting conditions to log file and to screen if msg_verbose is true.
	FileIO.logMsg(log_handle,"m_rate %s"%(m_rate *(children_per_gen*L)),msg_verbose)
	FileIO.logMsg(log_handle,"t_size %s"%(t_size),msg_verbose)
	FileIO.logMsg(log_handle,"num_elites %s"%(num_elites),msg_verbose)
	FileIO.logMsg(log_handle,"cross_type %s"%(cross_type),msg_verbose)
	FileIO.logMsg(log_handle,"max_generations %s"%(max_gens),msg_verbose)
	FileIO.logMsg(log_handle,"max_time %s"%(max_time),msg_verbose)
	FileIO.logMsg(log_handle,"Weight Dictionary %s "%(weight_dict),msg_verbose)
	#Create and load our data maps
	for variable, traits_dict in weight_dict.iteritems():
		maps_list.append(Maps(variable, traits_dict,header_lines))
	
	
	
	
	#Does the hotstart setup.
	if hot_start:
		hot_file = rip.kwdGetValue('hstart_filename')
		hot_file_full = rip.kwdGetValue('hstart_all_filename')
		hot_pres_header_dict,hot_pres_data = FileIO.loadFile(hot_file)
		hot_pres_headers = hot_pres_data[0]
		hot_pres_mem_list = hot_pres_data[1:]
		hot_full_header_dict,hot_full_data = FileIO.loadFile(hot_file_full)
		hot_full_headers = hot_full_data[0]
		hot_full_mem_list = hot_full_data[1:]	
		createHotPop(hot_pres_mem_list,hot_pres_headers,weight_dict,maps_list,pres_gen_list)
		createHotPop(hot_full_mem_list,hot_full_headers,weight_dict,maps_list,all_members_list)
		for member in all_members_list:
			add_visited(previous_pop_dict,member)	
	#Goes to normal random start with models and weights.  
	else:
		createMemberList(gen_pop_size,pres_gen_list, weight_dict,previous_pop_dict)	
	
	#calculate the cdMatricies and fitness of the first generation parents
	calcCDMatrix(model_type,pres_gen_list, maps_list,con_model_user_file,grid_filename,log_handle,\
			msg_verbose,log, one_dim_map )
	calcFitness(pres_gen_list, genetic_dist_ary,gen_counts_ary,ed_dist_ary,test_type,R_lib_home,test_nperms)



	#Initialize values for number of generations and the top fitness tracking variable
	gen_num = 0
	top_fitness = 0.0
		
	while gen_num < max_gens and top_fitness < 1.0:
		gen_num +=1
		FileIO.logMsg(log_handle,"Now working on generation %s "%(gen_num),msg_verbose)
			
		new_children = 0
		children_parent_list = pres_gen_list
		
		log.profile('Start of children generation loop',verboseOverride =True)
		stop_watch = log.processTimer
		stop_watch.stop()
		tot_time = stop_watch.getTime().seconds
		
		#Start child creation while loop
		#Children are created until reaching a certain number or the maximum time
		#specified by the user is reached.
		while new_children < children_per_gen and  int(tot_time) < int(max_time):
			stop_watch.stop()
			tot_time = stop_watch.getTime().seconds

			new_child = False
			#create a child object and give it a set of
			#parents for crossover
			child = trnySelect(pres_gen_list, t_size)
			
			crossover(child,weight_dict,cross_type)
			child.list_weights()
			
			for gen_member in children_parent_list:
				if has_visited(previous_pop_dict,child):
					new_child = True
					break

			if not new_child == True:	
				previous_pop_dict = add_visited(previous_pop_dict,child)
				new_children += 1
				children_parent_list.append(child)
		# end child creation while loop
		
		#Find the number of immigrants to produce, update the log
		num_immigrants = children_per_gen - new_children
		FileIO.logMsg(log_handle,"Immigrants %s "%(num_immigrants),msg_verbose)
		FileIO.logMsg(log_handle,"Total time for child loop %s"%(tot_time),msg_verbose)
		log.profile('End of children generation loop',verboseOverride =True)
		
		#create the new member list combining parents, children and immigrants into a single list
		createMemberList(num_immigrants, children_parent_list, weight_dict,previous_pop_dict)	
		
		
		# start fitness recalc loop
	 	for i,member in enumerate(pres_gen_list):
			FileIO.logMsg(log_handle,"Checking member %s for fitness recalc "%(i + 1),msg_verbose)
			
			if member.fitness == None:
				FileIO.logMsg(log_handle,"Recalculation fitness on individiual %s "%(i + 1),msg_verbose)
				
				#a list of objects of member objects is required, so throw the member into a 
				#list before sending 
				calcCDMatrix(model_type,[member], maps_list,con_model_user_file,grid_filename,\
				log_handle,msg_verbose,log, one_dim_map )
				calcFitness([member], genetic_dist_ary,gen_counts_ary,ed_dist_ary, test_type,R_lib_home,test_nperms)
		# end fitness recalc loop
		
		#mutation is a self-contained function that does all recalculations inside the function
		mutation(children_parent_list,m_rate,weight_dict,previous_pop_dict,all_members_list,\
			maps_list,con_model_user_file,grid_filename,log_handle,msg_verbose,\
			genetic_dist_ary, gen_counts_ary,ed_dist_ary, test_type,R_lib_home,test_nperms,model_type,log, one_dim_map )
		
		#do proportional fitness selection on the full member list to create the next generation of
		#models
		next_gen_list,previous_pop_dict,all_members_list = propFitSelection(pres_gen_list, num_elites, gen_pop_size,previous_pop_dict,all_members_list)
	
		#Find the top member in the model list to report and output its associated resistance map, 
		#and report some data about its characteristics
		top_member = sorted(next_gen_list, key=attrgetter('fitness'))[-1]
		FileIO.logMsg(log_handle,"Outputting resistance map for top member.",msg_verbose)
		res_map = buildRMap(top_member, maps_list)
		FileIO.outputGrid(out_top_filename,res_map,header_dict = maps_list[0].header_dict,keys = keys)
		FileIO.logMsg(log_handle,"Fitness of the top member: %s \n and weight dictionary of the top population member \n %s" \
		%(top_member.fitness, top_member.wgt_dict),msg_verbose)
		
		#Set the present generation list to the next generation list to start the loop over
		pres_gen_list = next_gen_list
		
		#write to fit file some basic information about this generation
		write_fit_file(output_fit_info,test_type,top_member,num_immigrants,next_gen_list,log_handle,msg_verbose,\
				out_fit_filename)

		#write to all member list
		out_R_format(all_members_out,keywords,all_members_list)
		#write to current member list
		out_R_format(gen_members_out,keywords,pres_gen_list)
		
		#set the top fitness variable to current top member and restart the loop
		top_fitness = top_member.fitness	
	
	# end main while loop()
	log.close()	
if __name__ == '__main__':
	main(*sys.argv[1:])
		

		
		
		
		
