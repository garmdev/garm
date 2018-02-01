# -------------------------------------------------------------------------------------------------------------
# MantelTester.py
# Author(s): Brian Hand			     
# originally created  : Oct 2011
# GOAL: A script to use the zt win test to calculate Mantel and Partial Mantel tests
# Keywords: Mantel Tests, Partial
# ------------------------------------------------------------------------------------------------------
appRele = ""
appName = "Mantel Tester"
appVers = "1.0"
authorNames = "Brian Hand"


#File absolute paths for importing functions
UTILITIES_PATH =  "../utilities/"
UNICOR_PATH = "../unicor/"
GARM_PATH = "../garm/"

# Import Modules with Except/Try statements
# For debugging
import pdb

# Platform and system functions
try:
	import os, sys, random, math
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
GARM_folder = os.path.dirname(os.path.abspath(GARM_PATH+"GARM"))

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
if GARM_folder not in sys.path:
     sys.path.insert(0,GARM_folder)

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

# Utility functions
try: 
	import FileIO
	import CSVParse 
	import Log
	import RandomFun
	from RipMgr import mgrParseStanzaInputs    # Joe's Rip Manager for inputs
	from RipMgr import *
	import GARM as GARM
except ImportError as eMsg:
	print("ImportError: (%s) FileIO.py and CSVParse.py are required."%(eMsg))
	sys.exit(-1)

import copy


test_type = 'partial_pvalue'
genetic_filename = 'C:/Users/brian/Documents/SchoolWork/PhDWork/PaperDrafts/BirdsProject/Data/St/500GensOnly/1318760945batchrun0mcrun0_inputvariables_50000invsqNoOverlapSt_tree_100_1_10.asc_St_XY.csv.cdmatrix.csv/Gdmatrix500.csv'
edmatrix_filename = 'C:/Users/brian/Documents/SchoolWork/PhDWork/PaperDrafts/BirdsProject/Data/EDMatricies/Edmatrix_St_XY.csv'
cdmatrix_filename = 'C:/Users/brian/Documents/SchoolWork/PhDWork/PaperDrafts/BirdsProject/Data/EDMatricies/St_tree_100_1_10.asc_St_XY.csv.cdmatrix'

header_dict,data = FileIO.loadFile(genetic_filename,delimiter=',')
genetic_dist_ary = np.asarray(data,dtype='float')

header_dict,data = FileIO.loadFile(cdmatrix_filename,delimiter=',')
cdmatrix_dist_ary = np.asarray(data,dtype='float')


if test_type == 'partial_R' or test_type == 'partial_pvalue':
	header_dict,ed_data = FileIO.loadFile(edmatrix_filename ,delimiter=',')
	ed_dist_ary = np.asarray(ed_data,dtype='float')
	
	
new_member = GARM.PopMember()
new_member.cd_matrix = cdmatrix_dist_ary

GARM.calcFitness([new_member], genetic_dist_ary,ed_dist_ary,test_type)



