# RipMgr.py
modVers= 'v3.02-2011-03-30T12:11MST'
# (c) 2010,2011 Lupine Logic Inc.
#
# SYNOPSIS: RipMgr is a symbol table manager for Python applications, emphasizing the
#  management of runtime inputs. Variables in an application may optionally be registered
#  to a RipMgr instance 
#  The collective RipMgr facility is comprised of two cooperative Python classes:
#   -- RipObj() enscapsulates one keyword name-value pair, or, one field's properties.
#   -- RipMgr() encapsulates a collection of RipObj() as used in either of the above use-contexts or use cases
# STATUS of features for v3.02
# 1) Integer lists if comma delimited are now supported.
# 2) Float lists if comma delimited are now supported.
# 3) Stanza syntax is correctly parsed if the stanza label is all upper case and if the
#    stanza terminates with END_STANZA_LBL. Though parsed, stanza syntax it is NOT yet functional
#    in applications.
# STATUS of features FOR V3.01
# 1) add specific Exceptions in try,except blocks which lacked these (most did).
# 2) continue testing, identify any bugs, and correct these
# 3) add ValidateVar('myVarName') operator, which uses internal {min,max,valids-list} if
#    available to verify the correctness of a given variable.
# 4) refine parser to be more robust, via refactoring of parseLine() and classifyType()
# 5) enhance stanza parser to allow for multiple stanzas per stanza RIP file.close
# 6) WARNING: Windows filepaths cannot contain embedded blanks at this time. Backslashes in
#    Windows OS filepaths are normally converted to forward slashes internally here.
# History
# v3.02-2011-03.30T12:11MST(HP311,jmg)-- for Brian, removed messages about manifest file
# that aren't yet implemented and changed default doUnitTest to False and Verbose to False
#
# v3.01-2010-11-11T13:57:01MST --added uuid for RipMgr instances, and e
#  dded behavior-policy global variable
#  convertStrBooleans (a boolean), which if true, converts boolean values in string
#  form (e.g. TRUE,True,FALSE,False,ENABLE,DISABLE,YES,NO) to their native
#  datatype boolean Python equivalent {True|False}.
#  Lists in the RIP file may also now be entered as '1,2,3' or explicitly as [1,2,3]
#  Strip off all \n and \r end-of-line characters on RIP file lines, now by default
#  Also, added better exception handling throughout, though more work on exceptions still
#  needs to be done.
# v3.0--corrected a bug in __setitem__() where failed if dKey not yet present (now a new dKey causes a
#   new entry to be constructed in tDict{} and ripSetValue() succeeds.
#   --added __setitem__(),__getitem__() method overrides and allow for [] style assignment
#   and retrieves, and also fixed bugs in ripSetDescript(),kwdGetDescript()
#   deprecated ripAddItem() for ripSetItem, and deprecated kwdAddItem() for kwdSetItem()
# v.2.9.9--added TypeClassifier()
# v.2.9.8--unified storage of list objects in same location as all scalar objects.
#   and totally revample ripRole to embrace new semantics, replacing old distinction between
#   single|multiple with better classification scheme.
#   deprecated ripRole at v2.9.8 and higher,no longer needed or useful.
# v.2.9.7-2010.05.11-11:22:09 -- deprecated separate ripoUserObject and modified ripSetValue(), kwdSetValue() 
#       to allow for creation of a new name (dKey) if it does not exist, followed by newValue assignment.
#  so that we can create+assign value in one call
# v.2.9.7-2010.05.10-17:22:09 -- added Groups classifier facility, added logic-filter facility
# v.2.9.6-2010.05.10-11:31:02, corrected bugs in kwdSetUserObject() kwdAttachUserObject() etc.
# v.2.9.6 -- 2010.05.06-10:32:04, added ripAttachUserObject(),ripDetachUserObject(), and
#        next plan to add ripAddGroupTag(),ripGetGroupTag() to facilitate classification of
#        of parameter by stanza,group or other categorical descriptor.
#
# v2.9.5 -- removed all SQL functionality, removing tblXXX sqlXXX members 
#         and added kwdSetDescrip() kwdGetDescript()
#        Note: mgrStore()and mgrRestore() using pickle are not yet implemented, should be soon.
#        Note: the max no. of keywords supported per RipMgr object is now 65535. This may be raised as
#        necessary later.
# v.2.9.4
# v2.9.3-2009.01.19-18:30  -- added asserts() in fldSetOrdPos() and fldSetMode(), fldGetMode(),fixed SqlCreateTable etc.
# v2.9.2-2009.01.02-16:08  -- more acceptance testing, added mgrLoadItemList() mgrFetchItemList()
# v2.9.0-2009.01.02-02:23 -- recoding lots of changes lost when only edited source file was
#        accidentally overwritten by a botched up ptags.py I hacked.
#
# v.2.9.0-2009.01.01-16:32 -- performed major consolidation on FldObj, changed name FldObj to RipObj, and
# cleaned up namespace problems throughout both the new RipObj and RipMgr classes.
#
# v.2.8.0-2008.12.31-14:27,changed tOutputEntries to tShowKeywords() and tShowFieldMeta()
# v2.7.1-2008.12.26-18:16 changed appVers to modVers to avoid conflicts, deprecated use of dict.has_key('a')
#    in favor of new Python 2.6 and 3.0 usage (if 'a' in dict...), and added documentation.
# v2.7.0-2008.12.12-10:17 Added npKeywordTidyUp() function and added documentation.
# v2.6-2008.12.08-17:03 (changed name from TableDict() to RipMgr.py
# v2.5-2008.11.14.12:49-MonTEC
# v2.4-2008.11.14-10:32 home,mg -- debugged, now working for all test cases
# v2.3-2008.11.12-11:56 MonTEC,jmg -- added minimal file IO for testing keyword management, lots of new tests
# v2.2-2008.11.12-11:56 MonTEC,jmg -- test mode, left with failed KeywordItems functionality, went on to v2.3 
#                       implement keywordItems as a list instead of as a dict
# v2.1-2008.11.12-02:29 Home vaio,jmg
# v2.0-2008.11.11-15:48 DHC,jmg
# v1.2-2008.11.11-11:36 MonTEC,jmg
# v1.1-2008.11.11-10:43 jmg
# v1.0-2008.11.10-08:22pm joe glassy
#
# Built-in exception hierarchy is:
# BaseException
#  +-- SystemExit
#  +-- KeyboardInterrupt
#  +-- GeneratorExit
#  +-- Exception
#       +-- StopIteration
#       +-- StandardError
#       |    +-- BufferError
#       |    +-- ArithmeticError
#       |    |    +-- FloatingPointError
#       |    |    +-- OverflowError
#       |    |    +-- ZeroDivisionError
#       |    +-- AssertionError
#       |    +-- AttributeError
#       |    +-- EnvironmentError
#       |    |    +-- IOError
#       |    |    +-- OSError
#       |    |         +-- WindowsError (Windows)
#       |    |         +-- VMSError (VMS)
#       |    +-- EOFError
#       |    +-- ImportError
#       |    +-- LookupError
#       |    |    +-- IndexError
#       |    |    +-- KeyError
#       |    +-- MemoryError
#       |    +-- NameError
#       |    |    +-- UnboundLocalError
#       |    +-- ReferenceError
#       |    +-- RuntimeError
#       |    |    +-- NotImplementedError
#       |    +-- SyntaxError
#       |    |    +-- IndentationError
#       |    |         +-- TabError
#       |    +-- SystemError
#       |    +-- TypeError
#       |    +-- ValueError
#       |         +-- UnicodeError
#       |              +-- UnicodeDecodeError
#       |              +-- UnicodeEncodeError
#       |              +-- UnicodeTranslateError
#       +-- Warning
#           +-- DeprecationWarning
#            +-- PendingDeprecationWarning
#            +-- RuntimeWarning
#            +-- SyntaxWarning
#            +-- UserWarning
#            +-- FutureWarning
# 	   +-- ImportWarning
# 	   +-- UnicodeWarning
# 	   +-- BytesWarning
#
# Details:
#
#  --The RipObj() object encapsulates individual elements, and the RipMgr() class organizes collections of RipObj()
#  instances.
#  --Users only need access to RipMgr() methods, as these wrap all functionality in RipObj()
#  --The RipMgr() class includes two sets of member functions, one targeting keyword name-value pair
#  services, and one set of member functions targeting field data-dictionary services.
#  Overall, the RipMgr() capability set provides these services:
#  1) read, store keyword name value pair tokens out of a runtime input parameter text file (or)
#  2) read, store field data dictionary entries, storing a complete set of database field properties
#     for each field managed.
#  Capability Notes
#  1) when used in the keyword name-value pair role ''
#  Method Notes
#  ------------
#  Within RipObj() member function names dedicated to field-data dictionary services have a 'f' prefix
#    and function names providing keyword name-value pair services have a 'k' prefix.
#  Within RipMgr(), member function names providing field-data dictionary services have a 't' prefix
#    and function names providing keyword name-value pair services have a 'np' prefix (name-pair)
#    (We may clean this namespace up soon, replacing the 't' prefix with 'f' to be clearer and match the
#    RipObj() conventions, or otherwise get these in sync better.
#
# NEXT STEPS:
# Documentation goal: verify that each member function has a proper Python docstring entry, and that the
# format of these entries more of less match each other.
#
# Implement and test a setup.py packaging interface to allow this RipMgr() class to be a
# first class Python citizen e.g. a stand alone module, callable from anywhere.
# Perhaps use eGenix mx toolsibrary or other commercial library methods as a model for this.
#
# IMPLEMENTATION NOTES:
#
# The RipObj() class and its instances offer direct field-level access to various field properties.
# Typically however, the RipObj() member functions are accessed only through wrapper functions
# exposed by an instance of RipMgr().
#
# Note that internally a given RipObj() instance may be used in one of two contexts:
#   FldMode='field' -- to define and manage a field data dictionary entity
#   FldMode='keyword'--to define and manage a name-pair "keyword with value(s)" entity
# Note that special provisions allow a single instance of FldObj (or RipMgr) to store objects
# of either type at the same time, but this is not generally encouraged. It is recommended
# that users define a separate RipMgr() instance for use as a field data dictionary, and
# another RipMgr() instance for use in managing collections of name-pair entities.
#
# RipMgr() a set of member functions illustrate how interior RipObj() attributes
#  may correctly be assigned or retrieved via service functions present in instances of
# the higher level object, RipMgr()
#
import os,sys, time
import copy
import random
import cPickle
import hashlib
import math
# added 2010-11-09T13:18:01 jmg
import uuid
#
# only enable this when testing under ActiveState on Windows
# import win32traceutil
# -------------------------------------------------------------------------------
# GLOBAL SYMBOLS
# -------------------------------------------------------------------------------
# POLICY-FLAG: if True, converts string tokens "TRUE","ENABLE","YES" to actual
#   Python booleans. If this is False, leaves string tokens
#   "TRUE","ENABLE","YES" as literal strings
convertStrBooleans = True
# POLICY-FLAG: if convertStrScalars is True, converts string literals like "10.343" to float and
#   "2010" to integer, storing these in their native Python datatype int() or float()
#   upon parsing.
convertStrScalars = True
# POLICY-FLAG: to convert Windows path specs to universal Linux/Unix path with forward slashes
convertWinPathToUnix = True
#
# set True for development diagnostics, otherwise for production set false
verbose           = False
doMultiStanzaTest = False
doFullDiagTest    = False
# if testing, enable one or the other, as these two next are mutually exclusive
doShowPropertiesTest = False
doShowSummaryTest    = False
# if doUnitTest enabled, must have access to separate module RipMgrUnitTest.py
doUnitTest        = False
doLog             = True
defaultStanzaName = "DEFAULT_STANZA"
#
# eventually, this hard limit will likely be deprecated.
RipMaxKeyWords = 131070
#
typeForms = ['dict','list','tuple','scalar','class','stanza']
typeData  = ['integer','double','string','boolean']
#
# SIMPLE FUNCTIONS first...
# -------------------------------------------------------------------------------
def TestIfIter(obj):
	"""
	TestIfIter() --boolean function that returns True if obj is a list,tuple,
	dic or other iterable, and otherwise returns False. I suspect there is a
	better way to implement this, will review methods later.
	"""
	testResult= False
	if type(obj) == type(list()) or type(obj)== type(tuple()) or type(obj)==type(dict()):
		testResult = True
	else:
		testResult = False
	return testResult
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
class tokenObj(object):
	"""
	tokenObj - instance stores one complete name-value-pair object, as a result of
	parsing it from a delimited text token on a line that consists of
	one keyword, a type-classifier, a value, OR a list of values, or a string etc.
	In the future, additional syntax's may support inclusion of a name, datatype,
	min,max,default, range(low,high,incrementBy) on one line.
	"""
	def __init__(self,tokKeyword,tokType=None):
		self.tkKinds=['integer','double','string','intList','doubleList','intRange','doubleRange']
		self.tkType ='string'
		self.tkHasIncrement = False
		if tokType is not None:
			self.tkType = tokType
		self.tkForm    = 'scalar'
		self.tkKeyword = tokKeyword
		self.tkValue   = 0
		self.tkCount   = 0
		self.tkRngLow  = 0
		self.tkRngHigh = 0
		self.tkRange   = 0
		self.tkRngIncr = 0
	# end::__init__()
	def getProperties(self):
		"""
		getProperties() --retrieve a collection of name-value-pair object properties
		all bundled into a tuple. These are:
		 tkKeyword --the text keyword identifying the object
		 tkValue -- the value of the name-value-pair object
		 tkForm -- a form code e.g. 'scalar','list','dict','tuple',etc
		 tkCount --the number of elements if tkValue is a list,dict,set,or tuple
		 tkType  --the scalar datatype: ('integer','float','boolean','string')
		"""
		tupl= (self.tkKeyword,self.tkValue,self.tkType,self.tkCount,self.tkForm)
		return(tupl)
	#
	def setKeyword(self,newKeyword):
		"""
		setKeyword() --assign a keyword to this single token instance
		"""
		self.tkKeyword = newKeyword
		# end::setKeyword()
	def setType(self,newType):
		self.tkType = newType
	def setForm(self,newForm):
		self.tkForm = newForm
	def setValue(self,newValue):
		"""
		setValue() --assign a value to this name-value-pair token instance
		"""
		self.tkValue = newValue
	def setCount(self,newCount):
		self.tkCount = newCount
	def setRangeTuple(self,lo,hi,incr=None):
		"""
		setRangeTuple() --assign low and high value, and optionally, an
		 increment factor as well.
		"""
		self.tkRngLow = lo
		self.tkRngHigh= hi
		self.tkRange = self.tkRngHigh - self.tkRngLow
		if incr is not None:
			self.tkRngIncr= incr
			self.tkHasIncrement = True
		# end:: setRangeTuple()
	def setRangeLow(self,lo):
		"""
		setRangeLow() --assign a range's low value via 'lo'
		"""
		self.tkRngLow = lo
		# end::setRangeLow()
	def setRangeHigh(self,hi):
		"""
		setRangeHigh() --assign a range's high value via 'hi'
		"""
		self.tkRngHigh= hi
		self.tkRange = self.tkRngHigh - self.tkRngLow
		# end::setRangeHigh()
	def setRangeIncr(self,incr):
		self.tkRngIncr = incr
		self.tkHasIncrement = True
		# end::setRangeIncr() 
	# GET accessors
	def getKeyword(self):
		return self.tkKeyword
		# end::getKeyword()
	def getType(self):
		return self.tkType
		# end::getType()
	def getForm(self):
		return self.tkForm
		# end::getForm()
	def getValue(self):
		return self.tkValue
		# end::getValue()
	def getCount(self):
		return self.tkCount
		# end::getCount()
	def getRangeLow(self):
		return self.tkRngLow
		# end:getRangeLow()
	def getRangeHigh(self):
		return self.tkRngHigh
		# end::getRangeHigh()
	def getRangeIncr(self):
		return self.tkRngIncr
		# end::getRangeIncr()
	def getRangeTuple(self):
		"""
		getRangeTuple() --retrieves a name-value-pair's range factors as
		 low and high, and if incrementBy is defined, also returns this,
		 all bundled in a tuple (lo,hi,{incrmentBy})
		"""
		tupl = (self.tkRngLow,self.tkRngHigh,self.tkRngIncr)
		if self.tkHasIncrement:
			tupl=(self.tkRngLow,self.tkRngHigh,self.tkRngIncr)
		else:
			tupl=(self.tkRngLow,self.tkRngHigh)
		return tupl
		# end::getRangeTuple()
	# end:: tokenObj()

# -------------------------------------------------------------------------------
def classifyType(valu):
	"""
	classifyType() --classifies a name-value-pair's 'value' object as to 'form'
	(scalar,list,tuple,dict,or CLASS), and then, by its native Python datatype
	(integer,float,string,boolean). There may be other native types
	added later.
	Input  : one instance of a Python (object)
	Outputs: tForm,dType couplet, returned in a tuple (tForm,dType).
	Returns: a 2-element tuple as (formIdentifier,datatypeIdentifier)
	NOTE: this classifyType() function works somewhat independently, but for its fullest
	 use, is coupled with specific logic in parseLine() to work effectively.
	NOTE: this classifyType() treats OS directory paths just like strings - no
	  special handling or conversions of backslashs to forward slashs
	  is performed at this time; these conversions if done are performed in
	  parseLine().
	NOTE: tForm returned as: ('scalar','list','dict','tuple','range','literal')
	NOTE: if 'valu' is a object of tForm 'dict','tuple','list' or 'set' and this
	 is an empty collection, dType can be very hard to empirically determine since
	 we use the first element of a collection to classify the objects atomic datatype.
	 NOTE: a 'literal' is NOT treated like a collection, but is treated as "one item".
	We also assume that if an object is a collection, that all members are the same
	scalar type as the first element [0].append
	Although this is a dangerous assumption in Python, we make it nonetheless and is a
	documented limitation of RipMgr.
	"""
	isClass = False
	tForm = 'undetermined'
	dType = 'undetermined'
	evalInst = None
	#
	if type(valu) == type(str()) and valu.count(":") >= 1:
		# verify this is NOT a MSDOS path with a colon in it!!
		if valu.count("/")>= 1 or valu.count("\\") >=1:
			tForm = 'scalar'  # encountered a PATH most likely...
			dType = 'string'
		else:
			rLis = valu.split(":")
			tForm = 'range'
			evalInst = rLis[0]
			if type(evalInst)==type(int()):
				dType = 'integer'
			elif type(evalInst)== type(float()):
				dType = 'float'
			elif type(evalInst)==type(str()):
				dType = 'string'
				didConvertStrToScalar = False
				if evalInst.count(".") == 1:
					try:
						asFloat = float(evalInst)
						didConvertStrToScalar = True
						dType = 'float'
					except (ValueError,TypeError) as e:
						pass;
				else:
					try:
						asInt = int(evalInst)
						didConvertStrToScalar = True
						dType = 'integer'
					except (ValueError,TypeError) as e:
						pass;
				# end block if potentially a int or float
				
	else:
		if type(valu) == type(type):
			isClass = True
			tForm = 'class' 
			# retrieve NAME of class and use that 
			dType = valu.__name__
		elif type(valu) == type(dict()):
			tForm = 'dict'
			if len(valu) >= 1:
				evalInst = valu[valu.keys()[0]]
		elif type(valu) ==type(list()):
				tForm = 'list'
				if len(valu)>= 1:
					evalInst = valu[0]
		elif type(valu) == type(tuple()):
				tForm = 'tuple'
				if len(valu)>= 1:
					evalInst = valu[0]
		elif type(valu) == type(str()) and (valu.startswith('"') and valu.endswith('"')):
			# address use-case of explictly quoted string with embedded blanks
			# print("found-literal [%s]"%(valu))
			dType = 'string'
			tForm = 'literal'
	# done with collection types, literal-string, (or range)

	# first, attempt to evaluate entity-instance as a string that
	# MAY be representing a numeric quantity
	# NOTE: First attempt to coerce string into a float
	# this multi-stage logic is tricky...
	if tForm == 'undetermined':
		if type(valu) == type(str()) and valu.count('.')==0:
			try:
				tstAsInt = int(valu)
				tForm = 'scalar'
				dType = 'integer'
			except (TypeError,ValueError) as e:
				pass
		elif type(valu) == type(str()):
			try:
				tstAsFloat = float(valu)
				tForm = 'scalar'
				dType = 'float'
			except (TypeError,ValueError) as e:
				dType = 'string'
				tForm = 'scalar'
		elif type(valu) == type(int()):
			try:
				tstAsInt = int(valu)
				tForm = 'scalar'
				dType = 'integer'
			except (TypeError,ValueError) as e:
				pass
		elif type(valu) == type(float()) :
			try:
				tstAsFloat = float(valu)
				tForm = 'scalar'
				dType = 'float'
			except (TypeError,ValueError) as e:
				pass
		elif type(valu) == type(bool()):
			dType = 'boolean'
			tForm = 'scalar'

	if type(valu) == type(str()):
		qualifyVal = valu.upper()
		t1 = (qualifyVal == 'TRUE'   or qualifyVal =='FALSE')
		t2 = (qualifyVal == 'ENABLE' or qualifyVal =='DISABLE')
		t3 = (qualifyVal == 'YES'    or qualifyVal == 'NO')
		t4 = (qualifyVal == 'Y'      or qualifyVal == 'N' )
		t5 = (qualifyVal == 'T'      or qualifyVal == 'F' )
		# only convert special strings to booleans if the policy setting
		# indicted by 'convertStrBooleans' flag is True..
		if t1 or t2 or t3 or t4 or t5:
			tForm = 'scalar'
			if convertStrBooleans:
				dType = 'boolean'
			else:
				dType = 'string'

	# to get this far, tForm must have been {dict,list,tuple,set} and
	# valu must have been coded as a NATIVE int,float,boolean
	if evalInst is not None and tForm != 'scalar' and dType == 'undetermined':
		if type(evalInst)== type(int()):
			dType = 'integer'
		elif type(evalInst) == type(float()):
			dType = 'float'
		elif type(evalInst) == type(bool()):
			dType = 'boolean'
		elif type(evalInst) == type(str()):
			dType = 'string'
			try:
				iChk = int(evalInst)
				dType = 'integer'
			except (TypeError,ValueError) as e:
				try:
					fChk = float(evalInst)
					dType = 'float'
				except (TypeError,ValueError) as e:
					pass
	# if unexpected use-case slipped through all this, do a last-ditch
	# classify...
	if dType =='undetermined' and type(valu)==type(str()):
		dType = 'string'
	if tForm == 'undetermined':
		tForm = 'scalar'
	#
	return (tForm,dType)
	# classifyType()

# -------------------------------------------------------------------------------
class RipObj(object):
	"""
	RipObj() -- a class encapsulating a dictionary managed symbol table management object.
	Specific support is offered for managing name-value pairs originating in runtime input
	files, parsed as either vertically oriented (stanza form) sets or CSV style horizontally
	oriented sets.  RipMgr objects support the creation, naming, assignment, retrieval, and
	 removal of single-valued symbols of datatypes (integer,double,char,string) as well as
	 multiple-valued entities of the same datatype, bound to a name.  
	 Simple Example:  rp = RipMgr("NameOfCollection")
	  rp.kwdCreate("N_YEARS",'double')
	  rp.kwdSetValue("N_YEARS",12)
	  nYears = rp.kwdGetValue("N_YEARS")
	  rp.kwdRemove("N_YEARS")	
	"""
	def __init__(self,newName=None,newDatatype=None,newObjInstance=None):
		"""
		__init__() default RipObj constructor. Each RipObj represents one variable but this
		variable may be any Python type (scalar,tuple,list,dict,class-instance,etc). The
		value is stored in ripoObj.
		"""
		self.ripoModVers = modVers
		self.ripoName = 'no-name-assigned'
		self.ripoTime0 = time.time()
		self.ripoDescript = 'no-description-assigned'
		# deprecated self.ripoRole at v2.9.8 and higher
		self.ripoType = 'no-type-assigned'
		self.ripoDefaultValue = 'no-default-assigned'
		# boolean, used in conjunction with tGroupDict{} group-wise filtering of various properties
		# (basic visibility, whether to pickle or not, whether to include in a selection action, etc)
		self.ripoLogicFilter = True
		# length,precision,record ordinal position all intentionally default to
		# nonsense value of -1 to indicate not yet explicitly assigned
		self.ripoLength    = 1
		self.ripoPrecision = 0
		self.ripoOrdPos = -1
		self.ripoColOffset = -1
		self.ripoMinValue = -1
		self.ripoMaxValue = -1
		# new unified container for all Python objects,lists, scalar etc.
		self.ripoObj   = None # one value per name-pair entity, as:{scalar,tuple,list,dict,class etc}
		self.ripoCurrGroup = 'DefaultGroup'
		# RipObj valid 'form' ripoForm forms are:
		# ('scalar','tuple','list','dict','array','class')
		self.ripoForm = 'scalar'
		#
		if newDatatype != None:
			self.ripoType = newDatatype.lower()

		# v.2.3 and higher, - re-implemented Keyworditems as a list not a dict
		# v2.9.8 and  higher-unified storage of lists in same location as scalars so eliminated
		self.ripoItemCount = 0
		if newName != None :
			self.ripoName = newName
		# unified interface on 5.14.2010 under one roof, 
		if newObjInstance != None:
			self.ripoObj = newObjInstance
	# end::__init__() constructor for RipObj
	
	def ripSetName(self,newName):
		"""
		ripSetName --assign a name to keyword.
		"""
		self.ripoName = newName
		# end::ripSetName()
		
	def ripGetName(self):
		"""
		ripGetName() returns the current scalar object name as assigned previously.
		"""
		return self.ripoName

	def ripSetDatatype(self,newType):
		"""
		ripSetDatatype() assigns a datatype string code for object and
		stores string descriptor as all lowercase:: ('integer','double','boolean','string',...)
		"""
		self.ripoType = newType.lower()
		# end::ripSetDatatype()
		
	def ripGetDatatype(self):
		"""
		ripGetDatatype() returns a datatype string code 'Integer','String' etc as assigned earlier
		"""
		return self.ripoType
		# end::ripGetDatatype()
		
	def ripSetDescript(self,newDesc):
		"""
		ripSetDescript() -- assigns a longer-description string, which is a user defined,
		free form metadata string to provide more context to how this specific instance
		is used in an application.
		"""
		self.ripoDecript = newDesc
		
	def ripGetDescript(self):
		"""
		ripGetDescript() retrieve the longer-description string
		"""
		return self.ripoDescript
		
	def ripSetOrdPos(self,newOrdPos):
		"""
		ripSetOrdPos(newOrdinalPosition) -- assign a ordinal (sequence) index (1..N) to this
		RipObj instance, corresponding to the physical left-to-right (or top to bottom) order
		of this 'field' relative to all other fields.
		"""
		self.ripoOrdPos = newOrdPos
		# end::ripSetOrdPos() 
		
	def ripGetOrdPos(self):
		"""
		ripGetOrdPos() returns the integer (base0) ordinal position of this object
		among the set of peer objects, as assigned earlier.
		"""
		return self.ripoOrdPos

	def ripSetColOffset(self,newColOffset):
		"""
		ripSetColOfripSet - assigns a base1 column ofripSet associated with the original source data
		field used to populate this in-memory field. This is relevant only if the source data was
		a fixed format text file (a rare but possible use-case).
		"""
		self.ripoColOffset = newColOffset
	def ripGetColOffset(self):
		"""
		ripGetColOffset() --returns an optionally set column offset (where one is appropriate)
		"""
		return self.ripoColOffset

	def ripSetLength(self,newLen):
		"""
		ripSetLength(newLength) -- assign a variable or keyword length (in character units) 
		for this RipObj instance. Note that currently, this assigned length is not enforced 
		or checked against the reality of how long the variables VALUE entity may really be.
		This property is more directly useful within a GUI where we may limit user input to
		this many chars etc.
		"""
		self.ripoLength = newLen
		
	def ripGetLength(self):
		"""
		ripGetLength() --return the length of the fields value entity, if it matters.
		Note that currently, this returned length is not enforced or checked against the reality
		of how long the fields VALUE entity may really be.
		"""
		return self.ripoLength

	def ripSetPrecision(self,newPrecision):
		"""
		ripSetPrecision(newPrec) -- assign a variables decimal precision factor  (in character units)
		for this RipObj() instance.
		"""
		self.ripoPrecision = newPrecision
		# end::ripSetPrecision()
		
	def ripGetPrecision(self):
		"""
		ripGetPrecision() returns the precision of the object, as appropriate, relevant if a
		Real or Float. This is included more for compatibility with xBase interfaces.
		"""
		return self.ripoPrecision

	def ripFormatAsStr(self):
		"""
		ripFormatAsStr() -- uses preset variable properties (ripoLength,ripoPrecision) to format the current
		value of the RipObj as a formatted string for display, often in a GUI, menu, or other structured
		environment such as a report. Only very simplistic formatting is performed, and is relevant only to
		integer and double types. Native string variables are not intended to be used here.
		"""
		fmtStr    = ""
		resultStr = ""
		# mandatory: we require both the ripoLength and ripoPrecision to be defined,
		# though for integers, the ripoPrecsion may be 0 and is expected to be.
		if not self.ripoLength is None and not self.ripoPrecision is None: 
			# goal: %4d where 4 is the ripoLength
			if self.ripoType == 'integer' or self.ripoType == 'int':
				precStr = "%d"%(self.ripoLength) 
				fmtStr = "%" + precStr + "d"
			else:  # goal: %10.4f style relative-format spec
				lengStr= "%d." % (self.ripoLength) 
				precStr = "%d"%(self.ripoPrecision)
				fmtStr = "%" + lengStr + precStr + "f"
			#
			if type(self.ripoObj) == type(str() and not self.ripoObj is NoneType):
				resultStr = ((fmtStr) % (self.ripoObj))
			else:
				resultStr = ((fmtStr) % (self.ripoObj))
		else:
			print("ripFormatAsStr: RipObj length and precision properties NOT yet assigned and must be!")
		# finally, construct the formatted value statement
		return resultStr
		# end::ripFormatAsStr()
		
	def ripSetDefault(self,newDefault):
		"""
		ripSetDefault() assign a default value to RipObj
		"""
		self.ripoDefaultValue = newDefault
		# end::ripSetDefault()
		
	def ripGetDefault(self):
		"""
		ripGetDefault() return the default value currently assigned to RipObj
		"""
		return self.ripoDefaultValue

	def ripSetMinValue(self,newMinValue):
		"""
		ripSetMinValue() -- assign a minimum allowable (VALIDS) value for RipObj.
		This is NOT the same as calculating the minimum value that may be present in the set of
		RipObj list.
		"""
		self.ripoMinValue = newMinValue
		
	def ripGetMinValue(self):
		"""
		ripGetMinValue() -- return minimum allowable value defined for RipObj
		"""
		return self.ripoMinValue

	def ripSetMaxValue(self,newMaxValue):
		"""
		ripSetMaxValue -- assign a max allowable value for RipObj. Useful mostly for numeric types
		"""
		self.ripoMaxValue = newMaxValue
		
	def ripGetMaxValue(self):
		"""
		ripGetMaxValue() -- return max allowable value for RipObj. Useful mainly for numeric types
		"""
		return self.ripoMaxValue
		
	def ripSetValue(self,newObjInstance):
		"""
		ripSetValue(newValue) --assigns ANY single Python object to internal ripObj
		such that each ripObj has ONE name and ONE instance associated with that name.
		"""
		# check if a boolean submitted as a string, convert to a Python boolean if it is
		if newObjInstance is 'True' or newObjInstance is 'TRUE':
			newObjInstance = True
		if newObjInstance == 'False' or newObjInstance is 'FALSE':
			newObjInstance = False
		#
		self.ripoObj = newObjInstance
		# end::ripSetValue()

	def ripGetValue(self):
		"""
		ripGetValue() --return current value from RipObj
		"""
		return self.ripoObj
		# end::ripGetValue()

	def ripResetValue(self):
		"""
		ripResetValue() --resets RipObj's value to None
		"""
		self.ripoObj = None
		# end::ripResetValue()

	def ripSetItem(self,itemIndex,newValue):
		"""
		ripSetItem() -- assigns a new (list) value to a RipObj configured as a list. If this instance
		was NOT a list prior to this call, it is blindly converted to a list. If the list is empty,
		the newValue is appended to the start of the list. If itemIndex+1 exceeeds the length of the
		list, a warning is issued and no assignment is made.
		"""
		n =0
		if not type(self.ripoObj) == type(list()):
			self.ripoObj = []
		# NOTE: in the future, we may modify this logic so that if itemIndex is larger then
		# current length of list,
		# we auto-grow the list, filling with zeros, then insert the newValue at the indicated slot.
		if (itemIndex+1) < len(self.ripoObj):
			self.ripoObj[itemIndex] = newValue
		else:
			self.ripoObj.append(newValue)
			n = len(self.ripoObj)
			self.ripoCount = n
		return n
		# end::ripSetItem()

	def ripGetItem(self,itemIndex):
		"""
		ripGetItem() -- retrieve a multi-valued RipObj's item, as identified by a
						base0 itemIndex
		"""
		result = None
		n = len(self.ripoObj)
		if TestIfIter(self.ripoObj) and itemIndex < n:
			result = self.ripoObj[itemIndex]
		# return result which may be None, so caller must check value!
		return result
		# end::ripGetItem()

	def ripSortItems(self,whichOrder=None):
		"""
		ripSortItems() -- sort the current list of items in ripoObj[], by default in
		ascending order, or, if whichOrder is 'D'|'d', in descending order.
		NOTE: v3.1--this needs major refinement ,esp for list vs tuple identification
		and resulting action
		"""
		# if a list, 1st char of chk will be a '[' symbol.
		# if a tuple,1st char of chk will be a '(' symbol..
		chk = repr(self.ripoObj)
		# check if this object is a list and thus. is sortable
		if TestIfIter(self.ripoObj) and chk.startswith('['):
			if whichOrder is None or whichOrder.lower() == 'a':
				self.ripoObj.sort()
			else:
				cmd = whichOrder.lower()
				if cmd.startswith('d'):
					self.ripoObj.reverse()
				else:
					self.ripoObj.sort()
		return 0

	def ripUpdateItemCount(self):
		"""
		ripUpdateItemCount - updates and returns the count of a keywords item
		  collection as stored in its internal list, RipItems[], by evaluating len(RipItems).
		"""
		newCount = 0 
		if TestIfIter(self.ripoObj):
			newCount = len(self.ripoObj)
			self.ripoItemCount = newCount
		else:
			print("ripUpdateItemCount:: %s does NOT appear to be an interable!"%(self.ripGetName()))
		return newCount
		# end::ripUpdateItemCount()

	def ripGetItemCount(self):
		"""
		ripGetItemCount() -- returns count of multiple items maintained and owned by RipObj
		"""
		if TestIfIter(self.ripoObj):
			self.ripoItemCount = len(self.ripoObj)
		return self.ripoItemCount
		# end::ripGetItemCount()

	def ripGetItemList(self):
		"""
		ripGetItemList() --retrieves entire (list) of RipObj's items for this keyword
			returning the result as a native (list)
		"""
		# return native (list)of multiple keyword items for multi-valued keyWords
		if TestIfIter(self.ripoObj):
			rslt = self.ripoObj
		else:
			rslt = None
		return rslt 
		# end::ripGetItemList()

	def ripLoadItemList(self,myList,op=None):
		"""
		Inputs: myList 
				option op choices: {'replace'|'append'}
		ripLoadItemList() -- transfers a users list directly into the internal 
		ripObj.RipItems[] list. An op command may further specify if the copy is 
		destructive 'replace' or submitted as 'append' to add to whatever else (if anything)
		might have been previously held in the RipItems[] list.
		"""
		# Warning: force conversion to a list type, overwrites any previous value 
		self.ripoObj = [0]
		if TestIfIter(self.ripoObj):
			if op == 'replace':
				# TBD verify if just a copy.copy() is best here
				self.ripoObj = []
				self.ripoObj = copy.deepcopy(myList)
				self.ripoItemCount = len(self.ripoObj)
			elif op == 'append':
				for myItem in myList:
					self.ripoObj.append(myItem)
				self.ripoItemCount = len(self.ripoObj)
			else:
				# use soft copy semantics...
				self.ripoObj = []
				self.ripoObj = myList
		# end:: ripLoadItemList()

	def ripRemoveItem(self,itemIndex):
		"""
		ripRemoveItem() -- removes a (RipObj's) list element (value) at the base0 index
		                   specified as itemIndex. All following elements are moved UP in
		                   the list order after this (pop stack) operation.
		"""
		if TestIfIter(self.ripoObj):
			if itemIndex < len(self.ripoObj):
				self.ripoObj.pop(itemIndex)
		# end::ripRemoveItem()

	def ripResetItems(self):
		"""
		ripResetItems() -- unilaterally resets ALL list-items owned by a RipObj
		"""
		if TestIfIter(self.ripoObj):
			self.ripoObj = []
			self.ripoItemCount = 0
		# end::ripResetItems()

	def ripIsValueSet(self):
		"""
		ripIsScalarValueSet -- boolean, True if ripoScalarValue is not None
		"""
		result = False
		if self.ripoObj is None:
			result = False
		else:
			result = True
		return result
		# end::ripIsValueSet()
		
	def ripSetFilter(self,whichState=None):
		"""
		ripSetFilter() -- assigns a filter state (True|False) to this RipObjinstance.
		if whichState is not specified, the default is to set it True. 
		A RipObj's logic state is used to control selection-set-membership, various
		actions like which objects are pickled, or in some cases, masking a subset of
		RipObj.
		"""
		boolState = True
		if not whichState is None and (whichState is True or whichState is False):
			boolState = whichState
		#
		self.ripoLogicFilter = boolState
		# end::ripSetFilter()
		
	def ripGetFilter(self):
		"""
		ripGetFilter() -- retrieves the logic state of a RipObj instance.
		"""
		return self.ripoLogicFilter
	
	def ripAssignToGroup(self,thisGroupLabel):
		"""
		ripAssignToGroup() assigns a text group label to local instance of
		the current-group specifier for this RipObj instance. Use mgrAssignToGroup()
		in normal practice.
		"""
		self.ripoCurrGroup = thisGroupLabel
		# end::ripAssignToGroup()
		
	def ripGetGroup(self):
		"""
		ripGetGroup() retrieves the group classifer this RipObj is bound to. If user does not
		set this explicitly, it is set to DefaultGroup.
		"""
		return self.ripoCurrGroup
		# end::ripGetGroup()
	
	# end:: RipObj() object class

# -------------------------------------------------------------------------------
class RipMgr(object):
	"""
	RipMgr() -- a (Runtime-Input-Processor,RIP) class encapsulating a single generic
	service to process keyword name-value pairs, or, variables in a symbol table.
	Two classes are defined, where the RipMgr serves as the parent class (container)
	and RipObj serves as the child class.  RipMgr instances thus manage an array 
	of RipObj() class instances.
	--each RipMgr() object encapsulates or contains any number of RipObj child objects.
	Multiple RipMgr() objects may be easily defined within one application as needed.
	--RipObj instances may now be classified by "GROUP", where the user may create any number
	of group labels, and bind any RipObj instance to any group classifier. Classifying variables this
	way can be useful to control their visibililty, membership in selection sets (such as which
	variables ought to be pickled for persistent store), etc.  Variables NOT explicitly assigned
	to any other group are assigned to a 'DefaultGroup'.
	--each RipObj() instance encapsulates all metadata properties for a single variable
	or keyword name-value pair.
	--all RipObj() instances are referenced via their user defined text name.
	--RipObj() scalar values may be defined as 'int','integer', 'double', or 'string' datatypes
	--RipObj() may also be defined to hold a simple list of any Python object
	--RipObj() may also be defined to hold a single instance of any complex Python object such
	 as a dict{} or list[] nested to any level, etc.
	--In some cases, one RipObj() may represent several types of entity, but this is dangerous.
	--In this API, the term "keyword" connotes any variable with a name, so may be used
	interchangably with the term 'variable' since variables have a name and possess a value.
	--kwdAAAAA wrapper functions refer to 'keyword' (or variable) context functions
	--mgrAAAAA prefaced functions refer to RipMgr() system functions
	"""
	def __init__(self,newRipMgrName=None):
		self.RipMgrVersion   =  modVers  # inherit global instance
		self.MgrCreationSec     = time.time() # immutable creation seconds  stamp
		self.MgrCreationDatetime= time.ctime()# immutable creation datetime stamp
		self.MgrUUID_DNS = "" # specific DNS assigned by user, for UUID construction
		self.MgrUUID = 0  # immutable UUID assigned to this RipMgr instance as a whole
		#
		self.MgrSesStrtSec   = time.time()   # immutable session starting time in sec.
		self.nItems          = 0
		self.PersistPath = 'RipMgrPersist.pkl'  # name of cPickle or Shelve volume itself
		# mark with various timestamps
		self.MgrRevStr	     = time.ctime()
		self.MgrStartIO_Sec  = time.time()
		self.MgrEndingIO_Sec = time.time()
		self.MgrElapsed_Sec  = 0
		self.MgrLastRead = time.ctime()
		self.MgrLastWritten  = time.ctime()
		#
		self.RipMgrDefaultGroup = 'DefaultGroup'
		#
		if not newRipMgrName is None:
			self.RipMgrLabel = newRipMgrName
		else:
			self.RipMgrLabel = 'DefaultRipMgrInstance'
		# Dictionary for main object properties metadata
		self.tDict = {}
		# v2.9.6 and higher, added groups concept as a classifier for objects, allowing group-wise
		# symbol handling. As each entity instantiated, it may be bound to a new or pre-existing group
		self.tGroupDict = {}
		# v2.9.9 and higher, evaluating tStanzaDict, where to qualify, a matching END_{Blocktoken} must be
		# found for any {BlockToken}
		self.tStanzaDict ={}
		# added minimal session logging support 
		self.MgrUseLog = False
		# log message destinations: file=file only, both =screen and file, console=screen only
		self.MgrLogDest= "both"
		self.MgrLogf = 0
		self.MgrLogPath = ""
		# end:: __init__()
		
	def __setitem__(self,dKey,valu):
		"""
		__setitem__() override, for RipMgr dictionary item assignment.
		Inputs: dKey -- name of variable or keyword
				valu -- literal value or object to associate with this name,dKey
				Example: r3 = RipMgr('name')
				         r3['var1'] = 345 to assign,
				         getMyValue = r3['var1']
		"""
		# NOTE: may want to modify to use try, except KeyError trap instead..
		# for now, if dKey not yet present, add entry for minimal presence
		# then add the initial value
		if not dKey in self.tDict:
			self.kwdCreate(dKey)
		#
		self.tDict[dKey].ripSetValue(valu)
		# end::__setitem__()
		
	def __getitem__(self,dKey):
		"""
		__getitem__() override, for RipMgr dictionary item retrieval
			Example: r3 = RipMgr('name')
				r3['var1'] = 345 to assign,
				getMyValue = r3['var1']
		"""
		rslt = None
		try:
			rslt = self.tDict[dKey].ripGetValue()
		except KeyError as e:
			print("__getitem__ error %s"%(e))
		return rslt
	# end:: __getitem__() override

	# ---------------------------------------------------------------------
	def kwdCreate(self,newKeyword,dataType=None,groupLabel=None):
		"""
		kwdCreate() - create a name-pair entry uniquely identified by 'newKeyword'
		users are encouraged to use this rather than addFieldObj() to assure
		that all default constructor-time actions are performed.
		Inputs:
			newKeyword -- name of keyword
			dataType -- optional, datatype string for field, from:
			('Integer','Double','Real,'Text','String','VARCHAR')
			groupLabel -- if supplied, binds this new entity to a pre-existing 
			group classifier indicated by groupLabel.
		"""
		#
		# Assign new RipObj to a default (pooled) group in case no other option available
		dfltGroup = self.RipMgrDefaultGroup
		#
		dKey = newKeyword
		if newKeyword in self.tDict:
			kwObj = self.tDict[dKey]
		else: # call constructor, one of two ways...
			if dataType is not None:
				kwObj = RipObj(newKeyword,dataType)
			else:
				kwObj = RipObj(newKeyword)

		self.tDict[dKey] = kwObj
		ith = len(self.tDict)
		# as a convenience, auto-assign keywords ordinal position
		kwObj.ripSetOrdPos(ith)
		#
		# if user declared a GROUP to bind this RipObj to, implement the binding
		if groupLabel is None:
			groupLabel = dfltGroup
		#
		if groupLabel in self.tGroupDict:
			kwObj.ripAssignToGroup(groupLabel)
			self.tGroupDict[groupLabel].append(newKeyword)
		else:
			kwObj.ripAssignToGroup(groupLabel)
			# must first create the list, then can grow it
			self.tGroupDict[groupLabel] = []
			self.tGroupDict[groupLabel].append(newKeyword)
		#
		# Finally , as a convenience, return kwObj instance, even though it is kind of FAT
		return kwObj
		# end::kwdCreate()

	def kwdGetName(self,dKey):
		"""
		kwdGetName() -- retrieves an entities 'name' e.g. as a keyword identified by dKey
		In practice this is rarely used or needed since an objects
		dKey <---> keyword for this scheme!
		"""
		if dKey in self.tDict:
			keyWord = self.tDict[dKey].ripGetName()
		else:
			keyWord = 'no-keyword-assigned for key %s' % (dKey)
		return keyWord
		# end::kwdGetName()

	def kwdSetDatatype(self,dKey,newDataType):
		"""
		kwdSetDatatype() -- assigns a string datatype label to the (keyword) identified
		by dKey
		"""
		if dKey in self.tDict:
			self.tDict[dKey].ripSetDatatype(newDataType)
		else:
			print('kwdSetDatatype:: no key %s available!' % (dKey))
		# end::kwdSetDatatype()

	def kwdGetDatatype(self,dKey):
		"""
		npGetFieldDatatype() -- retrieves a string datatype label for the keyword
		identified by dKey
		"""
		if dKey in self.tDict:
			keywordDatatype = self.tDict[dKey].ripGetDatatype()
		else:
			keywordDatatype = 'no-type-assigned'
		return keywordDatatype
		# end::kwdGetDatatype()

	def kwdSetValue(self,dKey,newValue,keyValueType=None):
		"""
		kwdSetValue() assigns any Python variable bound to name (dKey)
		This function also optionally lets the user assign (or override a pre-assigned)
		datatype using 'keyValueType' to the variable, if a datatype has not yet been assigned.
		An optional permissive assignment protocol is now supported where if the name (dKey)
		does not already exist in this RipObj's namespace, it is added on-the-fly
		after which the assignment of newValue is performed. This facilitates
		an economical usage which omits an explicit prior r.kwdCreate(name) call to allow 
		r.kwdSetValue(name,value) to directly create 'name', then assign value in
		one call.
		
		Alternative Example:
			Example: r3 = RipMgr('name')
				    r3['var1'] = 345 to assign,
				    getMyValue = r3['var1']
		"""
		# check if a boolean submitted as a string, convert to a Python boolean if it is
		if newValue is 'True' or newValue is 'TRUE':
			newValue = True
		if newValue == 'False' or newValue is 'FALSE':
			newValue = False

		# use new permissive semantics where we create a new dKey instance
		# if one did not already exist
		if dKey not in self.tDict:
			self.kwdCreate(dKey)
		#
		try:
			self.tDict[dKey].ripSetValue(newValue)
		except UnboundLocalError as e:
			print("kwdSetValue: error (%s)"%(e))
		except KeyError as e:
			print("kwdSetValue: error (%s)"%(e))
		# next, set the datatype if possible
		if keyValueType is not None:
			self.tDict[dKey].ripSetDatatype(keyValueType)
		# end::kwdSetValue()

	def kwdGetValue(self,dKey):
		"""
		kwdGetValue(dKey) retrieves a value for a variable or keyword identified by
		a string name (e.g. dKey). If no value exists for this name, a None is returned.
		Alternative Example using dictionary style notation on RipMgr instance:
			Example: r3 = RipMgr('InputVariables')
				r3['var1'] = 345 ...use this to assign a value...
				getMyValue = r3['var1']   ...use this to retrieve a value
		"""

		'''
		2017_12_13. Ted: guessing "keyValues" should be "keysValue," meant to give the
		assigned var inside try, def-wide scope.

		'''
	#	keyValues = None
		keysValue = None

		if dKey not in self.tDict:
			print("kwdGetValue: warning! key(%s) is missing from this RipMgr instance"%(dKey))
		#
		try:
			keysValue = self.tDict[dKey].ripGetValue()
		except KeyError as e:
			print("kwdGetvalue missing-key-error %s"%(e))
		except UnboundLocalError as e:
			raise "kwdGetValue unbound local error", e
		#
		return keysValue
		# end::kwdGetValue()

	def kwdSetDefault(self,dKey,newDefault,keyDefaultType=None):
		"""
		kwdSetDefault() assigns an objects default value
		"""
		if dKey in self.tDict:
			self.tDict[dKey].ripSetDefault(newDefault)
		if keyDefaultType != None:
			self.tDict[dKey].ripSetDatatype(keyDefaultType)
		# end::kwdSetDefault()

	def kwdGetDefault(self,dKey):
		"""
		kwdGetDefault() returns an objects default value as assigned earlier.
		"""
		if dKey in self.tDict:
			defaultValue = self.tDict[dKey].ripGetDefault()
		else:
			defaultValue = None
		return defaultValue
		# end::kwdGetDefault()

	def kwdSetMinimum(self,dKey,newMin):
		"""
		"""
		if dKey in self.tDict:
			self.tDict[dKey].ripSetMinValue(newMin)
		# end::kwdSetMinimum()
		
	def kwdGetMinimum(self,dKey):
		"""
		kwdGetMinimum() retrieves a user assigned valid minimum associated with
		this symbol.
		"""
		thisMin = None
		if dKey in self.tDict:
			thisMin = self.tDict[dKey].ripGetMinValue()
		#
		return thisMin 
		# end::kwdGetMinimum()
	
	def kwdSetMaximum(self,dKey,newMax):
		"""
		kwdSetMaximum() assigns a user defined maximum value to a symbol, assuming it
		is a numeric quantity. For strings, this may be the lexically highest ordered
		entity.
		"""
		if dKey in self.tDict:
			self.tDict[dKey].ripSetMaxValue(newMax)	
		# end::kwdSetMaximum()
		
	def kwdGetMaximum(self,dKey):
		"""
		kwdGetMaximum() retrieves a user assigned valid maximum associated with
		this symbol.
		"""
		thisMax = None
		if dKey in self.tDict:
			thisMax = self.tDict[dKey].ripGetMaxValue()
		return thisMax

		# end::kwdGetMaximum()*/
		
	def kwdSetDescript(self,dKey,newDescrip):
		"""
		kwdSetDescript() --assigns a long descriptive string to a symbol
		"""
		if dKey in self.tDict:
			self.tDict[dKey].ripSetDescript(newDescrip)
		# end::kwdSetDescript()
		
	def kwdGetDescript(self,dKey):
		"""
		kwdGetDescript() retrieves a long descriptive string assoc with this symbol as
		referenced by its name.
		"""
		if dKey in self.tDict:
			desc = self.tDict[dKey].ripGetDescript()
		else:
			desc = None
		return desc
		# end::kwdGetDescript()

	def kwdRemove(self,dKey):
		"""
		kwdRemove() -- removes keyword or symbol identified by dKey
		"""
		if dKey in self.tDict:
			self.tDict.pop(dKey)
		msg = 'Entry named %s has been removed' % (dKey)
		return msg
		# end::kwdRemove()

	def kwdSetOrdpos(self,dKey,newOrdPos):
		"""
		npSetOrdpos() -assigns ordinal position of this keyword in the sequence of creation
		"""
		assert(newOrdPos >= 0 and newOrdPos <= RipMaxKeyWords)
		if dKey in self.tDict:
			self.tDict[dKey].ripSetOrdPos(newOrdPos)
		# end::kwdSetOrdpos()
		
	def kwdGetOrdpos(self,dKey):
		"""
		kwdGetordpos() --returns ordinal position of this keyword in the sequence of creation
		If no ordinal positional has been assigned this returns Pythons None.
		"""
		ordPos = None
		if dKey in self.tDict:
			ordPos = self.tDict[dKey].ripGetOrdPos()
			assert(ordPos >= 0 and ordPos <= RipMaxKeyWords)
		return ordPos
		# end:: kwdGetOrdpos()

	def kwdSetLength(self,dKey,newLen):
		"""
		kwdSetLength() assign a (char) length to a keywords value entity, if it matters
		This is tricky in a dynamic language, since 'length' can be ambiguous for types others
		then a char or string.
		"""
		if dKey in self.tDict:
			self.tDict[dKey].ripSetLength(newLen)
		# end::pwdSetLength()
		
	def kwdGetLength(self,dKey):
		"""
		npGetLength() --return length of keywords value item, where it matters.
		If no length has been assigned this returns None.
		"""
		thisLen = None
		if dKey in self.tDict:
			thisLen = self.tDict[dKey].ripGetLength()
		return thisLen
		# end:: npGetLength()
		
	def kwdSetPrecision(self,dKey,newPrecision):
		"""
		kwdSetPrecision() assigns a precision for formatting a RipObj scalar value as a string.
		   using kwdFormatAsStr()
		"""
		if dKey in self.tDict:
			self.tDict[dKey].ripSetPrecision(newPrecision)
		# end:: kwdSetPrecision()
		
	def kwdGetPrecision(self,dKey):
		"""
		kwdGetPrecision() retrieves the precision property for a RipObj scalar value, used for
		construction various string formatting transformations. kwdFormatAsStr() uses this
		automatically along with the variables length property.
		"""
		thisPrec = None
		if dKey in self.tDict:
			thisPrec = self.tDict[dKey].ripGetPrecision()
		return thisPrec
		
	def kwdSetItem(self,keyWord,itemIndex,newValue):
		"""
		kwdSetItem() --for key keyWord, assign a newValue at itemIndex. Normally the first item
		is submitted using kwdAddItem() and subsequent updates are handled via this function
		"""
		try:
			if not keyWord in self.tDict:
				print('kwdSetItem::Cannot add %s - for nonexistant key %s!' % (repr(newValue),keyWord))
			else:
				self.tDict[keyWord].ripSetItem(itemIndex,newValue)
		except (KeyError,TypeError) as e:
			print("kwdSetItem error (%s)"%(e))
		#
		# end::kwdSetItem()

	def kwdGetItem(self,keyWord,itemKey):
		"""
		kwdGetItem() --retrieves an item identified by itemKey owned by keyWord
		"""
		if keyWord in self.tDict:
			child = self.tDict[keyWord]
			if TestIfIter(child.ripoObj):
				keywordItem = self.tDict[keyWord].ripGetItem(itemKey)
		else:
			keywordItem = 'none'
		return keywordItem
		# end::kwdGetItem()

	def kwdGetItemCount(self,keyWord):
		resultCount = -1
		if keyWord in self.tDict:
			child = self.tDict[keyWord]
			if TestIfIter(child.ripoObj):
				resultCount = len(self.tDict[keyWord].ripoObj)
		return resultCount
		# end::kwdGetItemCount()

	def kwdFetchItems(self,keyWord):
		"""
		kwdFetchItems() returns a (list) of multi-valued items associated with a keyword
		whose role is 'multiple' and where all values are typically the same datatype.
		Note to test: this may require a shallow copy.copy() or a deep one copy.deepcopy()
		to work right - tbd.
		"""
		keyItems = None
		if keyWord in self.tDict:
			kwObj = self.tDict[keyWord]
			if TestIfIter(kwObj.ripoObj):
				keyItems = kwObj.ripoObj
		else:
			keyItems = []
		# return intact list keyItems[]
		return keyItems
		# end::kwdFetchItems()

	def kwdRemoveItem(self,keyWord,itemKey):
		"""
		kwdRemoveItem() remove an item at itemKey. This is expensive but useful when only a small
		number of items are to be removed.
		"""
		if keyWord in self.tDict:
			child = self.tDict[keyWord]
			if TestIfIter(child.ripoObj):
				self.tDict[keyWord].ripRemoveItem(itemKey)
			else:
				print("kwdRemoveKeywordItem() cannot operate on non-iterables!")
		# end::kwdRemoveKeywordItem()

	def kwdClearItems(self,dKey):
		"""
		kwdClearItems() -- for dKey, if a list, this zeros out list contents but leaves an empty list.
		"""
		for dKey in self.tDict:
			if TestIfIter(self.tDict[dKey].ripoObj):
				self.tDict[dKey].ripoObj = []
				self.tDict[dKey].ripoCount = 0 
		# end::kwdClearItems()

	def kwdSetFilter(self,dKey,whichState=None):
		"""
		kwdSetFilter() -- assigns a filter state (True|False) to this RipObjinstance.
		if whichState is not specified, the default is to set it True. 
		A RipObj's logic state is used to control selection-set-membership, various
		actions like which objects are pickled, or in some cases, masking a subset of
		RipObj.
		"""
		boolState = True
		if not whichState is None and (whichState is True or whichState is False):
			boolState = whichState
		#
		if dKey in self.tDict:
			self.tDict[dKey].ripSetFilter(boolState)
		
	def kwdGetFilter(self,dKey):
		"""
		kwdGetFilter() -- retrieves the logic (filter) state of a RipObj instance.
		"""
		boolResult = False
		if dKey in self.tDict:
			boolResult = self.tDict[dKey].ripGetFilter() 
		return boolResult
		# end::kwdGetFilter()
		
	# ---------------------------------------------------------------------
	# mgr prefaced functions serve higher level RipMgr() functionality
	# ---------------------------------------------------------------------
	def mgrSetUUID_DNS(self,newDNS):
		"""
		mgrSetUUID_DNS() -- assign a fully qualified domain bname (FQDN) such as
		www.nasa.gov or www.lupinelogic.com. Note that a URL typically adds one or more
		link qualifier to a FQDN, e.g. "www.mywebsite.com/this/is/a/url/"
		"""
		self.MgrUUID_DNS = newDNS
	# end::mgrSetUUID()

	def mgrGetUUID_DNS(self):
		"""
		mgrgetUUID_DNS -- retrieve the user defined fully qualified domain name (e.g.
		a FQDN) or absolute domain name,(e.g. www.lupinelogic.com), used to
		make the UUID when uuid.NAMESPACE_DNS qualifier is used. Note this is different
		from uuid.NAMESPACE_URL, which is a locator based input
		"""
		return self.MgrUUID_DNS
	# end:: mgrGetUUID_DNS()

	def mgrCreateUUID(self,seedInfo=None):
		"""
		Inputs: seedInfo --optional additional text material that is appended to
		the uuid.NAMESPACE_DNS part of the seed, to make the UUID.
		mgrCreateUUID() --uses a standardized internal DNS based recipe
		to create a UUID for this RipMgr instance.
		# DOC Extract: make a UUID using a SHA-1 hash of a namespace UUID and a name
		  >>> uuid.uuid5(uuid.NAMESPACE_DNS, 'python.org')
		  UUID('886313e1-3b8a-5372-9b90-0c9aee199e5d')
		"""
		# we also add a seconds count to further distiguish the uniqueness of this
		# uuid, based on uuid5 which is a SHA-1 hash...
		thisTimeStr = str(time.time())
		if seedInfo is not None:
			uniqueInfo = "%s_%s"%(repr(seedInfo),thisTimeStr)
		else:
			uniqueInfo = "%s"%(thisTimeStr)
		self.MgrUUID = uuid.uuid5(uuid.NAMESPACE_DNS,self.MgrUUID_DNS+uniqueInfo)
		return self.MgrUUID
		# end:: mgrCreateUUID()
		
	def mgrSetUUID(self,usersPredefinedUUID):
		"""
		mgrSetUUID -- assign RipMgr()'s uuid based on a users predefined UUID
		that was constructed elsewhere. it MUST be in the same form as the
		UUID that would have been constructed here, in a Python UUID object.
		"""
		self.MgrUUID = usersPrefinedUUID
		# end::mgrSetUUID()
		
	def mgrGetUUID(self):
		"""
		mgrGetUUID() -- fetch a printable UUID as held in MgrUUID. When stored as a
		official Python uuid object, the str() conversion "knows" how to correctly
		dereference this into a printable UUID, same as self.MgrUUID.get_hex() would.
		"""
		printableUUID = str(self.MgrUUID)
		return printableUUID
	
	def mgrGetVersion(self):
		"""
		mgrGetVersion() -- returns full version release tag string for this RipMgr software version.
		"""
		return self.RipMgrVersion
		# end::mgrGetVersion()
		
	def mgrGetLabel(self):
		"""
		mgrGetLabel() retrieves label or name of entire RipMgr instance itself
		"""
		return self.RipMgrLabel
		# end::mgrGetLabel()
		
	def mgrGetPlatform(self):
		"""
		mgrGetPlatform() --retrieves name of platform and architecture this app is running on.
		For example: 'platform: linux2 architecture: i386'
		"""
		# arch will be added later..avail on linux via uname -a etc etc.
		resultStr = 'platform: ' + sys.platform 
		return resultStr
		# end::mgrGetPlatform()

	def mgrLogOpen(self,logPath):
		"""
		mgrLogOpen() -- opens a session log file, named by logPath.
		"""
		self.MgrLogPath = ""
		try:
			self.MgrLogf = open(logPath,"wt")
		except:
			print("mgrLogOpen could not open log(%s)"%(logPath))
		# only store this path on success
		self.MgrLogPath = logPath
		# end::mgrLogOpen()

	def mgrLog(self,msg):
		"""
		mgrLog() -- writes a message to the session log file
		"""
		if self.MgrLogDest == "both":
			print(msg)
			self.MgrLogf.write(msg+'\n')
		else:
			self.MgrLogf.write(msg+'\n')
		# end::mgrLog()

	def mgrLogClose(self):
		"""
		mgrLogClose() --closes the session log file
		"""
		try:
			self.MgrLogf.close()
		except:
			print("mgrLogClose() could not close(%s)"%(self.MgrLogPath))
		# end::mgrLogClose()
		
	def mgrFetchIterator(self):
		"""
		mgrFetchIterator() -- returns a Python iterator object supporting
		the usual 'next' semantics, allowing caller to step linearly through the
		set of tDict{} elements one by one
		"""
		keyIter = iter(self.tDict)
		return keyIter
		# end::mgrFetchIterator()

	def mgrListAllVisible(self):
		"""
		mgrListAllVisible() -- generates a list of all RipObj instances whose logic filter state is True
		"""
		selSet = []
		# there is a good map-reduce (list comprehension) way to do this as well...
		for dKey in self.tDict:
			if self.tDict[dKey].ripGetFilter() is True:
				selSet.append(dKey)
		return selSet
		# end::mgrListAllVisible()
	
	def mgrSetAllFilters(self,whichState=None):
		"""
		mgrSetAllFilters() assigns ALL RipObj logic-filter states to either True or False
		as directed by whichState. If called without a directive (True or False) all logic 
		filter states are set to False by default.
		"""
		boolState = False
		if not whichState is None and (whichState is True or whichState is False):
			boolState = whichState
			for dKey in self.tDict:
				self.tDict[dKey].ripSetFilter(boolState)
		# end::mgrSetAllFilters()

	def mgrMaskByUserDict(self,boolDict):
		"""
		mgrMaskByUserSet() -- using a user defined dictionary where the RipObj name is the key
		and a logic filter setting (True or False) is the value, adjust the internal tDict
		logic filter settings for any RipObj named, using this {True|False} mask value.
		"""
		kTrue = 0
		for dKey in boolDict:
			if boolDict[dKey] is True:
				kTrue += 1
				self.tDict[dKey].ripSetFilter(True)
			else:
				self.tDict[dKey].ripSetFilter(False)
		# as a convenience, returns a count of those set true
		return kTrue
		# mgrMaskByUserDict()
	
	def mgrMaskByUserList(self,boolList,revSemantics=None):
		"""
		mgrMaskByUserList() -- this function by default assigns a logic filter state 
		of 'True' to any RipObj name appearing in the list (boolList) submitted. 
		Conversely, using the default semantics, any name NOT in the list (boolList) is assigned a l
		ogic filter state of False.
		
		Option: If the option revSemantics flag is set True, this reverses the default semantics as follows:
		the list of names now corresponds to all those RipObj whose logic filter state will be set to False,
		and those names NOT appearing in the boolList will have their logic filter set True.
		"""
		kTrue = 0
		if revSemantics is None:
			for dKey in self.tDict:
				if dKey in boolList:
					kTrue += 1
					self.tDict[dKey].ripSetFilter(True)
				else:
					self.tDict[dKey].ripSetFilter(False)
		else:
			for dKey in self.tDict:
				if dKey in boolList:
					kFalse += 1
					self.tDict[dKey].ripSetFilter(False)
				else:
					self.tDict[dKey].ripSetFilter(True)
			
		# as a convenience, returns a count of those set true
		return kTrue
		# mgrMaskByUserList()

	def mgrLoadItemList(self,myKey,myList,op=None):
		"""
		mgrLoadItemList() -- this is a wrapper for ripLoadItemList() that transfers a
		users list directly into ripObj.RipItems[] list.
		Inputs: myList -- user specified list
		"""
		if myKey in self.tDict:
			ripObj = self.tDict[myKey]
			# options are: 'replace' or 'append'
			ripObj.ripLoadItemList(myList,'replace')
		# end:: mgrLoadItemList()

	def mgrFetchItemList(self,myKey,op=None):
		"""
		mgrFetchItemList() -- retrieves an intact list for key myKey
		If myKey is not found, this returns None instead of a list.
		Inputs: myKey -- name of entry holding list
					 op --optional command 'sort' returns list in simple sorted order
					 (where no cmp() is used) otherwise returns list in un-sorted natural order.
		"""
		lis = None
		if myKey in self.tDict:
			ripObj = self.tDict[myKey]
			if TestIfIter(ripObj.ripoObj):
				lis = ripObj.ripGetItemList()
				# later add more explicit options for ascending,descending sort orders
				if not op is None and op == 'sort':
					lis.sort()
		return lis

	def mgrSesElapsedTime(self):
		"""
		mgrSesElapsedTime() returns the sessions elapsed time in decimal seconds since the
		initial invocation of the first RipMgr instance.
		"""
		elapSec = time.time() - self.MgrSesStrtSec
		return elapSec
		# mgrSesElapsedTime()
		
	def mgrTimeIO(self,whichOperation):
		"""
		mgrTimeIO -- general purpose IO timer. This needs to be refined, for storage
		in a named-event dictionary.
		 Inputs: operationString ('start','end')
		  when 'start' is used, a timer is started for a read or write operation
		  when 'end' is used, a timer is stopped and an elapsed time in sec is returned
		  if 'reset', the current event time is reset
		"""
		if whichOperation == 'reset':
			self.MgrElapsed_Sec	= 0.0
			self.MgrStartIO_Sec	= 0.0
			self.MgrEndingIO_Sec = 0.0
		elif (whichOperation == 'start' or whichOperation == 'begin'):
			self.MgrEndingIO_Sec = 0.0
			self.MgrElapsed_Sec = 0.0
			self.MgrStartIO_Sec = time.time()
		elif (whichOperation == 'end' or whichOperation == 'stop'):
			self.MgrEndingIO_Sec = time.time()
			self.MgrElapsed_Sec  = self.MgrEndingIO_Sec - self.MgrStartIO_Sec
			self.MgrStartIO_Sec  = 0.0
			self.MgrEndingIO_Sec = 0.0
		return self.MgrElapsed_Sec
		# end:: mgrTimeIO()

	def mgrDoesValueExist(self,checkKey):
		"""
		mgrDoesValueExist() --returns boolean as to whether a key identifying a
		'keyword' currently has a value defined (e.g. is not None)
		"""
		result = False
		if checkKey in self.tDict:
			ripObj = self.tDict[checkKey]
			if ripObj.ripoObj != None:
				result = True
			else:
				result = False
		return result
		# end::mgrDoesValueExist()

	def mgrFetchRipObj(self,objectKey):
		"""
		mgrFetchRipObj() -- retrieves a single RipObj() instance identified by
		  the objectKey passed in. This retrieval is on the basis of the variable name.
		"""
		ripObject = None  # dummy placeholder default to provide minimal existance
		if objectKey in self.tDict:
			ripObject = self.tDict[objectKey]
		return ripObject
		# end::mgrFetchObject()

	def mgrFetchObjByOrder(self,thisOrdPos):
		"""
		mgrFetchObjByOrder() -- fetch a whole RipObj() instance, by ordinal position
		Implementation is not yet as efficient as it should be.. TBD
		"""
		resultObj = None
		for d in self.tDict:
			fldObj = self.tDict[d]
			if thisOrdPos == fldObj.ripGetOrdPos():
				resultObj = fldObj
				break
		return resultObj
		# end::mgrFetchObjByOrder()

	def mgrRemoveObj(self,dKey):
		"""
		mgrRemoveObj() -- removes field identified by dKey from dictionary of table fields
		"""
		if dKey in self.tDict:
			self.tDict.pop(dKey)
		msg = 'Entry named %s has been removed' % (dKey)
		return msg
		# end::mgrRemoveObj()

	def mgrRebindObj(self,dKey,completedObj):
		"""
		mgrRebindObj() -- preemptively rebinds a presumed completed RipObj to
		an internal dictionary at entry identified by 'dKey'. In this case it does not
		matter if the completedObj was previously bound or not. This call is only very rarely
		required (and may be deprecated in the future)
		"""
		self.tDict[dKey] = completedObj
		# end::mgrRebindObj()

	def mgrGetCount(self):
		"""
		mgrGetCount() -- fetch count of the number of RipObj instances owned by this parent RipMgr()
		This is the count of keywords present in the collection.
		"""
		n = len(self.tDict)
		self.nItems = n
		return n

	def mgrGetRefCounts(self):
		"""
		mgrGetRefCounts() --retrieves reference counts of RipObj and RipMgr objects and returns
		these as a tuple (RipObjRefCount, RipMgrRefCount).
		"""
		refLis = []
		refLis.append(sys.getrefcount(RipObj))
		refLis.append(sys.getrefcount(RipMgr))
		return refLis
		# end::mgrGetRefCounts()

	def mgrListNames(self):
		"""
		mgrListNames() -- returns a sorted list of all field names or keywords present in the object
		"""
		namList = []
		# d (key) is the fieldName or Keyword identifying each RipObj() instance
		for d in self.tDict:
			namList.append(d)
		namList.sort()
		return namList
		# end::mgrListNames()

	def mgrListNameValues(self,asSorted=None):
		"""
		mgrListNameValues() -- produces a comma delimited list of all name value pairs,
		optionally sorted by name-key if asSorted argument is not None.
		"""
		lis = []
		resultLis =[]
		for d in self.tDict:
			rObj = self.tDict[d]
			rKey = d
			tupl = (rKey,rObj)
			lis.append(tupl)
		#
		if not asSorted is None:
			lis.sort()
		#
		for iObj in lis:
			rName = iObj[0]
			rObj  = iObj[1]
			outBuf = '%s,%s\n'%(rName,repr(rObj.ripGetValue()))
			resultLis.append(outBuf)
		return resultLis
		# end::mgrListNameValues()

	def mgrShowProperties(self,labels=None,sortBy=None):
		"""
		mgrShowProperties() -- returns a formatted list referencing all current field metadata entries
		labels .. if 'labels' argument passed in, produces formatted labeled output otherwise produces raw info
		'name' .. if sortby arg is 'name', sorts the output by name, or, can use 'order' to have sorted by ordinal pos
		"""
		rptLis = []
		idxLis = []
		if sortBy == None:
			sortBy = 'name'
		# NOTE 2009.01.02: this sorting logic needs a final-acceptance test and check
		# though it appears to be working ok.
		if sortBy == 'name' or sortBy == 'order':
			for d in self.tDict:
				if sortBy == 'name':
					# sort key is name here, 'd' is the key,stored at index 1.
					recd = (self.tDict[d].ripGetName(),d)
				else:
					# in this case, sort key is ordinal pos... 
					recd = (self.tDict[d].ripGetOrdPos(),d)
				idxLis.append(recd)
			# recd tuple (i,j) is now normalized such we always
			#	 sort on first element 'i' not 'j'
			idxLis.sort()
		#
		for idx in idxLis:
			# fetch the dictionary KEY out of the base0 tuples 2nd element, i=1
			dKey = idx[1]
			wrkObj = self.tDict[dKey]
			# fetch individual	interior RipObj() attributes
			nam	    = wrkObj.ripGetName()
			typ	    = wrkObj.ripGetDatatype()
			ord	    = wrkObj.ripGetOrdPos()
			fldLen  = wrkObj.ripGetLength()
			prec    = wrkObj.ripGetPrecision()
			colOrd  = wrkObj.ripGetColOffset()
			valu    = wrkObj.ripGetValue()
			# Format order is: Name,datatype,length,precision,ordinalPosition,colOffset
			if labels != None:
				fmtBuffer = 'Name %s Type %s Length %d Precision %d Order ([%s]=value)'% (nam,typ,fldLen,prec,ord,repr(valu))
			else:
				fmtBuffer = '%s %s %d %d %d %s' % (nam,typ,fldLen,prec,ord,repr(valu))
			# add latest entry to output list..
			rptLis.append(fmtBuffer)
			# return completed formatted output list to caller
		return rptLis
		# end::mgrShowProperties()

	def mgrShowSummary(self,labels=None):
		"""
		mgrShowSummary() -- returns a formatted list referencing a summary of properties for the fields or keywords
		if 'labels' argument passed in, produces formatted labeled output otherwise produces raw info.
		Note: since a populated list containing ALL values is returned, this can be memory intensive 
		for very large RipMgr collections.
		"""
		keywdLis = []
		idxLis = []
		# keyIdx here is a list, sorted by key value
		keyIdx = sorted(self.tDict)
		#
		for kIdx in keyIdx:
			wrkObj = self.tDict[kIdx]
			# fetch individual	interior RipObj() attributes
			nam	        = wrkObj.ripGetName()
			typ	        = wrkObj.ripGetDatatype
			keywdDflt   = wrkObj.ripGetDefault()
			keywdLen    = wrkObj.ripGetLength()
			keywdItemCount = wrkObj.ripGetItemCount()
			keywdItems     = wrkObj.ripGetItemList()
			keywdValu      = wrkObj.ripGetValue()
			# Format order is: Name,datatype,length,precision,ordinalPosition,colOffset
			if labels != None:
				fmtBuffer = 'Name %s Len %d Default %s Value %s ItemCount %d <%s>' %\
							(nam,keywdLen,keywdDflt,keywdValu,keywdItemCount,str(keywdItems))
			else:
				fmtBuffer = '%s %d %s %s %d <%s>' % \
							(nam,keywdLen,keywdDflt,keywdValu,keywdItemCount,str(keywdItems))
			# add latest entry to output list..
			keywdLis.append(fmtBuffer)
			# return completed formatted output list to caller
		return keywdLis
	# end::mgrShowSummary()

	def mgrShowRawSummary(self,labels=None):
		"""
		mgrShowRawSummary() -- returns a unsorted formatted list referencing a summary of
		properties for the fields or keywords
		if 'labels' argument passed in, produces formatted labeled output otherwise produces raw info.
		Note: since a populated list containing ALL values is returned, this can be memory intensive 
		for very large RipMgr collections.
		"""
		keywdLis = []
		ithName = 0
		for kName in self.tDict:
			ithName +=1
			wrkObj = self.tDict[kName]
			nam	        = wrkObj.ripGetName()
			typ	        = wrkObj.ripGetDatatype
			keywdDflt   = wrkObj.ripGetDefault()
			keywdLen    = wrkObj.ripGetLength()
			keywdItemCount = wrkObj.ripGetItemCount()
			keywdItems     = wrkObj.ripGetItemList()
			keywdValu      = wrkObj.ripGetValue()
			# Format order is: Name,datatype,length,precision,ordinalPosition,colOffset
			if labels != None:
				fmtBuffer = '%4.4d Name %s Len %d Default %s Value %s ItemCount %d <%s>' %\
				(ithName,nam,keywdLen,keywdDflt,keywdValu,keywdItemCount,str(keywdItems))
			else:
				fmtBuffer = '%4.4d %s %d %s %s %d <%s>' % \
				(ithName,nam,keywdLen,keywdDflt,keywdValu,keywdItemCount,str(keywdItems))
			# add latest entry to output list..
			keywdLis.append(fmtBuffer)
			# return completed formatted output list to caller
		return keywdLis
	# end::mgrShowRawSummary()


	def mgrResetAll(self):
		"""
		mgrResetAll() -- zero out memory used by all major dictionary elements and lists
		This is a minimal implementation, could be vastly improved to cover all the
		other fields in RipMgr()
		"""
		self.tDict       = {}
		self.tGroupDict  = {}
		self.tStanzaDict = {}
		# end::mgrResetAll()
		
	def mgrStore(self,ripStorePath=None):
		"""
		mgrStore() -- store the entire RipMgr object in a persistent form using pickle. 
		"""
		# Note: mgrStore() is implemented now but in a minimal fashion
		# Receipe:
		# self.tDict
		# self.tGroupDict
		# self.tStanzaDict
		# md5 on tDict
		if ripStorePath == None:
			thisPath = self.PersistPath
		else:
			thisPath = ripStorePath
		#
		pkf = open(thisPath,"wb")
		pkgDic ={}
		pkgDic['tDict']      = self.tDict
		sha1Obj = hashlib.sha1(str(self.tDict))
		pkgDic['tDict.sha1']   = sha1Obj.hexdigest()
		pkgDic['tGroupDict']  = self.tGroupDict
		pkgDic['tStanzaDict'] = self.tStanzaDict
		cPickle.dump(pkgDic, pkf,-1)
		pkf.close()
		print("Sorry! mgrStore() (persistent-store) is not yet implemented in this RipMgr version")
		 # end::mgrStore()

	def mgrRestore(self,ripStorePath=None):
		"""
		mgrRestore(path) -- restore a persistent-store instance of a RipMgr object from
		an external file. The storage format currently is a standard Python pickle instance.
		"""
		# Note: mgrRestore() implemented but only minimally.
		if ripStorePath is None:
			thisPath = self.PersistPath
		else:
			thisPath = self.ripStorePath
		#
		pkf = open(thisPath,"rb")
		pkgDic = {}
		pkgDic = cPickle.load(pkf)
		pkg.close()
		self.tDict       = pkgDic['tDict']
		sha1Obj = hashlib.sha1(str(self.tDict))
		chkDigest = sha1Obj.hexdigest()
		sha1digest = pkgDic['tDict.sha1'  ]
		# sha1digest should equal chkDigest here, but we do not yet enforce this.
		# more testing is needed first.
		self.tGroupDict  = pkgDic['tGroupDict' ]
		self.tStanzaDict = pkgDic['tStanzaDict']
		return 0
		# end::mgrRestore()


	def mgrStanzaAssignKeyword(self,stanzaLabel,ripName):
		"""
		mgrStanzaAssignKeyword() -- assign ripName (a keyword) to the stanza identified by stanzaLabel
		Recall that tStanzaDict is a dictionary keyed by stanzaLabels, that contains one list,
		where the list is the collection of keywords declared within the given stanza.
		To identify a single stanza vs. a multiple-stanza RipMgr instance, consult the
		number of stanzas via: nStanzas = len(self.tStanzaDict).
		"""
		if not stanzaLabel in self.tStanzaDict:
			self.tStanzaDict[stanzaLabel] = []
		# add this keyword (ripName) to the indicated stanza
		if ripName in self.tDict:
			self.tStanzaDict[stanzaLabel].append(ripName)
		return len(self.tStanzaDict)
	# end::mgrStanzaAssignKeyword()

	def mgrStanzaGetCount(self):
		"""
		mgrStanzaGetCount() --retrieve the count of stanzas defined in
		 a given RIP text file.
		"""
		return len(self.tStanzaDict)
		# end:: mgrStanzaGetCount()

	def mgrStanzaListAll(self):
		"""
		mgrStanzaListAll() --for each stanza,
			   lists ALL names (keywords) registered to given stanza.
			Returns: a list of tuples where the first element of any tuple
			is the stanzaLabel and the 2nd tuple element is the keyword.
		"""
		listOfTuples = []
		for stanzaName in self.tStanzaDict:
			kwdLis = self.tStanzaDict[stanzaName]
			# for this stanza-label 'stanzaName' list all keywords owned by it
			for kwdName in kwdLis:
				listOfTuples.append((stanzaName,kwdName))
		# return (sorted) list of tuples
		listOfTuples.sort()
		return listOfTuples
		# end:: mgrStanzaListAll()


	def mgrStanzaListNames(self,thisStanzaLabel):
		"""
		Inputs : thisStanzaLabel --name of one stanza that exists.
		Returns: a list of keywords owned by THIS ONE stanza.
		mgrStanzaListNames() -- for one stanza identified by thisStanzaLabel
		   return a list of keyword names owned only by THIS one stanza.
		"""
		nameLis = []
		for stzNam in self.tStanzaDict:
			if stzNam == thisStanzaLabel:
				kwdLis = self.tStanzaDict[thisStanzaLabel]
				for kwdName in kwdLis:
					nameLis.append(kwdName)
				break
		# return list of keywords owned JUST by thisStanzaLabel
		return nameLis
	# end ::mgrStanzaListNames()

	def mgrStanzaExtractToRipMgr(self,thisStanza):
		"""
		mgrStanzaExtractToRipMgr() -- for label 'thisStanza', separate out all keywords
		 owned by the particular stanza, bundle these in a separate RipMgr instance
		 named using the thisStanza label itself
		"""
		# extact JUST the keywords owned by thisStanza to a list, kLis[]
		kLis = self.mgrStanzaListNames(stanzaLabel)
		# for each qualifying keyword, transfer its value to a new kwd,value entity
		# inside the new RipMgr working instance, along with ALL OTHER RELEVANT
		# properties (min,max,default,categorySet, AND METADATA!! Note that we can't just make a
		# deepcopy of origRipMgr to rWrk, since we want to filter name-value-apairs to
		# JUST those belonging to thisStanza, and orig RipMgr (self) currently stores ALL
		# stanzas name-value-pair information, as well as log info and other metadata.
		# We do need to check that non tDict{} metadata is properly included in the new
		# RipMgr instance though, and check to propogate logging dynamics too etc.
		#
		# create a new RipMgr instance to store all of this stanzas assets in
		rWrk = RipMgr(stanzaLabel)
		# iterate through this stanzas assets, copying each to new RipMgr instance
		for kwd in kLis:
			rWrk.kwdCreate(kwd)
			rWrk.kwdSetValue(kwd,self.kwdGetValue(kwd))
			rWrk.kwdSetDatatype(kwd,self.kwdGetDatatype(kwd))
			rWrk.kwdSetOrdpos(kwd,self.kwdGetOrdpos(kwd))
			rWrk.kwdSetDefault(kwd,self.kwdGetDefault(kwd))
			rWrk.kwdSetMinimum(kwd,self.kwdGetMinimum(kwd))
			rWrk.kwdSetMaximum(kwd,self.kwdGetMaximum(kwd))
			rWrk.kwdSetLength(kwd,self.kwdGetLength(kwd))
			rWrk.kwdSetPrecision(kwd,self.kwdGetPrecision(kwd))
			rWrk.kwdSetDescript(kwd,self.kwdGetDescript(kwd))
			rWrk.kwdSetFilter(kwd,self.kwdGetFilter(kwd))
			if self.kwdGetItemCount() > 1:
				rWrk.mgrLoadItemList(kwd,self.mgrFetchItemList(kwd))
			# repeat for all other relevant properties and or metadata
		#
		#return the fresh populated RipMgr instance containing ALL of thisStanza's entities
		return rWrk
		# end ::mgrStanzaExtractToRipMgr()


	def mgrAddGroup(self,newGroupLabel,ripName=None):
		"""
		Inputs: newGroupLabel -- a text string identifying this group.
		        ripName -- optional, the name of one pre-existing keyword to
		        register to this new group.
		mgrAddGroup() -- define a new grouping identifier, 'newGroupLabel' to
		    group a set of keyword (variables) within.
		    Optionally, this function can also add the name of a pre-existing variable or keyword.
		Returns: a count of all groups defined.
		"""
		if ripName is not None and ripName in self.tDict:
			self.tGroupDict[newGroupLabel] = [ripName]
		else:
			self.tGroupDict[newGroupLabel] = []
		# return count of groups now defined
		return len(self.tGroupDict)
		# end::mgrAddGroup()

	def mgrIsNameInGroup(self,groupLabel,ripName):
		"""
		mgrIsNameInGroup() -- boolean function returning true if ripName is
		already in this groups membership set, otherwise returns False.
		"""
		boolResult = False
		if groupLabel in self.tGroupDict:
			objOwnedLis = self.tGroupDict[groupLabel]
			if ripName in objOwnedLis:
				boolResult = True
		# return True if ripName is present in this groups membership set
		return boolResult
		# end ::mgrIsNameInGroup()

	def mgrAssignToGroup(self,groupLabel,ripName):
		"""
		mgrAssignToGroup() -- assigns a pre-existing keyword (variable) name to a group label.
		If the group label does not previously exist, it is created prior to the bind operation.
		This is the inverse of mgrAddGroup(), so use whichever is most convenient.
		"""
		# create the group label if it does not previously exist (case-sensitive of course)
		if not groupLabel in self.tGroupDict:
			self.tGroupDict[groupLabel] = []   # init with empty list
		# Bug: ripName can be assigned multiple times, and refine this code such that
		# a given ripName can only be assigned to a group once.
		if ripName in self.tDict:
			self.tDict[ripName].ripAssignToGroup(groupLabel)
			# check if ripName ins list[] and only add it if not yet a member
			# this allows tGroupDict[group] to store multiple names but never more than
			# one instance of a given ripName (e.g. dKey).
			if not ripName in self.tGroupDict[groupLabel] :
				self.tGroupDict[groupLabel].append(ripName)
		else:
			print("%s does not yet exist - cannot be added to group %s"%(ripName,groupLabel))
		# end::mgrAssignToGroup()

	def mgrRemoveNameFromGroup(self,groupLabel,priorName):
		"""
		mgrRemoveNameFromGroup() removes a variable(name) previously residing as a member
		of a group (groupLabel) from that group.
		"""
		if groupLabel in self.tGroupDict:
			if priorName in self.tGroupDict[groupLabel]:
				self.tGroupDict[groupLabel].remove(priorName)
		# end::mgrRemoveNameFromGroup()

	def mgrRemoveGroup(self,groupLabel):
		"""
		mgrRemoveGroup() - removes one GROUP label from the set of groups. All names
		that may have been associated with this group are immediately left without this
		group association.
		"""
		if groupLabel in self.tGroupDict:
			self.tGroupDict.pop(groupLabel)
		# end::mgrRemoveGroup()

	def mgrListGroups(self,verbose=None):
		"""
		mgrListGroups() produces a list of groups defined, and if the 2nd arg is True (verbose) also
		displays for each group the list of names bound to that group. Finally, 
		A sorted list of groups, or if verbose is True, of groups+contents, is returned in either case.
		"""
		grpIdx = self.tGroupDict.keys()
		grpIdx.sort()   # design pattern: DSU...
		grpReport = []
		#
		if verbose is None:
			for grp in grpIdx:
				grpReport.append(grp)
				print("Group %s"%(grp))
		else:
			for grp in grpIdx:
				grpContents = self.tGroupDict[grp]
				grpContents.sort()
				grpReport.append(grp+','+repr(grpContents))
				print("Group %s Members %s"%(grp,grpContents))
		# end: mgrListGroups()
		#
		# return sorted list of groups defined 
		return grpReport

# end:: RipMgr() class
# -------------------------------------------------------------------------------
# NOTE: all subsequent code is NOT defined within the RipMgr() class but is
# supporting code.

# -------------------------------------------------------------------------------
def parseLine(lineStr):
	"""
	Inputs: lineStr -- one text line of delimited tokens to parse, populating a tokenObj instance.
	
	parseLine() -- this parser logic is used in mgrParseStanzaInputs(). It categorizes one text line
	into a keyword, and an associated collection of one or more referent values (value tokens).
	
	The value tokens may themselves be delimited using blanks,commas,tabs, or vert-bar |
	symbols, where most typically, a single value token is expected. A list of multiple value
	tokens may be specified, if it is comma delimited. Lists may also be explicitly specified
	in a RIP text file using a notation like [1,2,3].  Further, a value-token may consist
	of a numeric range specifier, using the colon delimited triplet notation
	   'start:end:incrementBy' triplet, or, if not incrementBy is known, using a couplet
	 notation StartValue:EndingValue (range) couplet. No embedded blanks are allowed in either
	 range specifier form.
	 NOTE: When a conformal range specifier is found, the resulting 2 or 3 values are
	 stored as a triplet tuple, or a couplet tuple when incrementBy is missing, in the
	 object's main 'value' member, accessible via tokObj.getValue() after parsing.
	 NOTE: this parseLine() routine works very closely with classifyType().
	 NOTE: single string literals, with embedded blanks, may be specified by bracketing these
	 with double quotes; where present, the information in these strings is protected and thus
	 not subjected to further parsing.  Single quote delimiting is NOT yet supported for
	 string literals.
	 NOTE: after a valid keyword and a delimiter, a "list" may be specified using
	 comma delimiting with no embedded blanks, e.g.  1,2,3,4,5. A 'list' of values may also
	 be explicitly denoted via square brackets, as in: [1,2,3,4,5] in the RIP file.

	Returns: a populated tokenObj class instance.
	"""
	# just in case parseLine() called independently, strip off end-of-line chars
	lineStr = lineStr.strip('\n').strip('\r')
	#
	dic = {}
	aTab = chr(9) # useful for identifying tabs explicitly
	foundRange = False
	dType = 'undetermined'
	tForm = 'undetermined'
	rngTuple  =  (0,0,0) # placeholder for possible range tuple
	# valid token delimiters are blank,comma,tab or vert bar
	delimSetLis = [" ",",",chr(9),"|"]
	for d in delimSetLis:
		dic[d] = 0
	# construct a list comprised of counts of occurences of dCh
	# if the delimiter used was any one of the qualifying delimiter chars as enumerated
	# in hard wired delimSetLis[] list.
	wrkLis =[lineStr.count(dCh) for dCh in delimSetLis]
	i =0
	# transfer these counts to a dictionary whose key is the
	# delimiter in question
	for chCount in wrkLis:
		dic[delimSetLis[i]] = chCount ; i += 1
	#
	keyWdChrs = []
	# First, partition each RIP text line into a'keyword' and the rest
	# of the line-string, accumulating chars into the keyword pool,
	# and, as soon as ANY given char matches the first valid char from the
	# delimiter dictionary, stop early.
	for ch in lineStr:
		if ch in dic:
			break
		else:
			keyWdChrs.append(ch)
	# glue the keywords parts together
	keyword = "".join(keyWdChrs)
	#
	# CREATE one Name-value-pair object for parsing, named after this keyword
	tokObj = tokenObj(keyword)
	# isolate all remaining text as a residual string
	residual = lineStr[len(keyword)+1:]
	# set default (fallback) properties of name-value-pair object
	tokObj.setType('integer')
	tokObj.setForm('scalar')
	tokObj.setValue(0)
	tokObj.setCount(0)
	tokObj.setRangeLow(-1)
	tokObj.setRangeHigh(-1)
	tokObj.setRangeIncr(-1)
	#
	# On 2010-11-03-22:44:01, revised to default setCount to 0, not 1
	# if length of string 'residual' here is at least one, delay registering its
	# item count...
	if len(residual)>= 1:
		tokObj.setCount(0)
	else:
		tokObj.setCount(1)
	# check FIRST if this is a string literal and protect it if it is!
	t1 = residual.startswith('"') and residual.endswith('"')
	t2 = residual.startswith("'") and residual.endswith("'")
	# we now allow either single or double quotes to demarcate a string literal.
	if t1 or t2:
		tForm = 'literal'
		dType = 'string'
		tokObj.setCount(1)
		tokObj.setType(dType)
		tokObj.setForm(tForm)
		# if a left most tab persists from (keyword vs resid) parse, strip it off
		residual = residual.lstrip(aTab)
		tokObj.setValue(residual) # placeholder, may refine as below...
		rangeIndicator = 0
	else:
		# can now safely test if internal values represent a
		#  range-specifer statement, as indicated by presence of a ":" symbol...
		rangeIndicator = residual.count(":")

	# a 2-element or 3-element LOW:HIGH:INCR range-syntax expression found,
	# so at least, extract low and high, and extract incrementBy factor if present
	if rangeIndicator >0 and residual.count("/") ==0 and residual.count("\\") == 0:
		# first set defaults
		foundRange = True
		residual = residual.lstrip(' ').rstrip(' ')
		# CLASSIFY DATA TYPE of RANGE elements (ok to re-classify here)
		tForm,dType = classifyType(residual)
		rngLis = residual.split(":")
		n = len(rngLis)
		# print("RANGE-POST-CLASSIFY dType(%s) tForm(%s)Lis[%s]"%(dType,tForm,str(rngLis)))
		tokObj.setType(dType)
		tokObj.setForm(tForm)
		# assert: tForm should be 'range' here!
		if dType == 'integer':
			tokObj.setRangeLow(int(rngLis[0]))
			tokObj.setRangeHigh(int(rngLis[1]))
			if n >2:
				tokObj.setRangeIncr(int(rngLis[2]))
		elif dType == 'float':
			tokObj.setRangeLow(float(rngLis[0]))
			tokObj.setRangeHigh(float(rngLis[1]))
			if n >2:
				tokObj.setRangeIncr(float(rngLis[2]))
		else:
			tokObj.setRangeTuple(int(rngLis[0]),int(rngLis[1]))
			if n>2:
				tokObj.setRangeIncr(int(rngLis[2]))
		# populate the RANGE TUPLE with (low,high,incrementBy)
		rngTuple=(tokObj.getRangeLow(),tokObj.getRangeHigh(),tokObj.getRangeIncr())
		tokObj.setValue(rngTuple)
	else:
		# to get this far, object was NOT a range style parameter...
		# If this is NOT a string literal, then isolate residual info on line,
		# replacing all blanks and tabs with normalization delimiter (a comma)
		if tForm == 'literal':
			# check for embedded TAB on left edge, remove it if present
			residual = residual.lstrip(aTab)
		else:
			# if NOT a literal, then allow obj to be treated as a list, so normalize
			# the collection delimiters to comma from whatevery they might have been.
			# Note: filepaths should never have any embedded blanks in them anyway
			# so this should normally NOT affect this use-case, except where Windows
			# filespecs with embedded blanks are used--these are NOT yet supported.
			residual = residual.replace(" ",",").replace(aTab,",")
		# in case a left blank was converted to a left comma, get rid of it
		residual = residual.lstrip(",")
		# if a explicit LIST was specified, as bracketed with '[' and ']' then remove.
		if residual.startswith('[') or residual.endswith(']'):
			residual = residual.lstrip('[').rstrip(']')
		#
		# convert path as policy-flag indicates..
		if convertWinPathToUnix and residual.count("\\")>0:
			residual = residual.replace("\\","/")
	
		# count no. of commas, to infer if symbol is an enumerated list
		nComma = residual.count(",")
		# NOTE: we rarely if ever expect to see a literal list of booleans
		# so these are NOT yet supported at this time.
		if nComma >0:
			# print("Found possible comma-delim list, or tuple. No. commas = %d"%(nComma))
			tokLis = [valTok for valTok in residual.split(",")]
			# CLASSIFY DATA TYPE
			tForm,dType = classifyType(tokLis)
			#
			tokObj.setForm(tForm)
			tokObj.setType(dType)
			wrkLis = []
			if dType == 'integer':
				for i in tokLis:
					targ = i.replace(" ","")
					if len(targ)>= 1:
						wrkLis.append(int(targ))
				tokLis = wrkLis
			elif dType == 'float':
				for x in tokLis:
					targ = x.replace(" ","")
					if len(targ)>= 1:
						wrkLis.append(float(targ))
				tokLis = wrkLis
			else:
				# at this point, treat list of items as a list of STRINGS
				for x in tokLis:
					targ = x.replace(" ","")
					if len(targ)>= 1:
						wrkLis.append(targ)
				tokLis = wrkLis
			# the COUNT variable stores the number of value-tokens held in tokLis[]
			# which is most often 1, but may be >1 where multiple values were encountered
			tokObj.setCount(len(tokLis))
			tokObj.setValue(tokLis)
		else: # to get here, likely we have a scalar, or
			# it might be a FILEPATH or a dict or class
			# CLASSIFY DATA TYPE
			tForm,dType = classifyType(residual)
			#
			# For scalar of dType 'boolean' process special case where
			# residual as a string was {'True'|'Enable'|'Yes'} convert string
			# to a Python boolean {True|False}on special case of dType as 'boolean'
			if dType == 'boolean' and type(residual)==type(str()) :
				tokn  = residual.upper()
				tForm = 'scalar'
				residual = False
				if(tokn =='TRUE' or tokn =='ENABLE' or tokn =='YES'):
					residual = True
				elif(tokn == 'T' or tokn =='Y'):
					residual = True
			#
			# assign rest of essential properties
			tokObj.setCount(1)
			tokObj.setForm(tForm)
			tokObj.setType(dType)
			# print("Setting type(%s) and setValue(%s)"%(dType,repr(residual)))
			# finally, assign the value back to the object
			tokObj.setValue(residual)
	# return populated tokenObjct
	return tokObj
	# end:: parseLine()

# -------------------------------------------------------------------------------
def mgrParseStanzaInputs(rtpPath,parentName=None,manifestPath=None):
	"""
	Inputs: rtpPath -- path of stanza'ized file to parse
	        parentName -- optional, name of a prior defined RipMgr instance, where if
	         no parent is specified, data is loaded into a default RipMgr() instance.
	        manifestPath --optional, name of controlled vocabulary file if available.
	mgrParseStanzaInputs() --this method parses a traditional vertically oriented 
	(e.g. stanza) keyword name-value pair text file structure.
	Lines beginning with the '#' symbol are ignored as comment lines. A dictionary entry is
	built for each valid variable found.
	Future: we may add a second runtime (manifest) list, 'KeywordControlledVocab.cv' that
	contains a full list of potential candidate variables to expect in this file,
	to facilitate a check that each and every variable expected and needed in
	a given context is actually present in the runtime input parameter file.
	"""
	manifestLis = []
	#
	# instantiate the single RipMgr() object that will hold all parsed information
	if not parentName is None:
		ripParent = RipMgr(parentName)
	else:
		ripParent = RipMgr("DefaultRipMgrParent")
	#
	if manifestPath is not None:
		manifestLis  =[]
		if os.path.exists(manifestPath):
			manif = open(manifestPath,'rt')
			manifestLis = [L.strip('\n') for L in manif]
			manifestLis.sort()
			manif.close()
		else:
			print("A manifest path (%s) was declared but could NOT be found!"%(manifestPath))
			sys.exit(-1)
	#
	if os.path.exists(rtpPath):
		inf = open(rtpPath,'rt')
	else:
		print('ERROR! parseRuntimeInputs: cannot find or open runtime input parameter file(%s)' % (rtpPath))
		sys.exit(-1)
	# 2) read whole input file into a buffer and parse out all variables
	runtimeInputLis = [ L.strip('\n') for L in inf ]
	if verbose:
		print("Input stanza(%s) contains %d lines"%(rtpPath,len(runtimeInputLis)))
	# close primary RIP input file itself
	inf.close()
	stanzaCount = 0 
	ithLine = 0
	# should name this in future the same as the parent or default RipMgr() instance
	tok = tokenObj('test')
	# NOTE: this parser intentionally employs a 2-pass scan, where the first scan
	# enumerates the RIP file's macro level (block, or multi-stanza-block) organization
	# via construction of a 'stanza' category.
	# NOTE: first do a whole list scan for END_{BlockTokens} and add any found to the
	# tStanzaDict
	ripParent.tStanzaDict = {}
	#
	for rLine in runtimeInputLis:
		rLine = rLine.strip(" ")
		if not rLine.startswith("#") and len(rLine) >= 1:
			# upper case to normalize, in case user typed "End_" or "end_"
			testBlockName = rLine.upper()
			# in future, may allow "{stanzalabel}_END" permissive syntax as well
			if testBlockName.startswith("END_"):
				# register this Blocktoken to the candidate qualifying dict
				offset = len("END_")
				# fetch and store the STANZA LABEL as blockLabel
				blockLabel = testBlockName[testBlockName.find("END_")+offset :]
				# store normalized block label in dictionary, then whenever a keyword
				# matches THIS block label, we suspect it is the block-intro label
				ripParent.tStanzaDict[blockLabel.upper()] = []
			# end:: if testBlockName END was found
		# end if a qualifying token
	# store count of discrete STANZAs found in this RIP text file
	nBlocks = len(ripParent.tStanzaDict)
	#
	if verbose:
		print("mgrParseStanzaInputs:: %d candidate block labels were registered"%(nBlocks))
		ithBlock =0
		if nBlocks >0:
			for bLbl in ripParent.tStanzaDict:
				ithBlock += 1
				print("%3d block-label (%s)"%(ithBlock,bLbl))

	stanzaDefaultName = "default_stanza"
	currentStanza = defaultStanzaName
	#
	# iterate through all text lines in RIP file...
	for rLine in runtimeInputLis:
		ithLine += 1
		rLine = rLine.strip(' ')
		# Linux requires stripping off \r, Windows requireds stripping off \n
		rLine = rLine.strip('\n').strip('\r')
		if not rLine.startswith('#') and len(rLine) >= 1 :
			valu = 'no-value'
			# return a class instance of token..
			tok = parseLine(rLine)
			keyword = tok.getKeyword()
			# print("DIAG-TOKEN-COUNT key(%s) tok-count %d"%(keyword,tok.getCount()))
			# only parse and store keywords that are NOT Block or GROUP labels
			# if this keyword matches a STANZA-START token, set the currentStanza marker to this
			if keyword in ripParent.tStanzaDict:
				currentStanza = keyword
			# keyword was NOT a STANZA-START token, so assign current stanza to fallback default.
			elif keyword.count(currentStanza)>0 and keyword.count("END_")>0:
				currentStanza = defaultStanzaName
			else:  # at this point, keyword is a 'regular' keyword, and was
				# not a STANZA-START or STANZA-END keyword
				# assert: a valid keyword was encountered, so create, register it
				ripParent.kwdCreate(keyword)
				count   = tok.getCount()
				tForm   = tok.getForm()
				dType   = tok.getType()
				loRng   = tok.getRangeLow()
				hiRng   = tok.getRangeHigh()
				incrRng = tok.getRangeIncr()
				valu    = tok.getValue()
				# print("DIAGCHK Key(%s)Value(%s) Type(%s)"%(keyword,str(valu),dType))
				ripParent.mgrStanzaAssignKeyword(currentStanza,keyword)
				ripParent.kwdSetDatatype(keyword,dType)
			#
			if valu is not 'no-value':
				ripParent.kwdSetValue(keyword,valu)
			else:
				if verbose:
					print("DIAG: key(%s) correctly ignored (%s)"%(keyword,valu))
			#
	# return fully populated RipMgr
	return ripParent
	# end::mgrParseStanzaInputs()

# -------------------------------------------------------------------------------
def mgrParseCsvInputs(rtpPath,manifestPath=None):
	"""
	mgrParseCsvInputs() --this method parses a traditional CSV style horizontally oriented 
	input file where the first line contains a comma delimited list of column headings,
	and subsequent lines contain a matching number of tokens per line.
	Lines beginning with the '#' symbol are ignored as comment lines. A dictionary entry is
	built for each valid variable found.
	NOTE: when completge, list 'lisMgr[] stores a collection of RipMgr instances, one per CSV line.
	FUTURE: we may add a second runtime (manifest) list, 'KeywordManifest.cfg' that contains
	 a full list of potential candidate variables to expect in this file, to facilitate a check
	 that each and every variable expected and needed in a given context is actually present
	 in the runtime input parameter file.
	"""
	manifestLis = []
	if not manifestPath is None:
		if os.path.exists(manifestPath):
			manif = open(manifestPath,'rt')
			manifestLis = [L.strip('\n').strip('\r') for L in manif]
			manifestLis.sort()
			manif.close()
	#
	if os.path.exists(rtpPath):
		try:
			inf = open(rtpPath,'rt')
		except (OSError,IOError) as e:
			print("mgrParseCsvInputs: (%s) cannot open(%s)"%(e,rtpPath))
			return -1
	else:
		print('ERROR! mgrParseCsvInputs: cannot find or open input CSV(%s)' % (rtpPath))
		sys.exit(-1)
	#
	# 2) read whole input file into a buffer and parse out all variables
	runtimeInputLis = [ L.strip('\n').strip('\r') for L in inf ]
	inf.close()
	stanzaCount = 0 
	ithLine = 0
	#
	# if true, expands compound value chains into parsable tokens
	# if false,ignores embedded delimiters and keeps metatokens intact
	compoundStrictExpand = True
	# if True, we expect first line to contain column label headings, one per column
	headerLabelPresent = True
	nColumns = -1
	nTokens  =  0
	# list lisMgr[] stores k-columns worth of RipMgr objects
	lisMgr = []
	capturedColumnLabels = False
	for rLine in runtimeInputLis:
		ithLine += 1
		rLine = rLine.strip(' ')  # skip over any totally blank lines
		dataTok= []
		if rLine.startswith('#') and verbose:
			print("COMMENT [%s]"%(rLine))
		elif len(rLine)>= 1:
			ValueContent = 'no-content-assigned'
			# before splitting, separate off the first token, and save entire
			# residual string segment, to protect paths with embedded blanks in them!
			if not capturedColumnLabels:
				columnLabels = rLine.split(",")
				nColumns = len(columnLabels)
				capturedColumnLabels = True
			else:
				# print("DATALINE [%s]"%(rLine))
				# create and populate a
				# new RipMgr for each line
				UniqID = "Trial"+str(ithLine).zfill(7)
				r = RipMgr(UniqID)
				dataTok  = rLine.split(",")
				nDataTok = len(dataTok)
				if nDataTok != nColumns:
					print("mgrParseCsvInputs: Warning! N %d column-labels do NOT match %d data values"%(nColumns,nDataTok))
					sys.exit(-1)
				else:
					# print("N datacolumns %d"%(nDataTok))
					for kFld in range(nDataTok):
						r.kwdCreate(columnLabels[kFld])
						r.kwdSetValue(columnLabels[kFld],dataTok[kFld])
						# FUTURE, NEXT-TO-ADD: r.kwdSetDatatype() etc and other properties
				#
				lisMgr.append(r)

	# return list of k columns worth of RipMgr instances
	return lisMgr
	# end::mgrParseCsvInputs()

# -------------------------------------------------------------------------------
def ripMgrStanzaTest(inputStanzaPath):
	"""
	ripMgrStanzaTest() -- exercises many parsing use-cases for stanza oriented 
	input parameter files.
	--Single Stanza use-case: loads all variables in a RIP file into one RipMgr instance.
	--Multi-Stanza use-case : for each stanza, populates a separate RipMgr instance as part
	  of a test sequence.
	"""
	r = mgrParseStanzaInputs(inputStanzaPath)
	if doMultiStanzaTest:
		nStanzas = r.mgrStanzaGetCount()
		if nStanzas > 1:
			"""
			de-reference this as a multi-stanza test file!
			"""
			for stanzaLabel in r.tStanzaDict:
				# unpack ALL the keywords owned by this stanza
				# into a stanza-wise RipMgr instance, named after the stanza itself
				rWrk = RipMgr(stanzaLabel)
				kLis = r.mgrStanzaListNames(stanzaLabel)
				# for each keyword, transfer its value to a new kwd,value entity
				# inside the new RipMgr working instance
				for kwd in kLis:
					valu = r.kwdGetValue(kwd)
					rWrk.kwdSetValue(kwd,valu)
					print("StanzaLabel (%s) KWD-TRANSFER of %s"%(stanzaLabel,kwd))
			#
	if doShowPropertiesTest:
		lisProp = r.mgrShowProperties()
		r.mgrLogOpen("rpShowProp.log")
		for p in lisProp:
			r.mgrLog(p)
		r.mgrLogClose()

	if doShowSummaryTest:
		r.mgrLogOpen("rpShowRawSummary.log")
		ssLis= r.mgrShowRawSummary()
		for lss in ssLis:
			r.mgrLog("raw-summary %s"%(repr(lss)))
		r.mgrLogClose()

	if doFullDiagTest:
		r.mgrTimeIO('start')
		r.mgrLogOpen("ripmgr.log")
		r.mgrLog("mgrListNameValues")
		lisAll = r.mgrListNameValues()
		for item in lisAll:
			r.mgrLog("mgrListNameValues %s"%(repr(item)))
		r.mgrLog("mgrGetRefCounts()")
		refCountLis = r.mgrGetRefCounts()
		for refCnt in refCountLis:
			r.mgrLog("mgrGetRefCounts %s"%(repr(refCnt)))
		#
		r.mgrLog("mgrShowProperties")
		lisProp = r.mgrShowProperties()
		for lp in lisProp:
			r.mgrLog("mgrShowProperties %s"%(repr(lp)))
		#
		r.mgrLog("mgrShowSummary")
		lis= r.mgrShowSummary()
		for lss in lis:
			r.mgrLog("mgrShowProperties %s"%(repr(lss)))
		elapSec = r.mgrTimeIO('stop')
		r.mgrLog("ripMgrStanzaTest required %.2f sec"%(elapSec))
		sLis = r.mgrListNamesInStanza()
		# dump list-of-tuples (stanzaID,keywordID)
		for sItem in sLis:
			tupl = sItem
			r.mgrLog("ripMgrStanzaTest  stanza %s name %s"%(tupl[0],tupl[1]))
		# close out the log
		r.mgrLogClose()
		return 0
	# end::ripMgrStanzaTest()


# -------------------------------------------------------------------------------
def ripMgrCsvTest(inputCsvPath):
	"""
	mgrCsvTest() -- exercise CSV hortizontal style parsing using a external
	input. mgrParseCsvInputs() returns a list of k lines worth of separate RipMgr
	objects fully populated.
	"""
	rTrialLis = mgrParseCsvInputs(inputCsvPath)
	print("ripMgrCsvTest: n %d run definitions"%(len(rTrialLis)))
	#
	for rMgr in rTrialLis:
		rMgr.mgrListNames()
	#
	if doShowPropertiesTest:
		ithRun = 0 
		for rMgr in rTrialLis:
			ithRun += 1
			rMgr.mgrLogOpen("ripMgrCsvTest_ShowProp.log")
			lisProp = rMgr.mgrShowProperties()
			for lp in lisProp:
				rMgr.mgrLog("ripMgrCsvTest RunDef %4d %s"%(ithRun,repr(lp)))
			rMgr.mgrLogClose()

	if doShowSummaryTest:
		ithRun =0
		for rMgr in rTrialLis:
			ithRun +=1
			rMgr.mgrLogOpen("ripMgrCsvTest_ShowSummary.log")
			lis = rMgr.mgrShowSummary()
			for lss in lis:
				rMgr.mgrLog("ripMgrCsvTest_ShowSummary %d %s"%(ithRun,repr(lss)))

	return rTrialLis
	# end::ripMgrCsvTest()

#--------------------------------------------------------------------------
# main()
#--------------------------------------------------------------------------
if __name__ == '__main__':
	"""
	main() -- test exercises for RipMgr() and RipObj() classes
	"""
	# The unit tests do not require any external input files
	if doUnitTest:
		from RipMgrUnitTest import *
		# there are currently 4 unit test sections, each testing a
		# different aspect of RipMgr functionality.
		print("UNIT-TEST version %s"%(unitTestVers))
		ripMgrUnitTest1()
		ripMgrUnitTest2()
		ripMgrUnitTest3()
		ripMgrUnitTest4()
		print("UNIT-TEST %s complete at %s"%(unitTestVers,time.ctime()))
		sys.exit(0)
	#
	if len(sys.argv)>= 2:
		inpPath = sys.argv[1]
		print('\nTest: %s started at %s' % (sys.argv[0],time.ctime()))
		print('RipMgr version %s' % (modVers))
		print("Test input file %s"%(inpPath))
		if inpPath.endswith('csv'):
			ripMgrCsvTest(inpPath)
		else:
			ripMgrStanzaTest(inpPath)
		# 
	else:
		print("RipMgr returning -- no external parsing tests specified")
	# done
