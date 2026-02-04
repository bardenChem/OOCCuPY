#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = EnergyRefinement.py

#================================
import pymp
from .commonFunctions import *
from pMolecule import *
from pMolecule.QCModel import *

from .MopacQCMMinput import MopacQCMMinput
import os, glob, sys, shutil
import numpy as np 

from pSimulation import *
from .QuantumMethods import *
#================================
#**********************************************************
class EnergyRefinement:
	'''
	Energy calculations for a set of structures using several QC methods and or External softwares
	'''
	def __init__(self,_refSystem,_trajFolder,_outFolder,_dims,_chg,_mult):
		'''
		Default constructor.
		Parameters:
			_refSystem : reference molecular information; System pDynamo class instance
			_trajFolder: folder path of the structures; string or path
			_outFolder : folder path where the results will be written; string or path
			_dims      : reaction coordinates size; list of integers
			_chg       : reference QC region charge; integer
			_multi     : reference QC region multiplicity; integer
		'''
		self.molecule 	 = _refSystem
		self.trajFolder  = _trajFolder  
		self.pureQCAtoms = []
		self.RC1 	 	 = []
		self.RC2 		 = []
		self.rc1CoordName= []
		self.rc2CoordName= []
		self.restart     = False
		self.xlen        = _dims[0]
		self.ylen        = _dims[1]
		self.charge 	 = _chg
		self.multiplicity= _mult
		self.text 		 = ""
		self.methods 	 = []
		self.fileLists   = []
		

		if hasattr(self.molecule,"qcState"): 
			self.pureQCAtoms = list(self.molecule.qcState.pureQCAtoms)
		i = 0
		self.baseName = _outFolder	
		if not os.path.exists(self.baseName): os.makedirs(self.baseName)
		if self.xlen > 1:
			_path = os.path.join( _trajFolder,"")
			self.fileLists = glob.glob(_path + "frame*.pkl")			
		elif self.xlen == 1:
			self.fileLists.append(_trajFolder+".pkl")		
		#----------------------------------------------------------------------
		if self.ylen  == 0:
			self.energiesArray = pymp.shared.array( (self.xlen) , dtype='float')
			self.indexArrayX   = pymp.shared.array( (self.xlen) , dtype='uint8')
			self.indexArrayY   = pymp.shared.array( (self.xlen) , dtype='uint8')
		else:
			self.energiesArray = pymp.shared.array( (self.xlen,self.ylen) , dtype='float')
			self.indexArrayX   = pymp.shared.array( (self.xlen,self.ylen) , dtype='uint8')
			self.indexArrayY   = pymp.shared.array( (self.xlen,self.ylen) , dtype='uint8')
		self.SMOenergies   = None
	
	
	#=====================================================================================
	def RunInternalSMO(self,_methods,_NmaxThreads):
		'''
		Run energy refinement with the semiempirical hamiltonians available wihthin pDynamo
		Parameters:
			_methods:     List of Hamiltoninas
			_NmaxThreads: Number of maximum threds to be used in the parallel section
		'''
		self.SMOenergies = {}
		self.methods 	 = _methods

		_qc_parameters = {  "active_system":self.molecule,
							"region":self.pureQCAtoms    , 
							"method_class":"SMO"         ,
							"Hamiltonian":"am1"          ,
							"multiplicity":self.multiplicity,
							"QCcharge":self.charge       }	
		#--------------------------------------------------------------------
		for smo in _methods:
			if VerifyMNDOKey(smo):
				with pymp.Parallel(_NmaxThreads) as p:
					for i in p.range( len(self.fileLists) ):
						_qc_parameters["Hamiltonian"] = smo
						qcSystem = QuantumMethods(_qc_parameters)
						qcSystem.Set_QC_System()	
						qcSystem.system.coordinates3 = ImportCoordinates3( self.fileLists[i], log=None )
						lsFrames= GetFrameIndex(self.fileLists[i][:-4])						
						if self.ylen > 0:
							try:  self.energiesArray[ lsFrames[0], lsFrames[1] ]   = qcSystem.system.Energy(log=None)
							except: self.energiesArray[ lsFrames[0], lsFrames[1] ] = self.energiesArray[0,0] + 1000
							self.indexArrayX[ lsFrames[0], lsFrames[1] ] = lsFrames[0]
							self.indexArrayY[ lsFrames[0], lsFrames[1] ] = lsFrames[1]							
						else:
							try: 	self.energiesArray[ lsFrames[0] ] = qcSystem.system.Energy(log=None)
							except: self.energiesArray[ lsFrames[0] ] = self.energiesArray[0] + 1000
							self.indexArrayX[ lsFrames[0] ] 		  = lsFrames[0]	
				#-----------------------------------------
				if self.ylen > 0:
					self.SMOenergies[smo] = self.energiesArray
					self.energiesArray = pymp.shared.array( (self.xlen,self.ylen) , dtype='float')	
				else:
					self.SMOenergies[smo] = self.energiesArray
					self.energiesArray = pymp.shared.array( (self.xlen) , dtype='float')		
			else:
				continue		
	#====================================================
	def RunInternalDFT(self,_functional,_basis,_NmaxThreads):
		'''
		Run energy refinement with the semiempirical hamiltonians available wihthin pDynamo
		Parameters:
			_methods:     List of Hamiltoninas
			_NmaxThreads: Number of maximum threds to be used in the parallel section
		'''
		self.SMOenergies = {}		
		self.methods.append(_functional)

		converger      = DIISSCFConverger.WithOptions( densityTolerance = 1.0e-10, maximumIterations = 550 )
		gridIntegrator = DFTGridIntegrator.WithOptions( accuracy = DFTGridAccuracy.Medium, inCore = True )
		qcModel        = None
		NBmodel        = self.molecule.nbModel

		if _functional == "hf": 
			qcModel = QCModelDFT.WithOptions(converger=converger,functional="hf", orbitalBasis=_basis)
		else                  :
			qcModel = QCModelDFT.WithOptions(converger=converger,functional=_functional,gridIntegrator=gridIntegrator, orbitalBasis=_basis)
			
		self.molecule.electronicState = ElectronicState.WithOptions(charge = self.charge)
		self.molecule.DefineQCModel( qcModel, qcSelection=Selection(self.pureQCAtoms) )		
		self.molecule.DefineNBModel( NBmodel )		
		
		#------------------------------------------------------------------------------
		with pymp.Parallel(_NmaxThreads) as p:
			for i in p.range( len(self.fileLists) ):				
				self.molecule.coordinates3 = ImportCoordinates3( self.fileLists[i],log=None )
				lsFrames= GetFrameIndex(self.fileLists[i][:-4])						
				if self.ylen > 0:
					self.energiesArray[ lsFrames[0], lsFrames[1] ] = self.molecule.Energy()
					self.indexArrayX[ lsFrames[0], lsFrames[1] ] = lsFrames[0]
					self.indexArrayY[ lsFrames[0], lsFrames[1] ] = lsFrames[1]
				else:
					self.energiesArray[ lsFrames[0] ] = self.molecule.Energy()
					self.indexArrayX[ lsFrames[0] ] = lsFrames[0]	
			#-----------------------------------------
			if self.ylen > 0:
				self.SMOenergies[ self.methods[0] ] = self.energiesArray
				self.energiesArray = pymp.shared.array( (self.xlen,self.ylen), dtype='float')	
			else:
				self.SMOenergies[ self.methods[0] ] = self.energiesArray
				self.energiesArray = pymp.shared.array( (self.xlen), dtype='float')

	#====================================================
	def RunMopacSMO(self,_methods,_keyWords):
		'''
		Create input for Mopac with its available Hamiltonians enabling the QC(QM)/MM potential
		Parameters:
			_methods: List of hamiltonians available in MOPAC 
			_NmaxThreads: Number of maximum threds to be used in the parallel section
		'''
		self.SMOenergies = {}
		self.methods     = _methods
		NBmodel          = self.molecule.nbModel	
		_mopacKeys       = _keyWords
		
		mop_pars = { 
			"active_system":self.molecule   ,
			"basename":self.baseName        ,
			"cood_name":"none"              ,
			"Hamiltonian":"am1"             ,
			"QCcharge":self.molecule.electronicState.charge,
			"multiplicity":self.multiplicity,
			"keywords":_keyWords            , 
		}

		for smo in _methods:
			for i in range( len(self.fileLists) ):				
				self.molecule.coordinates3 = ImportCoordinates3(self.fileLists[i],log=None)
				mop_pars["Hamiltonian"]    = smo 
				mop_pars["cood_name"]      = self.fileLists[i]
				mop = MopacQCMMinput(mop_pars)
				mop.CalculateGradVectors()
				mop.write_input(os.path.basename(mop_pars["cood_name"]))
				mop.Execute()				
				lsFrames = []
				if self.fileLists[i] == "single.pkl": lsFrames.append(0)
				else: lsFrames = GetFrameIndex(self.fileLists[i][:-4])		
				if self.ylen > 0:
					self.energiesArray[ lsFrames[0], lsFrames[1] ] = mop.GetEnergy()
					self.indexArrayX[ lsFrames[0], lsFrames[1] ] = lsFrames[0]
					self.indexArrayY[ lsFrames[0], lsFrames[1] ] = lsFrames[1]
				else:
					self.energiesArray[ lsFrames[0] ] = mop.GetEnergy()
					self.indexArrayX[ lsFrames[0] ] = lsFrames[0]					
			#----------------			
			if self.ylen > 0:
				self.SMOenergies[smo] = self.energiesArray
				self.energiesArray    = pymp.shared.array( (self.xlen,self.ylen) , dtype='float')	
			else:
				self.SMOenergies[smo] = self.energiesArray
				self.energiesArray    = pymp.shared.array( (self.xlen) , dtype='float')	
			
	#====================================================
	def RunDFTB(self,_NmaxThreads):
		'''
		Perform energy refinement using the interface available on the pDynamo with the DFTB+ software, enabling QC(QM)/MM potential.
		Parameters:
			_NmaxThreads: Number of maximum threds to be used in the parallel section
		'''
		self.methods.append("DFTB")
		for i in p.range(0, len(self.fileLists) ):				
			fle2 = os.path.basename(self.fileLists[i][:-4])
			_scratch = os.path.join(self.baseName, fle2)				
			if not os.path.exists(_scratch):
				os.makedirs(_scratch)
							
			self.molecule.electronicState = ElectronicState.WithOptions(charge       = self.charge 		, 
			                                                          	multiplicity = self.multiplicity )
			#-------------------------------------------------------------
			_QCmodel = QCModelDFTB.WithOptions( deleteJobFiles = False   ,
			                                	randomScratch  = True    ,
			                                 	scratch        = _scratch,
			                                 	skfPath        = skfPath ,
			                                 	useSCC         = True    )
			#-------------------------------------------------------------
			NBmodel = NBModelDFTB.WithDefaults()
			self.molecule.DefineQCModel( _QCmodel, qcSelection=Selection(self.pureQCAtoms) )
			self.cSystem.DefineNBModel( self.NBmodel ) # reseting the non-bonded model
			#--------------------------------------------------------------------
			self.molecule.qcModel.maximumSCCIterations=1200
			energy = self.cSystem.Energy()			
			self.molecule.coordinates3 = ImportCoordinates3( self.fileLists[i] )
				
			if self.ylen > 0:
				self.energiesArray[ lsFrames[0], lsFrames[1] ] = self.molecule.Energy()					
				self.indexArrayX[ lsFrames[0], lsFrames[1] ]   = lsFrames[0]
				self.indexArrayY[ lsFrames[0], lsFrames[1] ]   = lsFrames[1]
				tmpText = "{}".format(self.energiesArray[ lsFrames[0], lsFrames[1] ])
				tmpLog.write(tmpText)
				tmpLog.close()
			else:					
				self.energiesArray[ lsFrames[0] ] = self.molecule.Energy()
				self.indexArrayX[ lsFrames[0] ]   = lsFrames[0]
				tmpText = "{}".format(self.energiesArray[ lsFrames[1], lsFrames[0] ])
				tmpLog.write(tmpText)
				tmpLog.close()
		 
	#====================================================
	def SetRestart4Orca(self):
		'''
		Set the files to be run in the energy refinement for ORCA with the restart option.
		The function will read a files named frame*_**.eTmp written in the folder with the energy of the frame.
		If the Orca Refinement run terminate succesfully, these files will be removed and the entire log file will be written.
		'''
		_path = os.path.join(self.baseName,"") 
		tmpList = glob.glob(_path+"*.eTmp")
		for fle in tmpList:
			lf = GetFrameIndex(fle[:-5])
			File = open(fle,'r')
			energy = File.read()
			if self.ylen > 1:
				self.indexArrayX[lf[0],lf[1]] 	= lf[0]
				self.indexArrayY[lf[0],lf[1]] 	= lf[1]
				try:
					self.energiesArray[lf[0],lf[1]]	= float(energy)
				except:
					print(fle + " without energy written!")
					os.remove(fle)
			else:
				try:
					self.indexArrayX[lf[0]] 	= lf[0]
					self.energiesArray[lf[0]]	= float(energy)
				except:
					print(fle + " without energy written!")
					os.remove(fle)

		#-------------------------
		#remove files from list that already were calculated
		for fle in reversed(self.fileLists):			
			fle2 = os.path.basename(fle[:-4])
			_scratch = os.path.join(self.baseName, fle2, ".eTmp")
			if os.path.exists(_scratch):
				self.fileLists.remove(fle)			
		
	#====================================================
	def RunORCA(self,_method,_base,_NmaxThreads,_restart="no"):
		'''
		Perform energy refinement using the interface available on the pDynamo with the ORCA software, enabling QC(QM)/MM potential.
		Parameters:
		'''
		self.Print_Debug()
		self.methods.append(_method+_base)
		self.restart = _restart	
		if self.restart == "yes":
			self.SetRestart4Orca()	
		self.SMOenergies = {}	
		#self.Print_Debug()
		input()		
		#---------------------------------------------------------
		#Initiate parallel run
		with pymp.Parallel(_NmaxThreads) as p:
			#----------------------------------------
			#Initiate Loop			
			for i in p.range(0, len(self.fileLists) ):				
				fle2 = os.path.basename(self.fileLists[i][:-4])
				_scratch = os.path.join(self.baseName, fle2 )				
				_scratch2 = os.path.join(self.baseName, fle2,".eTmp" )				
				if not os.path.exists(_scratch):
					os.makedirs(_scratch)
				#----------------------------------------------
				lsFrames= GetFrameIndex(self.fileLists[i][:-4])
				#----------------------------------------------
				tmpLog = ""
				if self.ylen > 1: tmpPath = os.path.join( self.baseName,"frame{}_{}.eTmp".format(lsFrames[0],lsFrames[1]) )
				else: 			  tmpPath = os.path.join( self.baseName,"frame{}.eTmp".format(lsFrames[0]) )
				tmpLog  = open(tmpPath,'w')
				tmpText = ""
				#---------------------------------------------
				opt = ""
				options =  "\n% output\n"
				options +=  "print [ p_mos ] 1\n"
				options +=  "print [ p_overlap ] 5\n"
				options +=  "end # output\n"
				options +=  "!PrintBasis\n"
				
				#...............................................................................................
				self.molecule.electronicState = ElectronicState.WithOptions(charge       = self.charge 		, 
				                                                          	multiplicity = self.multiplicity )
				#...............................................................................................
				QCmodel = QCModelORCA.WithOptions( keywords        = [_method, _base, options], 
				                                   deleteJobFiles  = False                    ,
				                                   scratch         =_scratch                  )
				#...............................................................................................
				NBmodel = NBModelORCA.WithDefaults()
				self.molecule.DefineQCModel( QCmodel , qcSelection=Selection(self.pureQCAtoms) )
				self.molecule.DefineNBModel( NBmodel)
				self.molecule.coordinates3 = ImportCoordinates3( self.fileLists[i] )
				#---------------------------------------------------------------------------
				if self.ylen > 1:
					self.molecule.Energy()
					self.energiesArray[ lsFrames[0], lsFrames[1] ] = self.molecule.scratch.energyTerms["ORCA QC"]			
					self.indexArrayX[ lsFrames[0], lsFrames[1] ]   = lsFrames[0]
					self.indexArrayY[ lsFrames[0], lsFrames[1] ]   = lsFrames[1]
					tmpText = "{}".format(self.energiesArray[ lsFrames[0], lsFrames[1] ])
					tmpLog.write(tmpText)
					tmpLog.close()
				else:
					self.molecule.Energy()			
					self.energiesArray[ lsFrames[0] ] = self.molecule.scratch.energyTerms["ORCA QC"]	
					self.indexArrayX[ lsFrames[0] ]   = lsFrames[0]
					tmpText = "{}".format(self.energiesArray[ lsFrames[0] ])
					tmpLog.write(tmpText)
					tmpLog.close()
		#--------------------
		self.SMOenergies[self.methods[0]] = self.energiesArray
		#--------------------
		self.TreatOrcaFiles()
	#====================================================
	def RunPySCF(self,_method,_base,_SCF_type):
		'''
		'''
		self.methods.append(_method+_base)
		self.SMOenergies = {}		
		pySCF_pars = {"functional":_method,
					  "pySCF_method":_SCF_type,
					  "active_system":self.molecule,
					  "region":self.pureQCAtoms,
					  "QCcharge":self.charge,
					  "method_class":"pySCF",
					  "multiplicity":1,
					  "basis":_base}
		#---------------------------------------------------------
		#Initiate parallel run
		#----------------------------------------
		#Initiate Loop

		for i in range(0, len(self.fileLists) ):
			lsFrames= GetFrameIndex(self.fileLists[i][:-4])
			pySCF_pars["molden_name"] = os.path.join( self.baseName, os.path.basename(self.fileLists[i])[:-4] + ".molden") 
			qcmol = QuantumMethods(pySCF_pars)
			qcmol.system.coordinates3 = ImportCoordinates3(self.fileLists[i])
			qcmol.Set_QC_System()
			if self.ylen > 1: 
				self.energiesArray[ lsFrames[0], lsFrames[1] ] = qcmol.system.Energy(log=None)
				self.indexArrayX[lsFrames[0] , lsFrames[1] ] = lsFrames[0]
				self.indexArrayY[lsFrames[0] , lsFrames[1] ] = lsFrames[1]
			else:
				self.energiesArray[ lsFrames[0] ] = qcmol.system.Energy(log=None)
				self.indexArrayX[lsFrames[0]] = lsFrames[0]  
		
		self.SMOenergies[self.methods[0]] = self.energiesArray

	#====================================================
	def TreatOrcaFiles(self):
		'''
		Rename orca files on the scratch folder, bringing them to the base folder with the name related with the respective frames
		'''
		outFiles = glob.glob( self.baseName+"/frame*/"+"orcaJob.log" )
		for out in outFiles:
			outS = out.split("/")
			finalPath = os.path.join( self.baseName, outS[-2] + ".out" )
			shutil.move(out,finalPath)

	#====================================================
	def WriteLog(self):
		'''
		Write calculate energies to file.
		'''
		if self.ylen > 0:
			self.text += "x y Enrgy method\n"
			for smo in self.methods:
				for i in range(self.xlen):
					for j in range(self.ylen):
						self.text +="{} {} {} {}\n".format(self.indexArrayX[ i, j ],self.indexArrayY[ i,j ], self.SMOenergies[smo][i,j] - self.SMOenergies[smo][0,0], smo)
		else:
			self.text += "x y Enrgy method\n"
			for smo in self.methods:
				for i in range(self.xlen):
					self.text +="{} {} {}\n".format(self.indexArrayX[i], self.SMOenergies[smo][i] - self.SMOenergies[smo][0], smo)
		#--------------------------------------------------------------
		_filename = os.path.join(self.baseName,"energy.log")
		#----------------------------
		logFile = open(_filename,'w')
		logFile.write(self.text)
		logFile.close()
		return(_filename)
		#----------------------------
		#filesTmp = glob.glob( self.baseName+"/*.eTmp" )
		#for ftpm in filesTmp: os.remove(ftpm)	

	#================================================================
	def Print_Debug(self):
		'''
		'''
		print("="*40)
		print("Printing debug information of EnergyRefinement Class")
		print("Folder with trajectory: {}".format(self.trajFolder) )
		print("Xlen: {}".format(self.xlen))
		print("Ylen: {}".format(self.ylen))
		print("File lists len: {}".format(len(self.fileLists)) )
		print("="*40)
		


#================================================================================================#
#======================================END OF THE FILE===========================================#
#================================================================================================#

