 #!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = Analysis.py

#--------------------------------------------------------------
import os, glob, sys
import numpy as np
#--------------------
from pMolecule import *
from pBabel                    import *                                     
from pCore                     import *                                     
from pSimulation               import PruneByAtom
                           
from pMolecule.QCModel         import *
from pMolecule.NBModel import NBModelORCA


from addOns.pySCF import NBModelPySCF , \
                             QCModelPySCF
import MopacQCMMinput
import pyscf
#==============================================================
class QuantumMethods:
	'''
	Classe to set up quantum chemical method to the system.
	'''	
	def __init__(self,_parameters):
		'''
		Default constructor
		'''
		self.methodClass 	= "SMO" #  SMO, HF, DFT, ORCA, DFTB, MOPAC and PYSCF
		self.Hybrid  		= False  
		self.selection      = [] 
		self.systemBase     = None 
		self.QuantumSystem  = None 
		self.convergerLevel = "standard" # veryLoose, loose, standard, tight, veryTight 
		self.converger      = None
		self.qcSystem       = None 
		self.system         = None
		self.qcModel        = None
		self.pars           = None

		self.Check_Parameters(_parameters)
		self.system = _parameters["active_system"]
		NBmodel     = None
		if self.pars["region"]: self.Hybrid = True
		#---------------------------------------------
		if self.Hybrid:
			atomlist = []
			for sel in self.pars["region"]:
				if type(sel) == int:
					self.selection.append(sel)
				elif type(self.pars["region"]) == list:
					for i in range( len(sel) ):
						self.selection.append( sel[i] )
			self.selection = Selection.FromIterable(self.selection)

		newSelection           = AtomSelection.ByComponent(self.system,self.selection)
		newSystem    		   = PruneByAtom(self.system, Selection(newSelection) )	
		#ExportSystem( os.path.join( "qcSystem.pdb"), newSystem )
		try:
			new_charge  = self.GetQCCharge(newSystem)			
			if not new_charge == self.pars["QCcharge"]:
				if self.pars["correct_QMMM_charge"] == "yes":
					self.pars["QCcharge"] = new_charge
		except:
			pass

		
		#---------------------------------------------       	
		self.convergerLevel = self.pars["converger"]
		self.Set_Converger()
		self.system.electronicState = ElectronicState.WithOptions( charge = self.pars["QCcharge"],
																   multiplicity = self.pars["multiplicity"] )

	#--------------------------------------------------
	def Check_Parameters(self,_parameters):
		'''
		'''
		self.pars = {
			"method_class":"SMO",
			"region":None,
			"QCcharge":0,
			"multiplicity":1,
			"functional":"HF",
			"Hamiltonian":"am1",
			"basis":"sto3g",
			"fit_basis":"dgauss-a1-dftjfit",
			"scratch":os.getcwd(),
			"skfPath":"/home/igorchem/programs/pDynamo3/examples/dftbPlus/data/skf",
			"converger":"standard",
			"center_atom":-1,
			"new_radius_qc":0.0,
			"pySCF_method":"RHF",
			"NmaxThreads":1,
			"molden_name":"file.molden"
		}
		for key in _parameters.keys(): self.pars[key] = _parameters[key]
		self.methodClass = self.pars["method_class"]

	#------------------------------------------------------
	def Set_QC_System(self):
		'''
		'''
		if   self.methodClass == "SMO":     self.Set_SMO_internal()
		elif self.methodClass == "pySCF":   self.Set_pySCF() 
		elif self.methodClass == "abinitio":self.Set_Abinitio()
		elif self.methodClass == "ORCA":    self.Set_ORCA()

	#------------------------------------------------------
	def Set_SMO_internal(self):
		'''
		'''
		NBmodel             = self.system.nbModel
		self.system.nbModel = None
		self.qcModel = QCModelMNDO.WithOptions( hamiltonian = self.pars["Hamiltonian"],
												converger=self.converger )


		if self.Hybrid: 
			self.system.DefineQCModel( self.qcModel, qcSelection=self.selection )
			self.system.DefineNBModel( NBmodel, assignQCMMModels=self.Hybrid )
		else: self.system.DefineQCModel( self.qcModel )		
		
		self.system.DefineNBModel( NBmodel, assignQCMMModels=self.Hybrid )

	#--------------------------------------------------------
	def Set_pySCF(self):
		'''
		'''
		NBmodel = NBModelPySCF.WithDefaults( )			
		qcModel = QCModelPySCF.WithOptions( deleteJobFiles = False       ,
											functional     = self.pars["functional"],
                                            method         = self.pars["pySCF_method"],
                                            mf_kwargs      = { 'diis'    : pyscf.scf.ADIIS ( ) }, 
                                            mole_kwargs    = { 'verbose' : 0 , "molden_name":self.pars["molden_name"]} ,
                                            orbitalBasis   = self.pars["basis"] )

		if self.Hybrid: 
			self.system.DefineQCModel( qcModel, qcSelection=self.selection )
			self.system.DefineNBModel( NBmodel, assignQCMMModels=self.Hybrid )
			self.Export_QC_System()
		else:	
			self.system.DefineQCModel( qcModel )
		self.qcModel = qcModel

	#--------------------------------------------------------
	def Set_Abinitio(self):
		'''
		'''
		NBmodel             = self.system.nbModel
		self.system.nbModel = None

		_gridIntegrator = DFTGridIntegrator.WithOptions(accuracy = DFTGridAccuracy.Medium,
															inCore   = True                  )

		self.qcModel = QCModelDFT.WithOptions(converger   	 = self.converger 		  ,
												  functional     = self.pars["functional"],
												  orbitalBasis	 = self.pars["basis"] 	  ,
												  gridIntegrator = _gridIntegrator        ,
												  fitBasis       = self.pars["fit_basis"] )

		
		if self.Hybrid: 
			self.system.DefineQCModel( self.qcModel, qcSelection=self.selection )
		else          : self.system.DefineQCModel( self.qcModel )		
		
		self.system.DefineNBModel( NBmodel, assignQCMMModels=self.Hybrid )


	#.................................................................................
	def Set_ORCA(self):
		'''
		'''
		options = "\n% output\n"
		options += "print [ p_mos ] 1\n"
		options += "print [ p_overlap ] 5\n"
		options += "end # output\n"
		options += "!PrintBasis\n"
		#options +="%maxcore 1000\n"
		#options +="%pal\n"
		#options +="nprocs 2\n"
		#options +="end\n"
		_keyWords = [ self.pars["functional"],
					  self.pars["basis"],
					  #"PAL{}".format(self.pars["NmaxThreads"]),
					  options ]


		NBmodel  = NBModelORCA.WithDefaults()
		self.qcModel = QCModelORCA.WithOptions( keywords = _keyWords                  , 
                                            	deleteJobFiles  = False               ,
                                            	scratch         = self.pars["scratch"])
		
		if self.Hybrid: 
			self.system.DefineQCModel( self.qcModel, qcSelection=self.selection )
			self.system.DefineNBModel( NBmodel, assignQCMMModels=True )
		else: self.system.DefineQCModel( self.qcModel )		
			         
	#----------------------------------------------------------------------------
	def Export_QC_System(self,baseName = None):
		'''
		'''		
		if baseName == None: baseName = os.getcwd()
		self.qcSystem = PruneByAtom(self.system,self.selection)
		ExportSystem( os.path.join( baseName,"qcSystem.pdb"), self.qcSystem ) 
		ExportSystem( os.path.join( baseName,"qcSystem.pkl"), self.qcSystem )
	#----------------------------------------------------------------------------
	def Set_Converger(self):
		'''
		'''
		EnergyTolerance  = 3.0e-4
		DensityTolerance = 1.0e-8
		MaxIterations    = 1000 

		if self.convergerLevel == "Very_Loose":
			EnergyTolerance  = 1.0e-4
			DensityTolerance = 1.0e-6
			MaxIterations    = 3500
		elif  self.convergerLevel == "Loose":
			EnergyTolerance  = 3.0e-4
			DensityTolerance = 1.0e-7
			MaxIterations    = 2500
		elif self.convergerLevel == "standard":
			EnergyTolerance  = 3.0e-4
			DensityTolerance = 1.0e-8
			MaxIterations    = 1000
		elif self.convergerLevel == "Tight":
			EnergyTolerance  = 1.0e-5
			DensityTolerance = 1.0e-8
			MaxIterations    = 800
		elif self.convergerLevel == "Very_Tight":
			EnergyTolerance  = 3.0e-6
			DensityTolerance = 1.0e-8
			MaxIterations    = 500

		self.converger = DIISSCFConverger.WithOptions( energyTolerance   = EnergyTolerance ,
												  	   densityTolerance  = DensityTolerance,
												  	   maximumIterations = MaxIterations   )
	#-------------------------------------------------------------------
	def GetQCCharge(self,_system):
		'''
		'''
		qc_charge=0.0
		mmCharges = _system.mmState.charges
		for i in range( len(mmCharges) ): qc_charge += mmCharges[i]			
		return( round(qc_charge) )
	#--------------------------------------------------------------------
	def Change_QC_Region(self):
		'''
		Redefine QC selection from a given atomic coordinates with a certain radius
		Parameters:
			_centerAtom:
			_radius    :
		'''
		qcModel = None
		if type(self.pars["center_atom"]) == list:  
			_centerAtom = [ self.pars["center_atom"][0],
							self.pars["center_atom"][1],
							self.pars["center_atom"][2] ]
			atom_list = [] 		
			for i in self.system.atoms.items:
				x    = self.system.coordinates3[i.index, 0] 
				y    = self.system.coordinates3[i.index, 1]
				z    = self.system.coordinates3[i.index, 2]
				xd   = (x-_centerAtom[0])**2
				yd   = (y-_centerAtom[1])**2
				zd   = (z-_centerAtom[2])**2
				dist = np.sqrt( xd+yd+zd ) 
				if dist < _radius: atom_list.append(i.index)

			sel          = Selection.FromIterable(atom_list)
			newSelection = AtomSelection.ByComponent(self.system,sel)
			newSystem    = PruneByAtom(self.system, Selection(newSelection) )		
			self.charge  = self.GetQCCharge(newSystem)
			_pureQCAtoms = list(newSelection)	
		#-------------------------------------------------------------------
		else:
			sel              	  = Selection.FromIterable([_centerAtom])
			newSelection  	 	  = AtomSelection.Within(self.molecule,sel,_radius)
			newSelection 	 	  = AtomSelection.ByComponent(self.molecule,newSelection)
			newSystem    	 	  = PruneByAtom(self.system, Selection(newSelection) )		
			self.pars["QCcharge"] = self.GetQCCharge(newSystem)
			self.selection   	  = list(newSelection)	

		if self.system.qcModel == None: qcModel = QCModelMNDO.WithOptions( hamiltonian = "am1" )
		else: qcModel = self.system.qcModel
		#-------------------------------------------------------------------------------
		self.molecule.electronicState = ElectronicState.WithOptions( charge = self.pars["QCcharge"], multiplicity = self.pars["multiplicity"] )
		self.molecule.DefineQCModel(qcModel, qcSelection=Selection(self.selection) )
		#-------------------------------------
	

#==================================================================================
#=============================END OF FILE==========================================
#==================================================================================
