 #!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = Simulation.py

#--------------------------------------------------------------
import os, glob, sys
import numpy as np
#-----------------------------------------------------------------------------------------------------
#Loading own libraries
#-------------------------------------------------------------
from .EnergyAnalysis         import EnergyAnalysis
from .TrajectoryAnalysis 	import TrajectoryAnalysis
from .Analysis import Analysis
#-------------------------------------------------------------
from .GeometrySearcher 	    import GeometrySearcher
from .RelaxedScan 			import SCAN
from .MolecularDynamics  	import MD
from .UmbrellaSampling  	    import US
from .PotentialOfMeanForce   import PMF
from .ReactionCoordinate 	import ReactionCoordinate
from .EnergyRefinement	 	import EnergyRefinement
from .ScanRefinement			import ScanRefinement
#--------------------------------------------------------------
#loading pDynamo Libraries
from pBabel                    import *                                     
from pCore                     import *
#---------------------------------------                                     
from pMolecule                 import *                              
from pMolecule.MMModel         import *
from pMolecule.NBModel         import *                                     
from pMolecule.QCModel         import *
#---------------------------------------
from pScientific               import *                                     
from pScientific.Arrays        import *                                     
from pScientific.Geometry3     import *                                     
from pScientific.RandomNumbers import *                                     
from pScientific.Statistics    import *
from pScientific.Symmetry      import *
#--------------------------------------                                    
from pSimulation               import *
#=============================================================
class Simulation:
	'''
	Class to set up preset simulations to be perfomed
	'''
	def __init__(self,_parameters):
		'''
		Deafault constructor
		'''
		self.molecule   = _parameters["active_system"]
		self.parameters = None
		self.Initiate_Parameters(_parameters)
		self.baseFolder = _parameters["project_folder"]

		self.restart    = False
		self.adaptative = False
		if self.parameters["restart"] 	 == "yes": self.restart    = True
		if self.parameters["adaptative"] == "yes": self.adaptative = True

	#======================================================================
	def Initiate_Parameters(self,_parameters):
		'''		
		Initiate parameters attrubute dict and the deafult values of each entry.
		'''	
		
		self.parameters = {
			"restart":"not",
			"analysis_only":"not",
			"NmaxThreads":1,
			"temperature":300.15,
			"log_frequency":0,
			"sampling_factor":0,
			"save_format":".dcd",
			"save_frequency":0,
			"adaptative":"not",
			"trajectory_name":"trajectory.ptGeo",
			"seed":3029202049,
			"QCcharge":0,
			"charge":0,
			"multiplicity":1,
			"correct_QMMM_charge":"no",
			"pySCF_method":"RHF",
			#parameters energy refinement
			"xnbins":0,
			"ynbins":0,
			"methods_lists":["am1","rm1"],
			"mopac_keywords":["AUX","LARGE"],
			"orca_method":"HF",
			"functional":"HF",
			"Software":"internal",
			#parameters geometry opt
			"maxIterations":500,
			"save_pdb":False,
			"rmsGradient":0.1,
			"optmizer":"ConjugatedGradient",
			#scan parameters
			"rc_1":None,
			"rc_2":None,
			"ndim":1,
			"dincre_rc1":0.1,
			"dincre_rc2":0.1,
			"nsteps_rc1":0,
			"nsteps_rc2":0,
			"xlim":None,
			"ylim":None,
			#molecular dynamics parameters
			"MD_method":"LeapFrog",
			"pressure":1.0,
			"pressure_coupling":2000,
			"pressure_control":"False",
			"temperature_scale_option":"constant",
			"temperature_scale":10.0,
			"timeStep":0.001,
			"coll_freq":25.0,
			"start_temperature":20.0,
			"equilibration_nsteps":0,
			"heating_nsteps":0,
			"production_nsteps":0,
			"sampling_heating":0,
			"sampling_equilibration":0,
			"sampling_production":0,
			#restriction parameters
			"force_constants":[],
			#free energy parameters
			"relax":"False",
			"crd_format":"pkl",
			"optimize_US":False,
			"analysis_only":"False",
			#thermo parameters
			"cycles":10,
			"mode":0,
			"frames":10,
			#path finder parameters
			"NEB_bins":0,
			"init_coord":"False",
			"final_coord":None,
			"spring_force_constant":None,
			"fixed_terminal_images":None,
			"RMS_growing_intial_string":None,
			"reverse_rc1":"no",
			"reverse_rc2":"no",
			"contour_lines":15,
			"fig_size":[7,5],
			"crd_labels":[],
		}

		for key in _parameters.keys(): self.parameters[key] = _parameters[key]

	#=======================================================================
	def Execute(self):
		'''
		Function to call the class method to execute the preset simulation
		Mandatory keys:
			"simulation_type": Name of the simulation to execute
		'''		
		#-------------------------------------------------------------------------------------------------------------------------
		if 	 self.parameters["simulation_type"] == "Energy_Refinement": 			self.EnergyRefine()		
		elif self.parameters["simulation_type"] == "Geometry_Optimization":			self.GeometryOptimization()
		elif self.parameters["simulation_type"] == "Relaxed_Surface_Scan":	 		self.RelaxedSurfaceScan()
		elif self.parameters["simulation_type"] == "ScanRefinement":				self.ScanRefinement()
		elif self.parameters["simulation_type"] == "Molecular_Dynamics":			self.MolecularDynamics()	
		elif self.parameters["simulation_type"] == "Restricted_Molecular_Dynamics": self.RestrictedMolecularDynamics()
		elif self.parameters["simulation_type"] == "Umbrella_Sampling":				self.UmbrellaSampling()
		elif self.parameters["simulation_type"] == "Normal_Modes":					self.NormalModes()		
		elif self.parameters["simulation_type"] == "Delta_Free_Energy":				self.DeltaFreeEnergy()		
		elif self.parameters["simulation_type"] == "NEB":							self.ReactionSearchers()		
		elif self.parameters["simulation_type"] == "SAW":							self.ReactionSearchers()		
		elif self.parameters["simulation_type"] == "Baker_Saddle":					self.ReactionSearchers()
		elif self.parameters["simulation_type"] == "Steep_Path_Searcher":			self.ReactionSearchers()				
		elif self.parameters["simulation_type"] == "Simulating_Annealing":			self.SimulatingAnnealing()		
		elif self.parameters["simulation_type"] == "Steered_Molecular_Dynamics":	self.SMD()		
		elif self.parameters["simulation_type"] == "Monte_Carlo":					self.MonteCarlo()
		return(self.molecule.system)				
		
	#=================================================================================================================
	def EnergyRefine(self):
		'''
		Set up and execute energy refinement using a series of methods	
		'''
		dimensions    = [ self.parameters["xnbins"], 0 ] 
		nmaxthreads   =  self.parameters["NmaxThreads"]
		_trajfolder   = "single"
		_type         = "1DRef"
		if self.parameters["ynbins"] > 0: _type = "2DRef"
		if "ynbins"        in self.parameters: dimensions[1] = self.parameters["ynbins"]
		if "source_folder" in self.parameters: _trajfolder   = self.parameters["source_folder"] 
		#------------------------------------------------------------------
		ER = EnergyRefinement(self.molecule.system          ,
							  _trajfolder  		            ,
							  self.parameters["folder"]     ,
							  dimensions                    ,
							  self.parameters["QCcharge"]   ,
							  self.parameters["multiplicity"])
		
		#------------------------------------------------------------------
		if 	 self.parameters["Software"] == "pDynamo"   : ER.RunInternalSMO(self.parameters["methods_lists"],nmaxthreads)
		elif self.parameters["Software"] == "pDynamoDFT": ER.RunInternalDFT(self.parameters["functional"],self.parameters["basis"],nmaxthreads)
		elif self.parameters["Software"] == "DFTBplus"  : ER.RunDFTB()
		elif self.parameters["Software"] == "pySCF"     : ER.RunPySCF(self.parameters["functional"],self.parameters["basis"],_SCF_type=self.parameters["pySCF_method"])
		elif self.parameters["Software"] == "ORCA"		: ER.RunORCA(self.parameters["orca_method"],self.parameters["basis"],nmaxthreads,_restart=self.parameters["restart"])
		elif self.parameters["Software"] == "mopac" or self.parameters["Software"]=="MOPAC":
			_mopacKeyWords = ["AUX","LARGE"] 
			if "mopac_keywords" in self.parameters:
				for key in self.parameters["mopac_keywords"]: _mopacKeyWords.append(key)
			ER.RunMopacSMO(self.parameters["methods_lists"],_mopacKeyWords)
		#------------------------------------------------------------
		log_path = ER.WriteLog()		
		EA       = EnergyAnalysis(self.parameters["xnbins"],self.parameters["ynbins"],_type=_type)		
		EA.ReadLog(log_path)
		crd2_label = None
		#--------------------------------------------------------
		if len(self.parameters["crd_labels"]) == 1:
			crd1_label = self.parameters["crd_labels"][0]
		elif len(self.parameters["crd_labels"]) == 2:
			crd1_label = self.parameters["crd_labels"][0]
			crd2_label = self.parameters["crd_labels"][1]
		elif _type == "1DRef":
			try: crd1_label = self.molecule.reactionCoordinates[0].label
			except: crd1_label = "Reaction Path frames (n)"
		elif _type == "2DRef":
			try: crd1_label = self.molecule.reactionCoordinates[0].label
			except: crd1_label = "Reaction Path frames (n)"
			try: crd2_label = self.molecule.reactionCoordinates[1].label
			except: crd2_label = "Reaction Path frames (n)"
		
		_reverse_rc1 = False
		_reverse_rc2 = False
		if self.parameters["reverse_rc1"] == "yes": _reverse_rc1 = True 
		if self.parameters["reverse_rc2"] == "yes": _reverse_rc2 = True 
		if   _type == "1DRef": EA.MultPlot1D(crd1_label)
		elif _type == "2DRef": EA.MultPlot2D(self._parameters["contour_lines"],crd1_label,crd2_label,self.parameters["fig_size"],_reverse_rc1,_reverse_rc2)	
	#==================================================================
	def GeometryOptimization(self):
		'''
		Set up and execture the search of local minima for the system passed						
		'''
		_traj_name = None
		if "optmizer" 		 in self.parameters: _Optimizer = self.parameters["optmizer"]
		if "trajectory_name" in self.parameters: _traj_name = self.parameters["trajectory_name"]
		Gopt = GeometrySearcher(self.molecule.system,self.baseFolder,_trajName=_traj_name)		
		Gopt.ChangeDefaultParameters(self.parameters)
		Gopt.Minimization(self.parameters["optmizer"])
		Gopt.Finalize()
		self.molecule.system = Gopt.molecule

	#==================================================================
	def RelaxedSurfaceScan(self, plot = True):
		'''
		Set up and execute one/two-dimensional relaxed surface scans 
		By the defualt the PKLs were saved on a child folder from the base path passed in the parameters, named "ScanTraj.ptGeo"
		The trajectory can be saved as files of the formats allowed by pDynamo 3.0		
		'''
		X = self.parameters["nsteps_rc1"]
		Y = self.parameters["nsteps_rc2"]
		if X > 0: 
			_type = "1D"		
			scan = SCAN(self.molecule.system,self.baseFolder,self.parameters)
			crd2_label = None
			#--------------------------------------------------------------------
			scan = SCAN(self.molecule.system,self.baseFolder,self.parameters["optmizer"],self.parameters["adaptative"],self.parameters["restart"])
			scan.ChangeDefaultParameters(self.parameters)	
			#--------------------------------------------------------------------
			self.molecule.reactionCoordinates[0].SetInformation(self.molecule.system,self.parameters["dincre_rc1"])
			scan.SetReactionCoord(self.molecule.reactionCoordinates[0])
			if self.molecule.rcs == 2:
				print("Two dimensions relaxed surface scan!")
				self.molecule.reactionCoordinates[1].SetInformation(self.molecule.system,self.parameters["dincre_rc2"])
				scan.SetReactionCoord(self.molecule.reactionCoordinates[1])
				scan.Run2DScan(X, Y)
				_type = "2D"
				crd2_label = self.molecule.reactionCoordinates[1].label
			else: 
				print("One dimension relaxed surface scan!")
				scan.Run1DScan(self.parameters["nsteps_rc1"])
			log_path = scan.Finalize()

			EA = EnergyAnalysis( X, Y, _type=_type)		
			EA.ReadLog(log_path)
			#--------------------------------------------------------
			crd1_label = self.molecule.reactionCoordinates[0].label
			
			if   _type == "1D": EA.Plot1D(crd1_label)
			elif _type == "2D":					
				self.parameters["xlim"] = [scan.reactionCoordinate1[0] , scan.reactionCoordinate1[-1] ]
				self.parameters["ylim"] = [scan.reactionCoordinate2[0] , scan.reactionCoordinate2[-1] ]	
				EA.Plot2D(self.parameters["contour_lines"],crd1_label,crd2_label,_xlim=self.parameters["xlim"],_ylim=self.parameters["ylim"],_figS=self.parameters["fig_size"])		
	#==================================================================================	
	def ScanRefinement(self):
		'''
		'''
		self.parameters["folder"] = self.baseFolder
		scan = ScanRefinement(self.parameters)
		scan.SetReactionCoord(self.molecule.reactionCoordinates[0])
		scan.SetReactionCoord(self.molecule.reactionCoordinates[1])
		scan.RunRelaxedRefinement(self.parameters["functional"], self.parameters["basis"], self.parameters["pySCF_method"])

		log_path = scan.WriteLog()		
		EA       = EnergyAnalysais(scan.xsize,0,_type="1DRef")		
		EA.ReadLog(log_path)
		#--------------------------------------------------------
		crd1_label = "Reaction Path frames (n)"		
		EA.Plot1D(crd1_label)
	#==================================================================================	
	def MolecularDynamics(self):
		'''
		Set up and execute molecular dynamics simulations		
		'''				
		MDrun = MD(self.molecule.system,self.baseFolder,self.parameters)		

		if self.parameters["heating_nsteps"] 	   > 0: 
			MDrun.HeatingSystem(self.parameters["heating_nsteps"],self.parameters["sampling_heating"])
			if self.parameters["sampling_heating"] > 0:
				_trajAN = TrajectoryAnalysis(MDrun.trajectoryNameCurr,MDrun.molecule,MDrun.timeStep*MDrun.nsteps)
				_trajAN.CalculateRG_RMSD()
				_trajAN.PlotRG_RMS()
		if self.parameters["equilibration_nsteps"] > 0: 
			MDrun.RunProduction(self.parameters["equilibration_nsteps"],self.parameters["sampling_equilibration"],_equi=True)
			if self.parameters["sampling_equilibration"] > 0:
				_trajAN = TrajectoryAnalysis(MDrun.trajectoryNameCurr,MDrun.molecule,MDrun.timeStep*MDrun.nsteps)
				_trajAN.CalculateRG_RMSD()
				_trajAN.PlotRG_RMS()
		if self.parameters["production_nsteps"]    > 0: 
			MDrun.RunProduction(self.parameters["production_nsteps"],self.parameters["sampling_production"])
			if self.parameters["sampling_production"] > 0:
				_trajAN = TrajectoryAnalysis(MDrun.trajectoryNameCurr,MDrun.molecule,MDrun.timeStep*MDrun.nsteps)
				_trajAN.CalculateRG_RMSD()
				_trajAN.PlotRG_RMS()		
		#------------------------------------------------------------
		MDrun.Finalize()			
	#==================================================================
	def RestrictedMolecularDynamics(self):
		'''
		Set up and execute molecular dynamics simulations.
		'''
		#----------------------------------------------------------------
		restraints = RestraintModel()
		self.molecule.system.DefineRestraintModel( restraints )		
		rcs = []
		rc1 = self.molecule.reactionCoordinates[0]
		rc1.SetInformation(self.molecule.system,self.parameters["dincre_rc1"])
		rc1.GetRCLabel(self.molecule.system)
		rc1.SetInformation(self.molecule.system,0.0,)
		rcs.append(rc1)
		#-------------------------------------------------------------------
		restrainDimensions = self.parameters['ndim']
		forcK_1 = self.parameters["force_constants"][0]
		#-------------------------------------------------------------------
		
		nDims = self.parameters['ndim']
		rc2 = None
		if nDims == 2:
			rc2 = self.molecule.reactionCoordinates[1]
			rc2.SetInformation(self.molecule.system,self.parameters["dincre_rc2"])
			rc2.GetRCLabel(self.molecule.system)
			rc2.SetInformation(self.molecule.system,0.0)
			forcK_2 = self.parameters["force_constants"][1]
		#-------------------------------------------------------------------
		distance = rc1.minimumD
		rmodel = RestraintEnergyModel.Harmonic( distance, forcK_1 )
		if rc1.nAtoms == 3:				
			restraint = RestraintMultipleDistance.WithOptions( energyModel=rmodel, distances=[ [ rc1.atoms[1], rc1.atoms[0], rc1.weight13 ], [ rc1.atoms[1], rc1.atoms[2], rc1.weight31 ] ] ) 
		elif rc1.nAtoms == 2:				
			restraint = RestraintDistance.WithOptions( energyModel=rmodel, point1=rc1.atoms[0], point2=rc1.atoms[1] )
		elif rc1.nAtoms == 4:
			rmodel = RestraintEnergyModel.Harmonic( distance, forcK_1, period = 360.0 )
			restraint = RestraintDihedral.WithOptions( energyModel=rmodel, point1=rc1.atoms[0],point2=rc1.atoms[1],point3=rc1.atoms[2],point4=rc1.atoms[3] )
		restraints['M1'] =  restraint
		#-------------------------------------------------------------------
		if nDims == 2:
			distance = rc2.minimumD
			rmodel = RestraintEnergyModel.Harmonic( distance, forcK_2 )
			if rc2.nAtoms == 3:				
				restraint = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ rc2.atoms[1], rc2.atoms[0], rc2.weight13 ], [ rc2.atoms[1], rc2.atoms[2], rc2.weight31 ] ] ) 
			elif rc2.nAtoms == 2:				
				restraint = RestraintDistance.WithOptions( energyModel=rmodel, point1=rc2.atoms[0], point2=rc2.atoms[1] )
			elif rc2.nAtoms == 4:
				rmodel = RestraintEnergyModel.Harmonic( distance, forcK_2, period = 360.0 )
				restraint = RestraintDihedral.WithOptions( energyModel=rmodel, point1=rc2.atoms[0],point2=rc2.atoms[1],point3=rc2.atoms[2],point4=rc2.atoms[3] )	
			restraints['M2'] =  restraint	
			rcs.append(rc2)	
		#----------------------------------------------------------------
		traj_name = "trajectory"
		if "trajectory_name" in self.parameters: traj_name = self.parameters["trajectory_name"]
		
		MDrun = MD(self.molecule.system,self.baseFolder,self.parameters)
		if self.parameters["heating_nsteps"] 	   > 0: 
			MDrun.HeatingSystem(self.parameters["heating_nsteps"],self.parameters["sampling_heating"])
			if self.parameters["sampling_heating"] > 0:
				_trajAN = TrajectoryAnalysis(MDrun.trajectoryNameCurr,MDrun.molecule,MDrun.timeStep*MDrun.nsteps)
				_trajAN.CalculateRG_RMSD()
				_trajAN.PlotRG_RMS()
		if self.parameters["equilibration_nsteps"] > 0: 
			MDrun.RunProduction(self.parameters["equilibration_nsteps"],self.parameters["sampling_equilibration"],_equi=True)
			if self.parameters["sampling_equilibration"] > 0:
				_trajAN = TrajectoryAnalysis(MDrun.trajectoryNameCurr,MDrun.molecule,MDrun.timeStep*MDrun.nsteps)
				_trajAN.CalculateRG_RMSD()
				_trajAN.PlotRG_RMS()
		if self.parameters["production_nsteps"]    > 0: 
			MDrun.RunProduction(self.parameters["production_nsteps"],self.parameters["sampling_production"])
			if self.parameters["sampling_production"] > 0:
				_trajAN = TrajectoryAnalysis(MDrun.trajectoryNameCurr,MDrun.molecule,MDrun.timeStep*MDrun.nsteps)
				_trajAN.CalculateRG_RMSD()
				_trajAN.PlotRG_RMS()	
				_trajAN.DistancePlots(rcs)
				_trajAN.ExtractFrames()

	#=======================================================================
	def UmbrellaSampling(self):
		'''
		Set up and execute umbrella sampling simulations and Free energy calculations for reaction path trajectory.		
		'''
		#------------------------------------------------------------------
		rc1 = self.molecule.reactionCoordinates[0]
		rc1.GetRCLabel(self.molecule.system)
		rc1.SetInformation(self.molecule.system,0.0)		
		sampling   = self.parameters["sampling_production"]
		_crdFormat = self.parameters["crd_format"] 
		

		nDims = self.parameters['ndim']
		rc2   = None
		if nDims == 2:
			rc2 = self.molecule.reactionCoordinates[1]
			rc2.GetRCLabel(self.molecule.system)
			rc2.SetInformation(self.molecule.system,0.0)
		#---------------------------------------
		USrun = US(self.molecule.system  				  ,
			       self.baseFolder 						  ,
			       self.parameters["equilibration_nsteps"],
			       self.parameters["production_nsteps"]   ,
			       self.parameters["MD_method"]           ,
			       RESTART=self.restart                   ,
			       ADAPTATIVE=self.adaptative             ,
			       OPTIMIZE=self.parameters["optimize_US"]                     )
		#---------------------------------------
		USrun.ChangeDefaultParameters(self.parameters)
		USrun.SetMode(rc1)

		if self.parameters["analysis_only"] == "yes":
			self.parameters["active_system"] = self.molecule.system
			self.parameters["folder"] = self.baseFolder
			self.parameters["source_folder"] = self.baseFolder
			self.parameters["analysis_type"] = "PMF"
			WHAM = Analysis(self.parameters)
			WHAM.PMFAnalysis()
			USrun.Finalize()	
		else:		
			if self.parameters["ndim"]   == 1: 
				USrun.Run1DSampling(self.parameters["source_folder"],_crdFormat,sampling)
			elif self.parameters["ndim"] == 2:
				USrun.SetMode(rc2)
				USrun.Run2DSampling(self.parameters["source_folder"],_crdFormat,sampling)
			USrun.Finalize()		
			self.parameters["active_system"] = self.molecule.system
			self.parameters["folder"] = self.baseFolder
			self.parameters["source_folder"] = self.baseFolder
			self.parameters["analysis_type"] = "PMF"
			WHAM = Analyasis(self.parameters)
			WHAM.PMFAnalysis()		
	
	#==========================================================================
	def NormalModes(self):
		'''
		Simulation preset to calculate the normal modes and to write thr trajectory for a specific mode.			
		'''
		mode 		= 0
		temperature = 300.15
		Cycles 		= 10 
		Frames  	= 10 
		#-------------------------------
		if "temperature" in self.parameters: temperature = self.parameters["temperature"]
		if "cycles" 	 in self.parameters: Cycles 	 = self.parameters["cycles"]
		if "frames" 	 in self.parameters: Frames 	 = self.parameters["frames"]
		if "mode" 	 	 in self.parameters: mode   	 = self.parameters["mode"]
		#-------------------------------
		NormalModes_SystemGeometry ( self.molecule.system, modify = ModifyOption.Project )
		if mode > 0:
			trajectory = ExportTrajectory ( os.path.join (self.baseFolder, "NormalModes","ptGeo"), self.molecule.system )
			NormalModesTrajectory_SystemGeometry(	self.molecule.system		      ,
                                       			 	trajectory                ,
                                       				mode        = mode	      ,
                                       				cycles      = Cycles      ,
                                       				frames      = Frames 	  ,
                                       				temperature = temperature )
	#==========================================================================
	def DeltaFreeEnergy(self):
		'''
		Calculate the free energy difference between two configurations of the system using the 
		statistical thermodynamics partition functions from through the normal modes calculations
		Mandatory keys:
		Optional keys :
		'''		
		#initial Structure
		pressure       = 1.0
		temperature    = 300.15 
		symmetryNumber = 1

		if "pressure" in self.parameters: pressure = self.parameters["pressure"]

		self.molecule.system.coordinates3 = ImportCoordinates3(self.parameters["initial_coordinates"])
		e0 = self.molecule.system.Energy()
		NormalModes_SystemGeometry( self.molecule.system, modify = ModifyOption.Project )
		Gibbs = [] 
		tdics = ThermodynamicsRRHO_SystemGeometry ( self.molecule.system 							,
                                                    pressure       = pressure       	,
                                                    symmetryNumber = self.symmetryNumber 	,
                                                    temperature    = self.temperature    	)
		Gibbs.append( tdics["Gibbs Free Energy"] )
    	# Final struct
		self.molecule.system.coordinates3 = ImportCoordinates3(self.parameters["final_coordinates"])
		e1 = self.molecule.system.Energy()
		NormalModes_SystemGeometry ( self.molecule.system, modify = ModifyOption.Project )
    	 
		tdics = ThermodynamicsRRHO_SystemGeometry ( self.molecule.system 							,
                                                    pressure       = self.pressure       	,
                                                    symmetryNumber = self.symmetryNumber 	,
                                                    temperature    = self.temperature    	)
		Gibbs.append( tdics["Gibbs Free Energy"] )
	#=========================================================================	
	def ReactionSearchers(self):
		'''
		Class method to set up and execute Nudget Elastic Band simulations to generate a reaction path trajectory
		Mandatory keys in self.parameters:
			"NEB_bins"  : Number of points in the NEB trajectory
			"init_coord":
			"final_coord":

		Optional keys in self.parameters:
			"spring_force_constant"    :
			"fixed_terminal_images"    :
			"RMS_growing_intial_string":
			"refine_methods"           : 
			"crd1_label"               :
			"show"                     :
			"xlim_list"                :
		'''


		_traj_name = "ReactionPath"
		if "trajectory_name" in self.parameters: _traj_name = self.parameters["trajectory_name"]
		RSrun = GeometrySearcher(self.molecule.system,self.baseFolder,_trajName=_traj_name)		
		RSrun.ChangeDefaultParameters(self.parameters)		

		if   self.parameters["simulation_type"] == "NEB"                : RSrun.NudgedElasticBand(self.parameters)
		elif self.parameters["simulation_type"] == "SAW"                : RSrun.SelfAvoidWalking(self.parameters)
		elif self.parameters["simulation_type"] == "SteepDescent_path"  : RSrun.SteepestDescentPathSearch(self.parameters)
		elif self.parameters["simulation_type"] == "Baker_Saddle"       : RSrun.BakerSaddleOptimizer(self.parameters) 

		
		nmaxthreads = 1
		if "NmaxThreads" in self.parameters: nmaxthreads = self.parameters["NmaxThreads"]

		refMethod = []
		if "methods_lists" in self.parameters: refMethod = self.parameters["methods_lists"]
		if len(refMethod) > 0: 
			ER = EnergyRefinement(self.molecule.system  					        ,
								  RSrun.trajectoryName                      ,
								  self.baseFolder                           ,
								  [self.parameters["traj_bins"],0]          ,
								  self.molecule.system.electronicState.charge      ,
								  self.molecule.system.electronicState.multiplicity)
			ER.RunInternalSMO(refMethod,nmaxthreads)
			ER.WriteLog()
			crd1_label 	= "Reaction Coordinate #1"
			xlim 		= [ 0, self.parameters["traj_bins"] ]
			show  		= False
			#check parameters for plot
			if "crd1_label" in self.parameters: crd1_label = self.parameters["crd1_label"]
			if "xlim_list"  in self.parameters: xlim       = self.parameters["xlim_list"]
			if "show" 		in self.parameters: show       = self.parameters["show"]
			#------------------------------------------------------------				
			EA = EnergyAnalysis(self.parameters["traj_bins"],1,_type="1DRef")
			EA.ReadLog( os.path.join(ER.baseName,"energy.log") )
			EA.MultPlot1D(crd1_label,show)	
			RSrun.Finalize()
	#=========================================================================
	def MonteCarlo(self):
		pass
	
	#=========================================================================
	def SimulatingAnnealing(self):
		'''
		Set up and execute Simulate annealing simulations	
		Mandatory keys in self.parameters:
		Optional keys in self.parameters:	
		'''
		pass
	#=========================================================================
	def SMD(self):
		'''
		Set up and execute Steered Molecular Dynamics simulations
		Mandatory keys in self.parameters:
		Optional keys in self.parameters:
		'''
		pass
	#=========================================================================	
	def Print(self):
		'''
		Printing information of the simulations that will be run.
		'''
		print("Simulation Type: {}".format(self.parameters["simulation_type"]) )
		print("Working folder: {}".format(self.parameters["folder"]) )

#=============================================================================
#========================END OF THE FILE======================================
#=============================================================================
