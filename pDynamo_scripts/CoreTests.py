#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = Tests.py

#-------------------------------------------------------------
#Script file for the EasyHybrid 3.0 Core functionalities tests 

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#=================================================================
import os, glob, sys 

global NmaxThreads
NmaxThreads = 1 
#------------------------------------------------------
from pBabel             import *                                     
from pCore              import *                                     
from pMolecule          import *                    
from SimulationProject	import *
from ReactionCoordinate import *
from TrajectoryAnalysis import *
from MopacQCMMinput 	import *
from QuantumMethods		import *
#from WriteQMLog   		import * 
#-------------------------------------------------------------------
#path for the required files on the examples folder of EasyHynrid 3.0
easyhybrid   = os.path.join(os.getcwd())
ex_path      = os.path.join(easyhybrid, "examples")

scratch_path = os.path.join("TestsScratch")
timTop       = os.path.join(ex_path,"TIM","7tim.top")
timCrd       = os.path.join(ex_path,"TIM","7tim.crd")
balapkl      = os.path.join(ex_path,"bala","bAla.pkl")
meth         = os.path.join(ex_path,"pdb","methane.pdb")
#--------------------------------------------------------
if not os.path.exists(scratch_path):
	os.makedirs(scratch_path)

#=====================================================
def Scan1D_Dihedral(_nsteps,_dincre=10.0,name="Default"):
	#---------------------------------------------	
	_scanFolder = "SCAN1D_dihedral"
	if not name == "Default": _scanFolder = name
	proj=SimulationProject.From_PKL( os.path.join(balapkl)					, 
									 os.path.join(scratch_path,_scanFolder) )		

	proj.system.Summary()
	#---------------------------------------------
	#setting atoms for scan	
	atomsf = [ 4, 6,  8, 14] 
	#setting parameters
	parameters = { "ATOMS_RC1":atomsf     					,
				   "nsteps_RC1":_nsteps   					,
				   "dincre_RC1":_dincre   					, 
				   "rc_type_1" :"dihedral"					, 
				   "ndim":1               					,
				   "MC_RC1":True          					,
				   "save_format":".dcd"   					,
				   "log_frequency":50.0   					,
				   "Debug":True           					,
				   "simulation_type":"Relaxed_Surface_Scan" ,
				   "force_constant":100.0                   }
    #run the simulation
    #---------------------------------------------------------------------
	proj.Run_Simulation(parameters)
	proj.SaveSystem()
	proj.FinishRun()
#=====================================================
def Scan2D_Dihedral(_xnsteps,_ynsteps,_dincreX=10.0,_dincreY=10.0,name="Default"):
	#---------------------------------------------	
	_scanFolder = "SCAN2D_dihedral"
	if not name == "Default": _scanFolder = name
	
	proj=SimulationProject.From_PKL( os.path.join(balapkl)					, 
									 os.path.join(scratch_path,_scanFolder) )		
	proj.system.Summary()
	#---------------------------------------------
	#setting atoms for scan	
	atomsf = [ 4, 6,  8, 14] 
	atomss = [ 6, 8, 14, 16]
	#setting parameters
	parameters = { "ATOMS_RC1":atomsf     ,
				   "ATOMS_RC2":atomss     ,
				   "contour_lines":12     ,
				   "nsteps_RC1":_xnsteps  ,
				   "nsteps_RC2":_ynsteps  ,
				   "dincre_RC1":_dincreX  ,
				   "dincre_RC2":_dincreY  ,
				   "rc_type_1" :"dihedral", 
				   "rc_type_2" :"dihedral", 
				   "ndim":2               ,
				   "force_constant_1":20.0,
				   "force_constant_2":20.0,
				   "log_frequency":100    ,
				   "simulation_type":"Relaxed_Surface_Scan",
				   "NmaxThreads":8       }
    #run the simulation
    #---------------------------------------------------------------------
	proj.Run_Simulation(parameters)
	proj.SaveSystem()		
	proj.FinishRun()
#=====================================================
def FreeEnergyDihedral1D(nsteps):
	proj=SimulationProject.From_PKL(balapkl, _FolderName=os.path.join(scratch_path,"FE_dihedral_1D") )		
	#-------------------------------------------------
	_name = "SCAN1D_4FEcalculations_dihedral"
	_path = os.path.join( os.path.join(scratch_path,_name,"ScanTraj.ptGeo") )
	if not os.path.exists(_path ): Scan1D_Dihedral(20,18.0, name=_name )
	#-------------------------------------------------	
	atomsf = [4, 6,  8, 14]	
	rc1 = ReactionCoordinate(atomsf,False,0)
	rc1.GetRCLabel(proj.system)	
	rc1.SetInformation(proj.system,0.0)

	#-------------------------------------------------	
	parameters = { "ATOMS_RC1":atomsf			  ,
				   "ndim": 1 					  ,
				   "sampling_factor":nsteps/10	  ,
				   "rc_type_1":"dihedral"         ,
				   "equilibration_nsteps":nsteps/2,
				   "force_constant_1":30.0        ,
				   "production_nsteps":nsteps	  ,
				   "source_folder":_path 		  ,
				   "MD_method":"LeapFrog"		  ,
				   "simulation_type":"Umbrella_Sampling"}
	#-------------------------------------------------
	#RUN umbrella sampling
	proj.Run_Simulation(parameters)	
	#-------------------------------------------------
#=====================================================
def FreeEnergyDihedral2D(nsteps):
	'''
	'''
	proj=SimulationProject.From_PKL( balapkl, _FolderPath=os.path.join(scratch_path,"FE_2D_dihedral") )		
	
	atomsf = [ 4, 6,  8, 14 ] 
	atomss = [ 6, 8, 14, 16 ]
	
	rc1 = ReactionCoordinate(atomsf,False,0)
	rc1.GetRCLable(proj.system)	
	rc2 = ReactionCoordinate(atomss,False,0)
	rc2.GetRCLabel(proj.system)	
	
	_name = "SCAN2D_4FEcalculations_dihedral"
	_path = os.path.join( scratch_path,_name,"ScanTraj.ptGeo")
	if not os.path.exists(_path): Scan2D_Dihedral(12,12,name=_name)

	parameters = { "ATOMS_RC1":atomsf				,
				   "ATOMS_RC2":atomss				,
				   "ndim": 2 						,
				   "sampling_factor":nsteps/10		,
				   "equilibration_nsteps":nsteps/2 	,
				   "production_nsteps":nsteps		,
				   "source_folder":_path 			,
				   "force_constant_1":25.0          ,
				   "force_constant_2":25.0          ,
				   "rc_type_1":"dihedral"           ,
				   "rc_type_2":"dihedral"           ,
				   "MD_method":"LeapFrog"			,				   
				   "simulation_type":"Umbrella_Sampling",
				   "NmaxThreads":NmaxThreads		}
	
	proj.Run_Simulation(parameters)
	proj.FinishRun()

#===============================================================================
def Path_From_PES():
	'''
	'''
	if not os.path.exists( os.path.join(scratch_path,"QCMM_Scan2D_simple_distance") ):
		QCMMScan2DsimpleDistance(12,12,0.15,0.15)
	log_path = os.path.join(scratch_path,"QCMM_Scan2D_simple_distance","ScanTraj.log")
	parameters= {"xsize":12,"ysize":12,"type":"2D","log_name":log_path,
	"crd1_label":rc1_sd.label,"crd2_label":rc2_sd.label,
	"contour_lines":10,"analysis_type":"Energy_Plots","in_point":[0,0],"fin_point":[11,11]}
	proj.Run_Analysis(parameters)

#===============================================================================
def ReacCoordSearchers(_type):	
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts","7tim_#4.pkl") ):
		QCMM_optimizations()
	proj=SimulationProject( os.path.join(scratch_path,"ReactionPathsSearchers") )		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts","7tim_#4.pkl") )
	#-------------------------------------------------------------------------------
	#generate initial and final coordinates for NEB 
	#generate trajectory for SAW
	#-------------------------------------------------------------------------------
	_name = "SCAN1D_4NEB_and_SAW"
	_path = os.path.join( os.path.join(scratch_path,_name,"ScanTraj.ptGeo") )
	if not os.path.exists(_path):
		QCMMScanMultipleDistance(16,0.09,name=_name)
	init_path    = os.path.join( _path, "frame0.pkl")
	final_path   = os.path.join( _path, "frame15.pkl")
	saddle_coord = os.path.join( _path, "frame9.pkl") 

	#-------------------------------------------------------------------------------
	if _type == "NEB":
		parameters = {  "init_coord":init_path           					,
						"final_coord":final_path         					,
						"traj_bins":16                   					,
						"refine_methods":["rm1","am1","pm3","pddgpm3","pm6"],
						"RMS_growing_intial_string":1.0  					,
						"simulation_type":"NEB"          					,
						"spring_force_constant":500.0    					,
						"rmsGradient":0.10               					,
						"fixed_terminal_images":True    	                }
	#-------------------------------------------------------------------------------
	elif _type == "NEB_traj":
		parameters = {  "traj_source":_path,
						"traj_bins":16                   					,
						"refine_methods":["rm1","am1","pm3","pddgpm3","pm6"],
						"RMS_growing_intial_string":1.0  					,
						"simulation_type":"NEB"          					,
						"spring_force_constant":500.0    					,
						"rmsGradient":0.3               					,
						"fixed_terminal_images":True                       }
	#------------------------------------------------------------------------------
	elif _type == "SAW":
		parameters = {  "traj_source":_path,
						"traj_bins":16                   					,
						"refine_methods":["rm1","am1","pm3","pddgpm3","pm6"],
						"simulation_type":"SAW"          					,
						"rmsGradient":0.30                                  }
	#------------------------------------------------------------------------------					
	elif _type == "BakerSaddle":
		parameters = {  "saddle_coord":saddle_coord           			    ,
						"simulation_type":"Baker_Saddle"          			,
						"rmsGradient":0.10               					}
	#------------------------------------------------------------------------------
	elif _type == "SteepestDescent":
		parameters = {  "traj_source":_path                                  ,
						"refine_methods":["rm1","am1","pm3","pddgpm3","pm6"] ,
						"simulation_type":"Steep_Path_Searcher"          	 ,
						"rmsGradient":0.10    								 ,
						"function_step":0.025                                ,
						"mass_weighting":True                                ,
						"path_step":2.0          		              		}
	#------------------------------------------------------------------------------
	proj.Run_Simulation(parameters)
#================================================================================
def NEB_FreeEnergy():
	if not os.path.exists( os.path.join(scratch_path,"QCMMopts.pkl") ):
		QCMM_optimizations()
	proj=SimulationProject( os.path.join(scratch_path,"NEB_FreeEnergy") )		
	proj.LoadSystemFromSavedProject( os.path.join(scratch_path,"QCMMopts.pkl") )
	
	_name = "ReactionPathsSearchers"
	_path = os.path.join( os.path.join(scratch_path,_name,"ReactionPath.ptGeo") )

	_type = "NEB" 
	if not os.path.exists(_path): ReacCoordSearchers(_type)
	
	atom1 = AtomSelection.FromAtomPattern(proj.system,"*:LIG.*:C02")
	atom2 = AtomSelection.FromAtomPattern(proj.system,"*:LIG.*:H02")
	atom3 = AtomSelection.FromAtomPattern(proj.system,"*:GLU.164:OE2")	
	atom6 = AtomSelection.FromAtomPattern(proj.system,"*:LIG.*:O06")
	atom5 = AtomSelection.FromAtomPattern(proj.system,"*:HIE.94:HE2")
	atom4 = AtomSelection.FromAtomPattern(proj.system,"*:HIE.94:NE2")
	atomsf = [ atom1[0], atom2[0], atom3[0] ] 
	atomss = [ atom4[0], atom5[0], atom6[0] ]
	
	rc1 = ReactionCoordinate(atomsf,False,0)
	rc1.GetRCLabel(proj.system)	
	rc2 = ReactionCoordinate(atomss,False,0)	
	rc2.GetRCLabel(proj.system)

	USparameters = { "ATOMS_RC1":atomsf				,
				   "ATOMS_RC2":atomss				,
				   "ndim":2 						,
				   "sampling_factor":500		    ,
				   "equilibration_nsteps":2500 	    ,
				   "production_nsteps":5000		    ,
				   "source_folder":_path 			,
				   "MD_method":"LeapFrog"			,
				   "MC_RC1":True					,
				   "MC_RC2":True					,
				   "simulation_type":"Umbrella_Sampling",
				   "NmaxThreads":NmaxThreads		}

	proj.Run_Simulation(USparameters)

	_path = os.path.join( scratch_path, "NEB_FreeEnergy")	
	PMFparameters = { "source_folder":_path,
				   "xnbins":10           ,
				   "ynbins":10           ,
				   "ywindows":0          ,
				   "xwindows":16         ,
				   "crd1_label":rc1.label,
				   "crd2_label":rc2.label,
				   "oneDimPlot":True     ,
				   "simulation_type":"PMF_Analysis",
				   "temperature":300.15	 }
	#RUN WHAM, calculate PMF and free energy
	proj.Run_Simulation(PMFparameters)

#=====================================================
def write_qm_log():
	proj = SimulationProject( os.path.join(scratch_path, "QMlog") )
	if not os.path.exists(os.path.join(scratch_path, "QMlog") ):
		os.makedirs( os.path.join(scratch_path, "QMlog") )
	proj.LoadSystemFromSavedProject( balapkl )
	Pickle(balapkl[:-4]+"crd", proj.system.coordinates3)
	qcModel = QCModelMNDO.WithOptions( hamiltonian = "am1" )
	proj.system.DefineQCModel(qcModel)

	proj.system.Energy()
	test = WriteQMLog(proj.system,os.path.join(scratch_path, "test.log") )
	test.write()
	
	_mopacKeys = ["AUX", "LARGE"]	
	mop = MopacQCMMinput(proj.system,os.path.join(scratch_path, "QMlog"),balapkl[:-4]+"crd",_mopacKeys,"am1")
	mop.CalculateGradVectors()
	mop.write_input(0,1)
	mop.Execute()

#=====================================================
if __name__ == "__main__":	
	#------------------------------------
	if len(sys.argv) > 2:
		if sys.argv[2] == "-np": NmaxThreads = int(sys.argv[3])	
	#------------------------------------------------
	if sys.argv[1] == "help":
		help_text = "Options for testn\n"
		help_text+= "\t1:Molecular Dynamics Algorithms\n\t2:Molecular Dynamics heating protocol\n\t3:QC/MM energies"
		help_text+= "\n\t4:QC/MM DFTB Energy\n\t5:QC/MM ORCA Energy\n\t6:QC/MM optimizations\n\t7:QC/MM molecular dynamics\n\t8:QC/MM restricted molecular dynamics"
		help_text+= "\n\t9:Simple distance 1D scan\n\t10:Multiple distance 1D scan\n\t11:Simple distance 2D scan\n\t12:Mixed distance 2D scan\n\t13:Multiple distance 2D scan"
		help_text+= "\n\t14:Adaptative 2D scan\n\t15:Scan 1D dihedral\n\t16:Scan 2D dihedral\n\t17:Free Energy 1D simple distance\n\t18:Free Energy 1D multiple distance"
		help_text+= "\n\t19:Free Energy 1D dihedral\n\t20:Free Energy 1D dihedral with optimization\n\t21:Umbrella sampling restart test\n\t22:Free energy simple distance 2D"
		help_text+= "\n\t23:Free energy simple mixed distance 2D\n\t24:Free energy multiple distance 2D\n\t25:Semiempirical in pDynamo energy refinement\n\t26:Energy Plots analysis"
		help_text+= "\n\t27:Raction coordinates searching\n\t28:Semiempirical mopac energy refinement\n\t29:Semiempirical in pDynamo energy 2D\n\t30:Trajectory Analysis plots"
		help_text+= "\n\t31:pDynamo internal refinement\n\t32:ORCA Energy refinement\n\t33:Change QC Region tests"
		print(help_text)			
	#------------------------------------------------

	elif int(sys.argv[1]) == 0:	 SetMMsytem()
	elif int(sys.argv[1]) == 1:  MMMD_Algorithms()
	elif int(sys.argv[1]) == 2:	 MMMD_Heating()
	elif int(sys.argv[1]) == 3:	 QCMM_SMO_Energies()
	elif int(sys.argv[1]) == 4:  QCMM_DFTBplus()	
	elif int(sys.argv[1]) == 5:  QCMM_Orca()
	elif int(sys.argv[1]) == 6:  QCMM_optimizations()
	elif int(sys.argv[1]) == 7:  QCMM_MD()
	elif int(sys.argv[1]) == 8:  QCMM_MDrestricted()
	elif int(sys.argv[1]) == 9:  QCMMScanSimpleDistance(14,0.1)
	elif int(sys.argv[1]) == 10: QCMMScanMultipleDistance(14,0.1)
	elif int(sys.argv[1]) == 11: QCMMScan2DsimpleDistance(10,10,0.2,0.2)
	elif int(sys.argv[1]) == 12: QCMMScan2DmixedDistance(10,10,0.2,0.2)	
	elif int(sys.argv[1]) == 13: QCMMScan2DmultipleDistance(10,10,0.2,0.2)
	elif int(sys.argv[1]) == 14: QCMMScans2D_Adaptative(10,10,0.2,0.2)
	elif int(sys.argv[1]) == 15: Scan1D_Dihedral(18,20.0)
	elif int(sys.argv[1]) == 16: Scan2D_Dihedral(12,12,20.0,20.0)
	elif int(sys.argv[1]) == 17: FreeEnergy1DSimpleDistance(600)
	elif int(sys.argv[1]) == 18: FreeEnergy1DMultipleDistance(600)
	elif int(sys.argv[1]) == 19: FreeEnergyDihedral1D(5000)
	elif int(sys.argv[1]) == 20: FreeEnergy1DSimpleDistanceOPT(500)
	elif int(sys.argv[1]) == 21: UmbrellaSampling1Drestart(500)
	elif int(sys.argv[1]) == 22: FreeEnergy2DsimpleDistance(500)
	elif int(sys.argv[1]) == 23: FreeEnergy2DmixedDistance(500)
	elif int(sys.argv[1]) == 24: FreeEnergy2DmultipleDistance(500)
	elif int(sys.argv[1]) == 25: pDynamoEnergyRef_1D()
	elif int(sys.argv[1]) == 26: EnergyAnalysisPlots()
	elif int(sys.argv[1]) == 27: ReacCoordSearchers("NEB")
	elif int(sys.argv[1]) == 28: MopacEnergyRef()
	elif int(sys.argv[1]) == 29: pDynamoEnergyRef_2D()
	elif int(sys.argv[1]) == 30: TrajectoryAnalysisPlots()
	elif int(sys.argv[1]) == 31: pDynamoEnergyRef_abInitio()
	elif int(sys.argv[1]) == 32: ORCAEnergy_ref()
	elif int(sys.argv[1]) == 33: Change_QC_Region()
	elif int(sys.argv[1]) == 34: NEB_FreeEnergy()
	elif int(sys.argv[1]) == 35: write_qm_log()
	else: pass  


