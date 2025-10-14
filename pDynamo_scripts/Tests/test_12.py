#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamoWrapper import Scripts
import SimulationSystem 
import os, sys
#===================================
def info():
	print_message =  "OOCCuPy pDynamoWrapper Libray test #12:\n\t "
	print_message += "Testing the setting and run of bidimensional Umbrella Sampling.\n"

	print(print_message)
#----------------------------------
def Run_Test():
	'''
	Test umbrella sampling 2D
	'''
	info()
	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_05","qcmm_optam1","7tim_am1_opt_PF.pkl"),		
		"set_reaction_crd":2,	
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		"atoms_rc2":["*:HIE.94:NE2","*:HIE.94:HE2","*:LIG.*:O06"],
		"type_rc1":"Distance",
		"type_rc2":"Distance",
		"mass_constraint":["yes","yes"],
	}

	_path   = "test_06/Multiple_Distance_rm1/ScanTraj.ptGeo"
	
	US_parameters = {
				   "ndim": 2 					        ,
				   "sampling_production":100	   		,
				   "equilibration_nsteps":2000			,
				   "production_nsteps":20000   			,
				   "source_folder":_path 		       	,
				   "MD_method":"LeapFrog"		      	,
				   "simulation_type":"Umbrella_Sampling",
				   "NmaxThreads":50		  				,
				   "xnbins":10                          ,
				   "ynbins":10                          ,
				   "xwindows":20                        ,
				   "ywindows":20                        , 
				   "trajectory_name":"US_test"          ,
				   }
	#------------------------------------
	test_01 = Scripts("test_12")
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(US_parameters)
	test_01.SaveSystem()
	#-----------------------------------	
	
	
#===================================
if __name__ == '__main__': 
	if ( sys.argv[1] ) == "-print":	info()
	else: Run_Test()
