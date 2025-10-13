#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamo_Scripts import Scripts
import SimulationSystem 
import os, sys
#====================================================
def info():
	print_message =  "OOCCuPy pDynamo_Scripts Libray test #13:\t "
	print_message += "Testing the setting and run of Nudget Elastic Band routines.\n"

	print(print_message)
#----------------------------------
def Run_Test():
	'''
	'''
	info()
	_path = "test_05/Multiple_Distance_rm1/ScanTraj.ptGeo"

	init_path    = os.path.join( _path, "frame0.pkl")
	final_path   = os.path.join( _path, "frame19.pkl")
	saddle_coord = os.path.join( _path, "frame12.pkl") 

	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_05","qcmm_optam1","7tim_am1_opt_PF.pkl"),		
		"set_reaction_crd":2,	
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		"atoms_rc2":["*:LIG.*:O06","*:HIE.94:HE2","*:HIE.94:NE2"],
		"type_rc1":"Distance",
		"type_rc2":"Distance",
		"mass_constraint":["yes","yes"]
	}

	parameters_NEB = {  "init_coord":init_path           					,
						"final_coord":final_path         					,
						"traj_bins":12                   					,
						"refine_methods":["rm1","am1","pm3","pddgpm3","pm6"],
						"RMS_growing_intial_string":1.0  					,
						"simulation_type":"NEB"          					,
						"spring_force_constant":800.0    					,
						"rmsGradient":0.10               					,
						"fixed_terminal_images":"no"    	                }

	test_01 = Scripts("test_13_NEB")
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(parameters_NEB)
	test_01.SaveSystem()

	parameters_NEBT = { "traj_source":_path,
						"traj_bins":20                   					,
						"refine_methods":["rm1","am1","pm3","pddgpm3","pm6"],
						"RMS_growing_intial_string":1.0  					,
						"simulation_type":"NEB"          					,
						"spring_force_constant":800.0    					,
						"rmsGradient":0.3               					,
						"fixed_terminal_images":"no"                       }
	
	test_02 = Scripts("test_13_NEB_traj")
	test_02.Set_System(system_parameters)
	test_02.Run_Simulation(parameters_NEBT)
	test_02.SaveSystem()

	
	#------------------------------------------------------------------------------		
	parameters_steep = { "traj_source":_path                               ,
						"refine_methods":["rm1","am1","pm3","pddgpm3","pm6"],
						"simulation_type":"Steep_Path_Searcher"          	,
						"rmsGradient":0.10    								,
						"function_step":0.025                               ,
						"traj_bins":20                                      ,
						"mass_weighting":True                               ,
						"path_step":2.0          		              		}

	test_03 = Scripts("test_13_Steep")
	test_03.Set_System(system_parameters)
	#test_03.Run_Simulation(parameters_steep)
	test_03.SaveSystem()		
	#------------------------------------------------------------------------------
	parameters_saddle = {"saddle_coord":saddle_coord           			,
						"simulation_type":"Baker_Saddle"          		,
						"rmsGradient":0.50               				}
	test_04 = Scripts("test_13_Saddle")
	test_04.Set_System(system_parameters)
	#test_04.Run_Simulation(parameters_saddle)
	test_04.SaveSystem()				

#===================================
if __name__ == '__main__': 
	if ( sys.argv[1] ) == "-print":	info()
	else: Run_Test()
