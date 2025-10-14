#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamoWrapper import Scripts
import SimulationSystem 
import os, sys
#===================================
def info():
	print_message =  "OOCCuPy pDynamoWrapper Libray test #15:\t "
	print_message += "Testing the setting and run of bidimensional energy refinement with external softwares.\n"

	print(print_message)
#----------------------------------
def Run_Test():
	'''
	Test energy refinement
	'''
	info()
	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_05","qcmm_optam1","7tim_am1_opt_PF.pkl"),
		"set_reaction_crd":2,	
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		"atoms_rc2":["*:LIG.*:O06","*:HIE.94:HE2","*:HIE.94:NE2"],
		"type_rc1":"Distance",
		"mass_constraints":["yes","yes"],
	}

	_path   = "test_06/Multiple_Distance_rm1/ScanTraj.ptGeo"
	methods = ["am1","pm3","rm1","pm6"]
	
	simulation_parameters = { "xnbins":12			    ,
				   "ynbins":12                         ,
				   "source_folder":_path                , 
				   "folder":"test_15"                   ,
				   "QCcharge":-3	                    ,
				   "multiplicity":1 	                ,
				   #"correct_QMMM_charge":"yes",
				   "methods_lists":methods              ,					   
				   "NmaxThreads":12                     ,
				   "reverse_rc2":"yes",
				   "simulation_type":"Energy_Refinement",
				   "Software":"pDynamo"	}				  
					
	#------------------------------------
	test_01 = Scripts("test_15")
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(simulation_parameters)
	test_01.SaveSystem()
	#-----------------------------------
	methods.append("pm7")
	
	simulation_parameters["Software"] = "mopac"
	simulation_parameters["mopac_keywords"] = [] 
	simulation_parameters["folder"] = "test_15_mopac"

	test_02 = Scripts("test_15_mopac")
	test_02.Set_System(system_parameters)
	test_02.Run_Simulation(simulation_parameters)
	test_02.SaveSystem()
	
	simulation_parameters["Software"]    = "ORCA"
	simulation_parameters["folder"]      = "test_10_orca"
	simulation_parameters["orca_method"] = "b3lyp"   	                                         
	simulation_parameters["basis"]       = "6-31G*" 
	test_03 = Scripts("test_10_orca")
	test_03.Set_System(system_parameters)
	#test_03.Run_Simulation(simulation_parameters)
	test_03.SaveSystem()
	
#===================================
if __name__ == '__main__': 
	if ( sys.argv[1] ) == "-print":	info()
	else: Run_Test()