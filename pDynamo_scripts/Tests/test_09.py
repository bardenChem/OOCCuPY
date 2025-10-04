#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamo_Scripts import Scripts
import SimulationSystem 
import os, sys
#===================================
def Run_Test():
	'''
	Test internal energy refinement
	'''

	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_05","qcmm_optam1","7tim_am1_opt_PF.pkl"),
		"set_reaction_crd":1,	
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		#"atoms_rc2":["*:LIG.*:O06","*:HIE.94:HE2","*:HIE.94:NE2"],
		"type_rc1":"Distance",
		#"type_rc2":"Distance",
		"mass_constraints":["yes","yes"],
	}

	_path   = "test_05/Multiple_Distance_rm1/ScanTraj.ptGeo"
	methods = ["am1","pm3","rm1","pm6"]
	
	simulation_parameters = { "xnbins":20			    ,
				   "source_folder":_path                , 
				   "folder":"test_09"                   ,
				   "QCcharge":-3		                    ,
				   "multiplicity":1 	                ,
				   "methods_lists":methods              ,					   
				   "NmaxThreads":10                     ,
				   "simulation_type":"Energy_Refinement",
				   "Software":"pDynamo"	}				  
					
	#------------------------------------
	test_01 = Scripts("test_09")
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(simulation_parameters)
	test_01.SaveSystem()
	#-----------------------------------
	methods.append("pm7")
	
	simulation_parameters["Software"] = "mopac"
	simulation_parameters["mopac_keywords"] = [] 
	simulation_parameters["folder"] = "test_09_mopac"

	test_02 = Scripts("test_09_mopac")
	test_02.Set_System(system_parameters)
	test_02.Run_Simulation(simulation_parameters)
	test_02.SaveSystem()

	simulation_parameters["Software"] = "pySCF"
	simulation_parameters["folder"]   = "test_09_pyscf"
	simulation_parameters["pySCF_method"] = "RKS"
	simulation_parameters["functional"] = "b3lyp"   	                                         
	simulation_parameters["basis"]      = "6-31G*" 

	test_03 = Scripts("test_09_pyscf")
	test_03.Set_System(system_parameters)
	test_03.Run_Simulation(simulation_parameters)

	simulation_parameters["Software"]    = "ORCA"
	simulation_parameters["folder"]      = "test_09_orca"
	simulation_parameters["orca_method"] = "b3lyp"   	                                         
	simulation_parameters["basis"]       = "6-31G*" 
	test_04 = Scripts("test_09_orca")
	test_04.Set_System(system_parameters)
	test_04.Run_Simulation(simulation_parameters)
	test_04.SaveSystem()
	
#===================================
if __name__ == '__main__': Run_Test()