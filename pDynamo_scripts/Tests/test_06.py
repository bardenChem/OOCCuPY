#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamo_Scripts import Scripts
import SimulationSystem 
import os, sys
#==============================================
def info():
	print_message =  "OOCCuPy pDynamo_Scripts Libray test #06:\t "
	print_message += "Testing the setting and run of relaxed bidimensional scans.\n"

	print(print_message)
#-------------------------------------------------
def Simple_Distance2D(_hamiltonian):
	'''
	'''
	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_05","qcmm_opt"+_hamiltonian,"7tim_"+_hamiltonian+"_opt_PF.pkl"),		
		"set_reaction_crd":2,
		"atoms_rc1":["*:LIG.*:H02","*:GLU.164:OE2"],
		"atoms_rc2":["*:LIG.*:O06","*:HIE.94:HE2"],
		"type_rc1":"Distance",
		"type_rc2":"Distance",
		"mass_constraints":["no","no"]
	}		
	scan1_parameters = {
		"simulation_type":"Relaxed_Surface_Scan",
		"dincre_rc1":-0.1,
		"nsteps_rc1":16,
		"dincre_rc2":-0.1,
		"nsteps_rc2":16,
		"maxIterations":2200,
		"log_frequency":10,
		"optmizer":"SteepestDescent",
		"NmaxThreads":8,
		"force_constants":[1000.0,1000.0]
	}
	#test simple distance
	test_01 = Scripts("test_06/Simple_Distance2D_"+_hamiltonian)
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(scan1_parameters)
	test_01.SaveSystem("Simple_DistanceScan2D")


#-------------------------------------------------
def Mixed_Distance2D(_hamiltonian):
	'''
	'''
	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_05","qcmm_opt"+_hamiltonian,"7tim_"+_hamiltonian+"_opt_PF.pkl"),		
		"set_reaction_crd":2,
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		"atoms_rc2":["*:LIG.*:O06","*:HIE.94:HE2"],
		"type_rc1":"Distance",
		"type_rc2":"Distance",
		"maxIterations":2200,
		"mass_constraints":["yes","no"]
	}		
	scan1_parameters = {
		"simulation_type":"Relaxed_Surface_Scan",
		"dincre_rc1":0.1,
		"dincre_rc2":-0.1,
		"optmizer":"SteepestDescent",
		"maxIterations":2200,
		"nsteps_rc1":16,
		"nsteps_rc2":16,
		"log_frequency":10,
		"NmaxThreads":8,
		"force_constants":[1000.0,1000.0]
	}
	#test simple distance
	test_01 = Scripts("test_06/Mixed_Distance_"+_hamiltonian)
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(scan1_parameters)
	test_01.SaveSystem("Mixed_DistanceScan")


#-----------------------------------------------
def Multiple_Distance2D(_hamiltonian):
	'''
	'''
	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_05","qcmm_opt"+_hamiltonian,"7tim_"+_hamiltonian+"_opt_PF.pkl"),		
		"set_reaction_crd":2,
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		"atoms_rc2":["*:HIE.94:NE2","*:HIE.94:HE2","*:LIG.*:O06"],
		"type_rc1":"Distance",
		"type_rc2":"Distance",
		"maxIterations":2200,
		"mass_constraints":["yes","yes"]
	}		
	scan1_parameters = {
		"simulation_type":"Relaxed_Surface_Scan",
		"dincre_rc1":0.1,
		"dincre_rc2":0.1,
		"optmizer":"ConjugatedGradient",
		"maxIterations":2200,
		"nsteps_rc1":14,
		"log_frequency":10,
		"nsteps_rc2":14,
		"NmaxThreads":8,
		"force_constants":[1000.0,1000.0]

	}
	#test simple distance
	test_01 = Scripts("test_06/Multiple_Distance_"+_hamiltonian)
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(scan1_parameters)
	test_01.SaveSystem("Multiple_DistanceScan")
#-----------------------------------------------
def Run_Test():
	'''
	Test 
	'''
	info()
	if not os.path.exists( os.path.join("test_05","7tim.pkl") ):
		Prepare_MM_System()
	
	if not os.path.exists( os.path.join("test_05","7tim_optMM.pkl") ):
		Prepare_Prune_System()	

	Simple_Distance2D("am1")	
	Mixed_Distance2D("am1")	
	Multiple_Distance2D("am1")
		
	
#===================================
if __name__ == '__main__':
	if ( sys.argv[1] ) == "-print":	info()
	else: 						    Run_Test()