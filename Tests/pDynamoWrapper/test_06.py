#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamoWrapper import Wrapper
import os, sys

from config import get_config

config = get_config()
ooccupy_root = config.get_ooccupy_root()
folder05 = os.path.join(ooccupy_root, "Tests", "pDynamoWrapper", "test_05")
folder = os.path.join(ooccupy_root, "Tests", "pDynamoWrapper", "test_06")

#==============================================
def info():
	print_message =  "OOCCuPy pDynamoWrapper Libray test #06:\t "
	print_message += "Testing the setting and run of relaxed bidimensional scans.\n"

	print(print_message)
#-------------------------------------------------
def Simple_Distance2D(_hamiltonian):
	'''
	'''
	# Get OOCCuPY root from config
	

	if not os.path.exists( os.path.join(folder05,"qcmm_opt"+_hamiltonian,"7tim_"+_hamiltonian+"_opt_PF.pkl") ):
		os.system("python3 Tests/pDynamoWrapper/test_05.py")

	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join(folder05,"qcmm_opt"+_hamiltonian,"7tim_"+_hamiltonian+"_opt_PF.pkl"),		
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
		"nsteps_rc2":-1,
		"maxIterations":2200,
		"log_frequency":10,
		"restart":"yes",
		"optmizer":"SteepestDescent",
		"NmaxThreads":8,
		"force_constants":[1200.0,1200.0]
	}
	#test simple distance
	test_01 = Wrapper( os.path.join(folder,"Simple_Distance2D_"+_hamiltonian) )
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(scan1_parameters)
	test_01.SaveSystem("Simple_DistanceScan2D")


#-------------------------------------------------
def Mixed_Distance2D(_hamiltonian):
	'''
	'''
	if not os.path.exists( os.path.join(folder05,"qcmm_optm"+_hamiltonian+"_opt_PF.pkl") ):
		os.system("python3 test_05.py")
		
	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join(folder05,"qcmm_opt"+_hamiltonian,"7tim_"+_hamiltonian+"_opt_PF.pkl"),		
		"set_reaction_crd":2,
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		"atoms_rc2":["*:LIG.*:O06","*:HIE.94:HE2"],
		"type_rc1":"Distance",
		"type_rc2":"Distance",
		"maxIterations":2200,
		"mass_constraints":["no","no"]
	}		
	scan1_parameters = {
		"simulation_type":"Relaxed_Surface_Scan",
		"dincre_rc1":0.1,
		"dincre_rc2":-0.1,
		"optmizer":"SteepestDescent",
		"maxIterations":2200,
		"nsteps_rc1":-1,
		"nsteps_rc2":-1,
		"restart":"yes",
		"log_frequency":10,
		"NmaxThreads":8,
		"force_constants":[1200.0,1200.0]
	}
	#test simple distance
	test_01 = Wrapper( os.path.join( folder,"Mixed_Distance_"+_hamiltonian) )
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(scan1_parameters)
	test_01.SaveSystem("Mixed_DistanceScan")


#-----------------------------------------------
def Multiple_Distance2D(_hamiltonian):
	'''
	'''
	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join(folder05,"qcmm_opt"+_hamiltonian,"7tim_"+_hamiltonian+"_opt_PF.pkl"),		
		"set_reaction_crd":2,
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		"atoms_rc2":["*:HIE.94:NE2","*:HIE.94:HE2","*:LIG.*:O06"],
		"type_rc1":"Distance",
		"type_rc2":"Distance",
		"maxIterations":2200,
		"mass_constraints":["no","no"]
	}		
	scan1_parameters = {
		"simulation_type":"Relaxed_Surface_Scan",
		"dincre_rc1":0.1,
		"dincre_rc2":0.1,
		"optmizer":"ConjugatedGradient",
		"maxIterations":2200,
		"nsteps_rc1":-1,
		"restart":"yes",
		"log_frequency":10,
		"nsteps_rc2":-1,
		"NmaxThreads":8,
		"force_constants":[1200.0,1200.0]

	}
	#test simple distance
	test_01 = Wrapper( os.path.join(folder,"Multiple_Distance_"+_hamiltonian) )
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(scan1_parameters)
	test_01.SaveSystem("Multiple_DistanceScan")
#-----------------------------------------------
def Run_Test():
	'''
	Test 
	'''
	info()
	if not os.path.exists( os.path.join(folder05,"7tim.pkl") ):
		Prepare_MM_System()
	
	if not os.path.exists( os.path.join(folder05,"7tim_optMM.pkl") ):
		Prepare_Prune_System()	

	Simple_Distance2D("am1")	
	Mixed_Distance2D("am1")	
	Multiple_Distance2D("am1")
		
	
#===================================
if __name__ == '__main__':
	Run_Test()