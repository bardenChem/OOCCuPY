#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamoWrapper import Wrapper
import SimulationSystem 
import os, sys

folder = os.path.join("Tests","pDynamoWrapper","test_03")

#===================================
def info():
	print_message = "OOCCuPy pDynamoWrapper Libray test #03:\n\t "
	print_message += "Test geometry optimization algorithms.\n"

	print(print_message)
#------------------------------------
def Run_Test():
	'''
	Test geometry optimization algorithms
	'''
	info()
	algs = ["ConjugatedGradient",
			"LFBGS"             ,
			"SteepestDescent"   ,
			#"QuasiNewton"       ,
			"FIRE"              ]
	
	if not os.path.exists( os.path.join("Tests","pDynamoWrapper","test_01","7tim.pkl") ):
		try: os.system("python3 Tests/pDynamoWrapper/test_01.py")
		except: 
			print("There is no input file for this example! Run example #01!")
			return(False)

	_parameters = {
		"Input_Type":"pkl",
		"pkl_file":"Tests/pDynamoWrapper/test_01/7tim_pruned_and_fix.pkl",
		"simulation_type":"Geometry_Optimization",
		"save_format":".dcd",
		"save_frequency":20,
		"rmsGradient":0.1,
		"log_frequency":10,
		"maxIterations":2200,
	}
	#------------------------------------
	test_01 = Wrapper(folder)
	for alg in algs:
		test_01.Set_System(_parameters)
		_parameters["optmizer"]=alg
		_parameters["trajectory_name"]="7timMMopt_"+alg+".ptGeo"
		test_01.Run_Simulation(_parameters)
		test_01.SaveSystem("7tim_opt"+alg)
	#-----------------------------------
	#QC/MM optimization
	if not os.path.exists( os.path.join("Tests","pDynamoWrapper","test_18","7tim_qcmm_rm1_pruned") ):
		try: os.system("python3 Test/pDynamoWrapper/test_18.py")
		except: 
			print("There is no input file for this example! Run example #18!")
			return(False)

	_parameters["pkl_file"] = "Tests/pDynamoWrapper/test_18/7tim_qcmm_rm1_pruned.pkl"
	_parameters["optmizer"] = "ConjugatedGradient"
	_parameters["trajectory_name"]="7timqcmmm.ptGeo"
	test_02 = Wrapper(folder)
	test_02.Set_System(_parameters)
	test_02.Run_Simulation(_parameters)

	test_02.SaveSystem("7tim_qcmm_opt_PF")
	
	
#===================================
if __name__ == '__main__': 
	Run_Test()