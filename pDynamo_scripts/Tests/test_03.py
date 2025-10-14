#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamoWrapper import Scripts
import SimulationSystem 
import os, sys
#===================================
def info():
	print_message = "OOCCuPy pDynamoWrapper Libray test #03:\t "
	print_message += "Test geometry optimization algorithms."

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
			"QuasiNewton"       ,
			"FIRE"              ]
	
	if not os.path.exists( os.path.join("test_01","7tim.pkl") ):
		try: os.system("python3 test_01.py")
		except: 
			print("There is no input file for this example! Run example #01!")
			return(False)

	_parameters = {
		"Input_Type":"pkl",
		"pkl_file":"test_01/7tim_pruned_and_fix.pkl",
		"simulation_type":"Geometry_Optimization",
		"save_format":".dcd",
		"save_frequency":20,
		"rmsGradient":0.1,
		"maxIterations":2200,
	}
	#------------------------------------
	test_01 = Scripts("test_03")
	for alg in algs:
		test_01.Set_System(_parameters)
		_parameters["optmizer"]=alg
		_parameters["trajectory_name"]="7timMMopt_"+alg+".ptGeo"
		test_01.Run_Simulation(_parameters)
		test_01.SaveSystem("7tim_opt"+alg)
	#-----------------------------------
	#QC/MM optimization
	if not os.path.exists( os.path.join("test_01","7tim.pkl") ):
		try: os.system("python3 test_01.py")
		except: 
			print("There is no input file for this example! Run example #01!")
			return(False)

	_parameters["pkl_file"] = "test_02/7tim_qcmm_rm1_pruned.pkl"
	_parameters["optmizer"] = "ConjugatedGradient"
	_parameters["trajectory_name"]="7timqcmmm.ptGeo"
	test_02 = Scripts("test_03")
	test_02.Set_System(_parameters)
	test_02.Run_Simulation(_parameters)

	test_02.SaveSystem("7tim_qcmm_opt_PF")
	
	
#===================================
if __name__ == '__main__': 
	if ( sys.argv[1] ) == "-print":	info()
	else: 						    Run_Test()