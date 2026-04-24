#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamoWrapper import Wrapper
import os, sys

from config import get_config
config = get_config()
ooccupy_root = get_config().get_ooccupy_root()
folder = os.path.join(ooccupy_root, "Tests", "pDynamoWrapper", "test_03")

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
	
	# Get OOCCuPY root from config
	config = get_config()
	ooccupy_root = config.get_ooccupy_root()
	
	# Build paths using ooccupy_root
	folder = os.path.join(ooccupy_root, "Tests", "pDynamoWrapper", "test_03")
	test_01_folder = os.path.join(ooccupy_root, "Tests", "pDynamoWrapper", "test_01")
	test_18_folder = os.path.join(ooccupy_root, "Tests", "pDynamoWrapper", "test_18")

	algs = ["ConjugatedGradient",
			"LFBGS"             ,
			"SteepestDescent"   ,
			#"QuasiNewton"       ,
			"FIRE"              ]
	
	if not os.path.exists( os.path.join(test_01_folder, "7tim_pruned_and_fix.pkl") ):
		try: os.system("python3 Tests/pDynamoWrapper/test_01.py")
		except: 
			print("There is no input file for this example! Run example #01!")
			return(False)

	_parameters = {
		"Input_Type":"pkl",
		"pkl_file":os.path.join(test_01_folder, "7tim_pruned_and_fix.pkl"),
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
	if not os.path.exists( os.path.join(test_18_folder, "7tim_qcmm_rm1_pruned.pkl") ):
		try: os.system("python3 Tests/pDynamoWrapper/test_18.py")
		except: 
			print("There is no input file for this example! Run example #18!")
			return(False)

	_parameters["pkl_file"] = os.path.join(test_18_folder, "7tim_qcmm_rm1_pruned.pkl")
	_parameters["optmizer"] = "ConjugatedGradient"
	_parameters["trajectory_name"]="7timqcmmm.ptGeo"
	test_02 = Wrapper(folder)
	test_02.Set_System(_parameters)
	test_02.Run_Simulation(_parameters)

	test_02.SaveSystem("7tim_qcmm_opt_PF")
	
	
#===================================
if __name__ == '__main__': 
	Run_Test()