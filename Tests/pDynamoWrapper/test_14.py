#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamoWrapper import Wrapper
import os, sys

from config import get_config
config = get_config()
ooccupy_root = config.get_ooccupy_root()
folder = os.path.join(ooccupy_root,"Tests","pDynamoWrapper","test_14")


#====================================================
def info():
	print_message =  "OOCCuPy pDynamoWrapper Libray test #14:\n\t "
	print_message += "Testing the setting and run of Normal_Modes.\n"

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
		"set_initial_crd":init_path			
	}

	parameters_NM = {  "init_coord":init_path           		,						
						"simulation_type":"Normal_Modes"     	,
						"cycles":12,	
						"temperature":298.15,
						"mode":0}

	test_01 = Wrapper("folder")
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(parameters_NM)
	test_01.SaveSystem()

#====================================================
if __name__ == '__main__': 
	Run_Test()