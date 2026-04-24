#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamoWrapper import Wrapper
import os, sys

from config import get_config

config = get_config()
ooccupy_root = config.get_ooccupy_root()
folder05 = os.path.join(ooccupy_root, "Tests", "pDynamoWrapper", "test_05")
folder = os.path.join(ooccupy_root, "Tests", "pDynamoWrapper", "test_07")

#===================================
def info():
	print_message =  "OOCCuPy pDynamoWrapper Libray test #07:\t "
	print_message += "Testing the setting and run of Molecular dynamics simulations.\n"

	print(print_message)
#-------------------------------------------------
def Run_Test():
	'''
	Test molecular dynamics algorithms with qmmm
	'''
	info()
	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join(folder05,"qcmm_optam1","7tim_am1_opt_PF.pkl"),
	}

	simulation_parameters = {
				  "temperature": 315.15,
				  "simulation_type":"Molecular_Dynamics",
				  "equilibration_nsteps":5000,
				  "production_nsteps":10000,
				  "heating_nsteps":2000,
				  "sampling_equilibration":100,
				  "sampling_production":50,
				  "sampling_heating":50,
				  "log_frequency":10
				}
	
	#------------------------------------
	#protocol production
	test_01 = Wrapper(folder)
	test_01.Set_System(system_parameters)
	simulation_parameters["trajectory_name"]="7timQCMD"
	test_01.Run_Simulation(simulation_parameters)
	test_01.SaveSystem()
	#-----------------------------------
	
	
#===================================
if __name__ == '__main__': 
	Run_Test()