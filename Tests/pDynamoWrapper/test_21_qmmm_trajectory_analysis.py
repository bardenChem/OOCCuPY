#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from pDynamoWrapper import Wrapper
import os, sys

from config import get_config
config = get_config()
ooccupy_root = config.get_ooccupy_root()
folder = os.path.join(ooccupy_root, "Tests", "pDynamoWrapper", "test_21")
folder05 = os.path.join(ooccupy_root,"Tests","pDynamoWrapper","test_05")
folder08 = os.path.join(ooccupy_root,"Tests","pDynamoWrapper","test_08")


#===================================
def info():
	print_message =  "OOCCuPy pDynamoWrapper Libray test #21:\t "
	print_message += "Testing analysis of a run of Molecular dynamics hybrid simulations.\n"

	print(print_message)
#------------------------------------
def Run_Test():
	'''
	Test molecular dynamics analysis tools
	'''
	info()
	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join(folder05,"qcmm_optam1","7tim_am1_opt_PF.pkl"),
		"set_reaction_crd":2,	
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		"atoms_rc2":["*:LIG.*:O06","*:HIE.94:HE2","*:HIE.94:NE2"],
		"type_rc1":"Distance",
		"type_rc2":"Distance",
		"mass_constraints":["no","no"],
		"temperature": 293.15,
		"ndim":2,
		"analysis_type":"Trajectory_Analysis",
		"nsteps":2000,
		"source_folder":os.path.join( folder08, "7timQCMD_restrictedproduction.ptGeo"),
		"calculate_distances":"yes",
	}	
	#------------------------------------
	#protocol production
	test_01 = Wrapper(folder)
	test_01.Set_System(system_parameters)
	test_01.Run_Analysis(system_parameters)
	test_01.SaveSystem()
	#-----------------------------------
	
	
#===================================
if __name__ == '__main__': 
	Run_Test()
