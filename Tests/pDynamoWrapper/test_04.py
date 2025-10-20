#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamoWrapper import Wrapper
import SimulationSystem 
import os, sys


folder = os.path.join("Tests","pDynamoWrapper","test_04")

#===================================
def info():
	print_message =  "OOCCuPy pDynamoWrapper Libray test #04:\n\t "
	print_message += "Test molecular dynamics algorithms.\n"

	print(print_message)
#-----------------------------------
def Run_Test():
	'''
	Test molecular dynamics algorithms
	'''
	info()
	integrators = ["Verlet", "LeapFrog", "Langevin"]	

	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("Tests","pDynamoWrapper","test_03","7tim_optLFBGS.pkl"),
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
				}
	
	if not os.path.exists( os.path.join("Tests","pDynamoWrapper","test_03","7tim_optLFBGS.pkl") ):
		try: os.system("python3 Tests/pDynamoWrapper/test_03.py")
		except: 
			print("There is no input file for this example! Run example #03!")
			return(False)

	#------------------------------------
	#protocol production
	test_01 = Wrapper(folder)
	for alg in integrators:
		test_01.Set_System(system_parameters)
		simulation_parameters["MD_method"]=alg
		simulation_parameters["trajectory_name"]="7timMD_"+alg+".ptGeo"
		test_01.Run_Simulation(simulation_parameters)
		test_01.SaveSystem("7timMD_"+alg+".pkl")
	#-----------------------------------
	
	
#===================================
if __name__ == '__main__':
	Run_Test()