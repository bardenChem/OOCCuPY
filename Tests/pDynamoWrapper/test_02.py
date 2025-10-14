#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamoWrapper import Wrapper
import SimulationSystem 
import os,sys
#===================================
def info():
	print_message = "OOCCuPy pDynamoWrapper Libray test #02:\t "
	print_message += "Test the setting of all quantum methods classes available.\n"

	print(print_message)
#------------------------------------
def Run_Test():
	'''
	Test the setting of all quantum methods classes available.
	Internal semi-empirical methods in pDynamo
	DFTB+ 
	Orca 
	Mopac
	pySCF 
	'''
	info()
	SMOmodels = ["am1","am1dphot","pddgpm3","pm3","rm1","pm6"]		
	_parameters = {
		"Input_Type":"geometry",
		"crd_file":os.path.join("data","cyclohexane_single_frame.xyz"),
		"set_energy_model":"QM",
		"Hamiltonian":"am1",
		"method_class":"SMO"
	}
	for smo in SMOmodels:
		_parameters["Hamiltonian"] = smo
		test_01 = Wrapper()
		test_01.Set_System(_parameters)
	#------------------------------------	
	_parameters["method_class"]="ORCA"
	_parameters["functional"]  ="HF"
	_parameters["basis"]       ="6-31G*"
	test_03 = Wrapper()
	test_03.Set_System(_parameters)	
	input()
	#--------------------------------
	_parameters["method_class"]="pySCF"
	_parameters["functional"]  ="b3lyp"
	_parameters["basis"]       ="6-31G*"
	test_04 = Wrapper()
	test_04.Set_System(_parameters)
	input()
	#----------------------------------
	_parameters["method_class"]="abinitio"
	_parameters["basis"]       ="dgauss-dzvp"
	test_05 = Wrapper()
	test_05.Set_System(_parameters)
	#----------------------------------
	
#===================================
if __name__ == '__main__': 
	if ( sys.argv[1] ) == "-print":	info()
	else: 						    Run_Test()