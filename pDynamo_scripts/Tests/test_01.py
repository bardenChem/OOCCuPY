#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamo_Scripts import Scripts
import SimulationSystem 
import os, sys
#===================================
def info():
	print_message = "OOCCuPy pDynamo_Scripts Libray test #01:\t "
	print_message+= "Testing loading of force field files, spherical pruning and atom fixing.\n"

	print(print_message)
#-----------------------------------
def Run_Test():
	info()
	_parameters = {
		"Input_Type":"geometry",
		"crd_file":os.path.join("data","cyclohexane_single_frame.xyz"),	
	}
	#test load xyz
	test_01 = Scripts("test_01")
	test_01.Set_System(_parameters)
	test_01.SaveSystem()
	#test load gromacs topology and coordinate files 
	_parameters["Input_Type"] = "gromacs"
	_parameters["crd_file"] = os.path.join("data","1atp_peptide.gro")
	_parameters["top_file"] = os.path.join("data","1atp_peptide.top")	
	test_02 = Scripts("test_01")
	test_02.Set_System(_parameters)
	test_02.SaveSystem()
	#test load amber force field topology and coordinate files 
	_parameters["Input_Type"] = "amber"
	_parameters["crd_file"] = os.path.join("data","7tim.crd")
	_parameters["top_file"] = os.path.join("data","7tim.top")
	test_03 = Scripts("test_01")
	test_03.Set_System(_parameters)
	test_03.SaveSystem()
	#test load pkl and test spherical pruning and fixed atoms

	_parameters["Input_Type"] = "protein"
	_parameters["pdb_file"]   = "data/1l2y.pdb"
	test_04 = Scripts("test_01")
	test_04.Set_System(_parameters)
	test_04.SaveSystem()
	
	#test load pkl from amber FF and test spherical pruning and fixed atoms
	_parameters = {
		"Input_Type":"pkl",
		"pkl_file":"test_01/7tim.pkl",
		"spherical_prune":"*:LIG.248:C02",
		"spherical_prune_radius":25.0,
		"set_fixed_atoms":"*:LIG.248:C02",
		"free_atoms_radius":20.0,
		"set_reaction_crd":1,
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		"atoms_rc2":["*:LIG.*:O06","*:HIE.94:HE2","*:HIE.94:NE2"],
		"mass_constraints":["yes"],
		"type_rc1":"distance"	
	}	
	test_05 = Scripts("test_01")
	test_05.Set_System(_parameters)
	test_05.SaveSystem("7tim_pruned_and_fix")

#===================================
if __name__ == '__main__':
	if ( sys.argv[1] ) == "-print":	info()
	else: 							Run_Test()