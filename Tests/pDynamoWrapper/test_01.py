#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamoWrapper.pDynamoWrapper import Wrapper
from pDynamoWrapper.SimulationSystem import SimulationSystem
import os, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from config import get_config


#===================================
def info():
	print_message = "OOCCuPy pDynamoWrapper Libray test #01:\t "
	print_message+= "Testing loading of force field files, spherical pruning and atom fixing.\n"

	print(print_message)
#-----------------------------------
def Run_Test():
	info()

	# Get OOCCuPY root from config
	config = get_config()
	ooccupy_root = config.get_ooccupy_root()

	# Build paths using ooccupy_root
	folder = ooccupy_root / "Tests" / "pDynamoWrapper" / "test_01"
	data_dir = ooccupy_root / "data" / "pDynamoWrapper" / "data"

	_parameters = {
		"Input_Type":"geometry",
		"crd_file":os.path.join(ooccupy_root, "data", "cyclohexane_single_frame.xyz"),
	}
	#test load xyz
	test_01 = Wrapper(folder)
	test_01.Set_System(_parameters)
	test_01.SaveSystem()
	#test load gromacs topology and coordinate files 
	_parameters["Input_Type"] = "gromacs"
	_parameters["crd_file"] = os.path.join(ooccupy_root, "data", "1atp_peptide.gro")
	_parameters["top_file"] = os.path.join(ooccupy_root, "data", "1atp_peptide.top")	
	test_02 = Wrapper(folder)
	test_02.Set_System(_parameters)
	test_02.SaveSystem()
	#test load amber force field topology and coordinate files 
	_parameters["Input_Type"] = "amber"
	_parameters["crd_file"] = os.path.join(ooccupy_root, "data", "7tim.crd")
	_parameters["top_file"] = os.path.join(ooccupy_root, "data", "7tim.top")
	test_03 = Wrapper(folder)
	test_03.Set_System(_parameters)
	test_03.SaveSystem()
	#test load pkl and test spherical pruning and fixed atoms

	_parameters["Input_Type"] = "protein"
	_parameters["pdb_file"]   = os.path.join(ooccupy_root, "data", "1l2y.pdb")
	test_04 = Wrapper(folder)
	test_04.Set_System(_parameters)
	test_04.SaveSystem()
	
	#test load pkl from amber FF and test spherical pruning and fixed atoms
	_parameters = {
		"Input_Type":"pkl",
		"pkl_file": os.path.join(folder, "7tim.pkl"),
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
	test_05 = Wrapper(folder)
	test_05.Set_System(_parameters)
	test_05.SaveSystem( "7tim_pruned_and_fix" )

#===================================
if __name__ == '__main__':
	Run_Test()