#/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamoWrapper import Wrapper
import SimulationSystem 
import os
#===================================
def info():
	print_message = "OOCCuPy pDynamoWrapper Libray test #20:\t "
	print_message += "Test the setting of quantum methods with pySCF .\n"

	print(print_message)
def Run_Test():
	'''
	'''
	_parameters = {
		"Input_Type":"pkl",
		"set_energy_model":"QM",
		"functional":"b3lyp",
		"basis":"6-31G*",
		"QCcharge":-3,
		"set_reaction_crd":1,
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		"multiplicity":1,
		"mass_constraint":"True",
		"type":"Distance",
		"method_class":"pySCF"
	}
	_parameters["set_qc_region"]    = "yes"
	_parameters["residue_patterns"] = ["*:LIG.248:*","*:GLU.164:*","*:HIE.94:*","*:ASN.9:*"]
	_parameters["pkl_file"]         = "Test/pDynamoWrapper/test_05/qcmm_optam1/7tim_am1_opt_PF.pkl"

	_path   = "Test/pDynamoWrapper/test_05/Multiple_Distance_rm1/ScanTraj.ptGeo"
	if not os.path.exists(_path): os.system("python3 test_05")

	simulation_parameters = { "xnbins":20			    ,
				   "source_folder":_path                , 
				   "folder":"test_20"                   ,
				   "QCcharge":-3		                ,
				   "multiplicity":1 	                ,
				   "basis":"6-31G*"						,
				   "functional":"b3lyp"					,
				   "pySCF_method":"RKS"					,
				   "simulation_type":"Energy_Refinement",
				   "Software":"pySCF"					}	

	test_03 = Wrapper("test_20")
	test_03.Set_System(_parameters)
	test_03.Run_Simulation(simulation_parameters)
	
#===================================
if __name__ == '__main__': 
	Run_Test()