#/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamo_Scripts import Scripts
import SimulationSystem 
import os
#===================================
def info():
	print_message = "OOCCuPy pDynamo_Scripts Libray test #02:\t "
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
	_parameters["pkl_file"]         = "test_05/qcmm_optam1/7tim_am1_opt_PF.pkl"

	_path   = "test_05/Multiple_Distance_rm1/ScanTraj.ptGeo"

	simulation_parameters = { "xnbins":20			    ,
				   "source_folder":_path                , 
				   "folder":"test_02_pyscf"             ,
				   "QCcharge":-3		                    ,
				   "multiplicity":1 	                ,
				   "basis":"6-31G*"						,
				   "functional":"b3lyp"					,
				   "pySCF_method":"RKS"					,
				   "simulation_type":"Energy_Refinement",
				   "Software":"pySCF"					}	

	test_03 = Scripts("test_02_pyscf")
	test_03.Set_System(_parameters)
	test_03.Run_Simulation(simulation_parameters)
	


#===================================
if __name__ == '__main__': 
	if ( sys.argv[1] ) == "-print":	info()
	else: 						    Run_Test()