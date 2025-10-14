#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamoWrapper import Wrapper
import SimulationSystem 
import os
#===================================
def Run_Test():
	'''
	Test the setting of all quantum methods classes available.
	Internal semi-empirical methods in pDynamo
	DFTB+ 
	Orca 
	Mopac
	pySCF 
	'''
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
		test_01 = Wrapper("test_02")
		test_01.Set_System(_parameters)
		test_01.SaveSystem("cyclohexane"+smo+".pkl")
	
	#test QC/MM from gromacs
	_parameters["Input_Type"]      ="pkl"
	_parameters["pkl_file"]        ="test_01/1atp_peptide.pkl"
	_parameters["set_qc_region"]   ="yes"
	_parameters["residue_patterns"]=["*:ARG.19:*"]
	_parameters["QCcharge"]        = 1

	test_02 = Wrapper("test_02")
	test_02.Set_System(_parameters)
	test_02.SaveSystem("1atp_peptide_qmmm")

	#test QC/MM from AMBER
	_parameters["residue_patterns"] = ["*:LIG.248:*","*:GLU.164:*","*:HIE.94:*"]
	_parameters["pkl_file"]         = "test_01/7tim.pkl"
	
	test_03 = Wrapper("test_02")
	test_03.Set_System(_parameters)
	test_03.SaveSystem("7tim_qcmm_am1")

	_parameters["pkl_file"]        = "test_01/7tim_pruned_and_fix.pkl"
	test_03 = Wrapper("test_02")
	_parameters["Hamiltonian"]="rm1"
	test_03.Set_System(_parameters)
	test_03.SaveSystem("7tim_qcmm_rm1_pruned")

	#test QC/MM com AMBER carregado de pkl


	
#===================================
if __name__ == '__main__': Run_Test()