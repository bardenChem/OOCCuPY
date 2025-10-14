#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamoWrapper import Scripts
import SimulationSystem 
import os, sys
#==============================================
def info():
	print_message =  "OOCCuPy pDynamoWrapper Libray test #05:\t "
	print_message += "Testing the setting and run of relaxed unidimensional scans.\n"

	print(print_message)
#------------------------------------------------
def Prepare_MM_System():
	'''
	'''
	_opt_pars = {
		"Input_Type":"amber",
		"simulation_type":"Geometry_Optimization",
		"save_format":".dcd",
		"save_frequency":20,
		"rmsGradient":2.0,
		"maxIterations":2200,
		"folder":"test_03",
		"trajectory_name":"opt_full_tim.ptGeo"
	}

	#preparation of 7tim structure	
	#test load amber force field topology and coordinate files 
	_opt_pars["crd_file"] = os.path.join("data","7tim.crd")
	_opt_pars["top_file"] = os.path.join("data","7tim.top")
	test_03 = Scripts("test_05")
	test_03.Set_System(_opt_pars)
	test_03.Run_Simulation(_opt_pars)
	test_03.SaveSystem()

#-----------------------------------------------
def Prepare_Prune_System():
	'''
	'''
	_prune_pars= {
		"Input_Type":"pkl",
		"pkl_file":"test_05/7tim.pkl",
		"simulation_type":"Geometry_Optimization",
		"spherical_prune":"*:LIG.248:C02",
		"spherical_prune_radius":25.0,
		"set_fixed_atoms":"*:LIG.248:C02",
		"free_atoms_radius":20.0,
		"save_format":".dcd",
		"save_frequency":20,
		"rmsGradient":1,
		"maxIterations":2200,
		"trajectory_name":"opt_pruned_tim.ptGeo"		
	}

	test_03_b = Scripts("test_05")
	test_03_b.Set_System(_prune_pars)
	test_03_b.Run_Simulation(_prune_pars)
	test_03_b.SaveSystem("7tim_optMM")


#-----------------------------------------------
def Set_QC_MM(_hamiltonian="am1"):
	'''
	'''
	_qc_mmpars = {
		"Input_Type":"pkl",
		"pkl_file":"test_05/7tim_optMM.pkl",
		"set_energy_model":"QM",
		"Hamiltonian":_hamiltonian,
		"method_class":"SMO",
		"set_qc_region":"yes",
		"residue_patterns":["*:LIG.248:*","*:GLU.164:*","*:HIE.94:*","*:ASN.9:*"],
		"QCcharge":-3,
		"save_format":".dcd",
		"save_frequency":20,
		"rmsGradient":0.1,
		"maxIterations":2200,
		"trajectory_name":"opt_qcmm_tim_"+_hamiltonian+".ptGeo",
		"simulation_type":"Geometry_Optimization"
	}

	test_03_c = Scripts("test_05/qcmm_opt"+_hamiltonian)
	test_03_c.Set_System(_qc_mmpars)
	test_03_c.Run_Simulation(_qc_mmpars)
	test_03_c.SaveSystem("7tim_"+_hamiltonian+"_opt_PF")

def Set_QC_MM_multiple(methods):
	'''
	'''
	SMOmodels = methods	

	for smo in SMOmodels:
		Set_QC_MM(smo)

#-------------------------------------------------
def Simple_Distance(_hamiltonian):
	'''
	'''

	if not os.path.exists( os.path.join("test_05","qcmm_optm"+_hamiltonian+"_opt_PF.pkl") ):
		Set_QC_MM(_hamiltonian)

	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_05","qcmm_opt"+_hamiltonian,"7tim_"+_hamiltonian+"_opt_PF.pkl"),		
		"set_reaction_crd":1,
		"atoms_rc1":["*:LIG.*:H02","*:GLU.164:OE2"],
		"type_rc1":"Distance",
		"mass_constraints":["no"]
	}		
	scan1_parameters = {
		"simulation_type":"Relaxed_Surface_Scan",
		"dincre_rc1":-0.1,
		"nsteps_rc1":12,
		"maxIterations":2200,
		"optmizer":"SteepestDescent",
		"force_constants":[1200.0,1200.0]
	}
	#test simple distance
	test_01 = Scripts("test_05/Simple_Distance_"+_hamiltonian)
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(scan1_parameters)
	test_01.SaveSystem("Simple_DistanceScan")


#-------------------------------------------------
def Multiple_Distance(_hamiltonian):
	'''
	'''

	if not os.path.exists( os.path.join("test_05","qcmm_optm"+_hamiltonian+"_opt_PF.pkl") ):
		Set_QC_MM(_hamiltonian)

	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_05","qcmm_opt"+_hamiltonian,"7tim_"+_hamiltonian+"_opt_PF.pkl"),		
		"set_reaction_crd":1,
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		"type_rc1":"Distance",
		"maxIterations":2200,
		"mass_constraints":["yes"],
	}		
	scan1_parameters = {
		"simulation_type":"Relaxed_Surface_Scan",
		"dincre_rc1":0.1,
		"optmizer":"SteepestDescent",
		"maxIterations":2200,
		"nsteps_rc1":20,
		"log_frequency":10,
		"force_constants":[1200.0,1200.0]
	}
	#test simple distance
	test_01 = Scripts("test_05/Multiple_Distance_"+_hamiltonian)
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(scan1_parameters)
	test_01.SaveSystem("Multiple_DistanceScan")

#-----------------------------------------------
def Run_Test():
	info()
	if not os.path.exists( os.path.join("test_05","7tim.pkl") ):
		Prepare_MM_System()
	
	if not os.path.exists( os.path.join("test_05","7tim_optMM.pkl") ):
		Prepare_Prune_System()
	
	SMOmodels = ["am1","rm1","pm3","pm6","am1dphot","pddgpm3"]
	for smo in SMOmodels:		
		Simple_Distance(smo)		
		Multiple_Distance(smo)
	
#===================================
if __name__ == '__main__':
	if ( sys.argv[1] ) == "-print":	info()
	else: 						    Run_Test()