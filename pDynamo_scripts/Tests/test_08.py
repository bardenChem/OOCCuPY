#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamo_Scripts import Scripts
import SimulationSystem 
import os, sys
#===================================
def Run_Test():
	'''
	Test molecular dynamics algorithms with qmmm
	'''

	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_05","qcmm_optam1","7tim_am1_opt_PF.pkl"),
		"set_reaction_crd":2,	
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		"atoms_rc2":["*:LIG.*:O06","*:HIE.94:HE2","*:HIE.94:NE2"],
		"type_rc1":"Distance",
		"type_rc2":"Distance",
		"mass_constraints":["yes","yes"],
	}

	simulation_parameters = {
				  "temperature": 315.15,
				  "simulation_type":"Restricted_Molecular_Dynamics",
				  "equilibration_nsteps":5000,
				  "production_nsteps":10000,
				  "heating_nsteps":2000,
				  "sampling_equilibration":100,
				  "sampling_production":50,
				  "sampling_heating":50,
				  "force_constants":[300.0,300.0],
				  "log_frequency":10
				}
	
	#------------------------------------
	#protocol production
	test_01 = Scripts("test_08")
	test_01.Set_System(system_parameters)
	simulation_parameters["trajectory_name"]="7timQCMD_restricted"
	test_01.Run_Simulation(simulation_parameters)
	test_01.SaveSystem()
	#-----------------------------------
	
	
#===================================
if __name__ == '__main__': Run_Test()