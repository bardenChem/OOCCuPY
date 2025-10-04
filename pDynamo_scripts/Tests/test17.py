#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from pDynamo_Scripts import Scripts
import SimulationSystem 
import os, sys
#====================================================
def Run_Test():
	'''
	'''
	_path = "test_05/Multiple_Distance_rm1/ScanTraj.ptGeo"

	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("test_05","qcmm_optam1","7tim_am1_opt_PF.pkl"),		
		"set_reaction_crd":2,	
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		"atoms_rc2":["*:LIG.*:O06","*:HIE.94:HE2","*:HIE.94:NE2"],
		"type_rc1":"Distance",
		"type_rc2":"Distance",
		"mass_constraints":["yes","yes"],
		"analysis_type":"Split_Traj",
		"break_point":5,
		"trajectory_name":_path
	}

	test_02 = Scripts("test_17")
	test_02.Set_System(system_parameters)
	test_02.Run_Analysis(system_parameters)
	test_02.SaveSystem()

#===================================
if __name__ == '__main__': Run_Test()