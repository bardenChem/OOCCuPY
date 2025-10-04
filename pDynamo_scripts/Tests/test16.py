#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamo_Scripts import Scripts
import SimulationSystem 
import os, sys
#==============================================
def Run_Test():


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

	_path   = "test_06/Multiple_Distance_rm1/ScanTraj.ptGeo"

	analysis_parameters = {
		"analysis_type":"Energy_Plots",
		"log_name":"test_06/Multiple_Distance_rm1/ScanTraj.log",
		"retrieve_path":_path,
		"xsize":12,
		"ysize":12,
		"contour_lines":14,
		"folder":"test_16",
		"type":"2D",
		"xlim_list": [-1.4,-0.4]   ,
		"ylim_list": [-1.2,-0.2]   ,
		"in_point" :[1,1]   ,
		"fin_point":[10,10]    
	}

	test_01 = Scripts("test_16")
	test_01.Set_System(system_parameters)
	test_01.Run_Analysis(analysis_parameters)
	test_01.SaveSystem()

#===================================
if __name__ == '__main__': Run_Test()