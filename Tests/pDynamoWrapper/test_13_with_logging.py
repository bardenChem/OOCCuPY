#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ENHANCED test_13.py WITH DEBUG LOGGING
========================================

This is an enhanced version of test_13.py that demonstrates how to enable
comprehensive logging for debugging the NEB with trajectory source issue.

The logging will help identify where the segmentation fault occurs on the
remote machine by capturing detailed execution information.

Usage:
------
python test_13_with_logging.py

Output:
-------
- Console: Real-time status messages
- test_13/GeometrySearcher.log: Main log file
- test_13/GeometrySearcher_DEBUG.log: Detailed debug file (verbose)
- test_13_trj/GeometrySearcher.log: Main log file for trajectory test
- test_13_trj/GeometrySearcher_DEBUG.log: Detailed debug file for trajectory test
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pDynamoWrapper import Wrapper
import os, sys

from config import get_config
config = get_config()
ooccupy_root = config.get_ooccupy_root()

folder = os.path.join(ooccupy_root, "Tests", "pDynamoWrapper", "test_13")
folderTRJ = os.path.join(ooccupy_root, "Tests", "pDynamoWrapper", "test_13_trj")
folder05 = os.path.join(ooccupy_root, "Tests", "pDynamoWrapper", "test_05")


#====================================================
def info():
	print_message =  "OOCCuPy pDynamoWrapper Library test #13 (ENHANCED WITH LOGGING):\n\t "
	print_message += "Testing the setting and run of Nudged Elastic Band routines.\n"
	print_message += "\nWith comprehensive logging to help debug segmentation fault issues.\n"
	print(print_message)


#----------------------------------
def Run_Test_WithLogging():
	'''
	Enhanced test with full logging enabled.
	This helps debug the NEB with trajectory source issue.
	'''
	info()

	_path   = os.path.join( folder05, "Multiple_Distance_am1","ScanTraj.ptGeo" )
	init_path    = os.path.join( _path, "frame0.pkl")
	final_path   = os.path.join( _path, "frame19.pkl")
	saddle_coord = os.path.join( _path, "frame12.pkl") 

	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join(folder05,"qcmm_optam1","7tim_am1_opt_PF.pkl"),				
	}

	# ===== TEST 1: NEB from coordinates =====
	print("\n" + "="*70)
	print("TEST 1: NEB from Initial and Final Coordinates")
	print("="*70)
	
	parameters_NEB = {  
		"init_coord":init_path,
		"final_coord":final_path,
		"traj_bins":12,
		"refine_methods":["rm1","am1","pm3","pddgpm3","pm6"],
		"RMS_growing_intial_string":1.0,
		"simulation_type":"NEB",
		"spring_force_constant":800.0,
		"rmsGradient":0.10,
		"fixed_terminal_images":"no",
		# LOGGING PARAMETERS - Enable detailed logging
		"enable_debug_file": True,  # Creates GeometrySearcher_DEBUG.log
		"debug_verbosity": "DEBUG",  # Set to DEBUG, INFO, WARNING, ERROR
	}

	test_01 = Wrapper(folder)
	test_01.Set_System(system_parameters)
	
	print("\n>>> Starting NEB from coordinates...")
	print("    Check: {}/GeometrySearcher.log".format(folder))
	print("    Check: {}/GeometrySearcher_DEBUG.log (detailed)".format(folder))
	
	test_01.Run_Simulation(parameters_NEB)
	test_01.SaveSystem()
	
	print("\n✓ TEST 1 COMPLETED")
	print("  Output logs created in: {}".format(folder))

	# ===== TEST 2: NEB from trajectory (trajectory source) =====
	# THIS IS WHERE THE SEGMENTATION FAULT LIKELY OCCURS ON REMOTE MACHINE
	print("\n" + "="*70)
	print("TEST 2: NEB from Trajectory Source")
	print("="*70)
	print("\n!!! THIS TEST USES traj_source WHICH CAUSES SEGFAULT ON REMOTE MACHINE !!!")
	print("    Enable logging to diagnose the issue.\n")

	parameters_NEBT = { 
		"traj_source":_path,
		"refine_methods":["rm1"],
		"RMS_growing_intial_string":1.0,
		"simulation_type":"NEB",
		"spring_force_constant":800.0,
		"rmsGradient":0.1,
		"fixed_terminal_images":"no",
		# LOGGING PARAMETERS - CRITICAL for debugging this failing case
		"enable_debug_file": True,   # Creates GeometrySearcher_DEBUG.log with full details
		"debug_verbosity": "DEBUG",  # Maximum verbosity to capture all details
	}
	
	test_02 = Wrapper(folderTRJ)
	test_02.Set_System(system_parameters)
	
	print("\n>>> Starting NEB from trajectory source...")
	print("    Source trajectory: {}".format(_path))
	print("    Check: {}/GeometrySearcher.log".format(folderTRJ))
	print("    Check: {}/GeometrySearcher_DEBUG.log (detailed)".format(folderTRJ))
	
	print("\n>>> DETAILED DEBUG NOTES <<<")
	print("If this crashes, check GeometrySearcher_DEBUG.log for:")
	print("  1. Last successful checkpoint before crash")
	print("  2. Trajectory loading confirmation")
	print("  3. Frame count verification")
	print("  4. NEB optimization progress")
	print("  5. Finalize method execution")
	
	try:
		test_02.Run_Simulation(parameters_NEBT)
		test_02.SaveSystem()
		print("\n✓ TEST 2 COMPLETED")
		print("  Output logs created in: {}".format(folderTRJ))
	except Exception as e:
		print(f"\n✗ TEST 2 FAILED WITH EXCEPTION: {type(e).__name__}")
		print(f"  Message: {str(e)}")
		print(f"\n  CRITICAL: Check {folderTRJ}/GeometrySearcher_DEBUG.log")
		print(f"  This file contains the full execution trace leading to the crash.")
		raise

	# ===== ANALYSIS =====
	print("\n" + "="*70)
	print("LOG ANALYSIS RECOMMENDATIONS")
	print("="*70)
	print("""
To analyze these logs for the segmentation fault issue:

1. Compare both log files:
   diff {}/GeometrySearcher.log {}/GeometrySearcher.log

2. Look for error checkpoints in TEST 2 debug log:
   grep "ERROR\\|CHECKPOINT" {}/GeometrySearcher_DEBUG.log

3. Find the last successful operation:
   grep "CHECKPOINT.*OK" {}/GeometrySearcher_DEBUG.log | tail -5

4. Check trajectory loading details:
   grep "Trajectory\\|Frame" {}/GeometrySearcher_DEBUG.log

5. On the remote machine, if crash occurs, send:
   - {}/GeometrySearcher_DEBUG.log
   - System info: uname -a, free -h, etc.
   This will help identify machine-specific issues.
	""".format(folder, folderTRJ, folderTRJ, folderTRJ, folderTRJ, folderTRJ))


if __name__ == "__main__":
	Run_Test_WithLogging()
	print("\n" + "="*70)
	print("ALL TESTS COMPLETED")
	print("="*70)
	print("\nLog files created in:")
	print(f"  - {folder}/GeometrySearcher*.log")
	print(f"  - {folderTRJ}/GeometrySearcher*.log")
	print("\nFor remote debugging, copy these files and compare with local results.")
