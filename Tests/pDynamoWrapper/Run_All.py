#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
List of Tests 

test_01_system_setup.py -- Set/save systems with Molecular Mechanics Force Field; Pruning and atoms fixing.
			requires files from "data" folder: 1atp_peptide.gro,1atp_peptide.top,7tim.crd,7tim.top,1l2y.pdb.
test_02_quantum_methods.py -- Set/save systems with Quantum Mechanics methods.
			requires files from "data" folder: cyclohexane_single_frame.xyz
test_03_geometry_optimization.py -- Geometry Optimizations algoritims and Analysis.
			requires files from tests 01 and 02.
test_04_molecular_dynamics.py -- Molecular Dynamics Algoritims and Analysis (test parallel version, AKA sampling format).
			requires files from test 03.
test_05_relaxed_scan_1d.py -- Unidimensional Relaxed Scans.
			requires files from "data" folder: 7tim.crd,7tim.top.
test_06_relaxed_scan_2d.py -- Bidimensional Relaxed Scans.
			requires files from test 05.
test_07_qmmm_molecular_dynamics.py -- QCMM molecular dynamics
test_08_restricted_qmmm_molecular_dynamics.py -- QCMM restricted molecular dynamics
test_09_energy_refinement_1d.py -- Energy Refinement
test_10_energy_refinement_variable_qm_region.py -- Energy Refinement with changing quantum region
test_11_umbrella_sampling_1d.py -- Unidimensional Umbrella Sampling + WHAM with and without step optimization
test_12_umbrella_sampling_2d.py -- Bidimensional Umbrella Sampling + WHAM
test_13_nudged_elastic_band.py -- Reaction Path Algorithms
test_14_normal_modes.py -- Normal modes


#extra
test_xx.py -- Monte Carlo Simulation (pDyamo examples)
test_xx.py -- Simulated Annealing (pDyamo examples)
test_xx.py -- Set/save protein system
test_xx.py -- Modeling protein system ( pDynamo examples )
test_xx.py -- Dihedral Relaxed Scans. 

'''
import subprocess
import sys
from pathlib import Path

test_list = [
	"test_01_system_setup.py",
	"test_02_quantum_methods.py",
	"test_03_geometry_optimization.py",
	"test_04_molecular_dynamics.py",
	"test_05_relaxed_scan_1d.py",
	"test_06_relaxed_scan_2d.py",
	"test_07_qmmm_molecular_dynamics.py",
	"test_08_restricted_qmmm_molecular_dynamics.py",
	"test_09_energy_refinement_1d.py",
	"test_10_energy_refinement_variable_qm_region.py",
	"test_11_umbrella_sampling_1d.py",
	"test_12_umbrella_sampling_2d.py",
	"test_13_nudged_elastic_band.py",
	"test_14_normal_modes.py",
	"test_15_external_energy_refinement_2d.py",
	"test_16_reaction_path_analysis.py",
	"test_17_split_trajectory.py",
	"test_18_qmmm_setup.py",
	"test_19_orca_quantum_methods.py",
	"test_20_pyscf_quantum_methods.py",
	"test_21_qmmm_trajectory_analysis.py",
]


#----------------------------------------------------------------
def RUN_ALL_TESTS():

	base_dir = Path(__file__).resolve().parent
	for test in test_list:
		subprocess.run([sys.executable, str(base_dir / test)], check=True)

#=================================================================
if __name__ == '__main__':
	RUN_ALL_TESTS()
