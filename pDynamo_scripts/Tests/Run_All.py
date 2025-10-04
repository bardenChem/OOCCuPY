#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
List of Tests 

test_01.py -- Set/save systems with Molecular Mechanics Force Field; Pruning and atoms fixing.
			requires files from "data" folder: 1atp_peptide.gro,1atp_peptide.top,7tim.crd,7tim.top,1l2y.pdb.
test_02.py -- Set/save systems with Quantum Mechanics methods.
			requires files from "data" folder: cyclohexane_single_frame.xyz
test_03.py -- Geometry Optimizations algoritims and Analysis. 
			requires files from tests 01 and 02.
test_04.py -- Molecular Dynamics Algoritims and Analysis (test parallel version, AKA sampling format).
			requires files from test 03.
test_05.py -- Unidimensional Relaxed Scans.
			requires files from "data" folder: 7tim.crd,7tim.top.
test_06.py -- Bidimensional Relaxed Scans. 
			requires files from test 05.
test_07.py -- QCMM molecular dynamics
test_08.py -- QCMM restricted molecular dynamics
test_09.py -- Energy Refinement
test_10.py -- Energy Refinement with changing quantum region 
test_11.py -- Unidimensional Umbrella Sampling + WHAM with and without step optimization
test_12.py -- Bidimensional Umbrella Sampling + WHAM 
test_13.py -- Reaction Path Algorithms


#extra
test_xx.py -- Monte Carlo Simulation (pDyamo examples)
test_xx.py -- Simulated Annealing (pDyamo examples)
test_xx.py -- Set/save protein system
test_xx.py -- Modeling protein system ( pDynamo examples )
test_xx.py -- Dihedral Relaxed Scans. 

'''


test_list = ["test_01.py",
			 "test_02.py",
			 
