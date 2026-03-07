#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""MOPAC quantum chemistry input file generation module.

This module provides tools for generating MOPAC input files for semi-empirical
quantum mechanical calculations. Supports various charge states and calculation
methods (PM7, PM6, etc.). Can process protein structures and generate multiple
input files for neutral, cation, and anion states.

Typical usage example:
    mopac = mopac_inp('structure.xyz', charge=0, multi=1, 
                      inpnam='output.mop', mozyme=True, method='PM7')
    mopac.write_mop()
"""

#=======================================================================

#load modules

import os
from xyz_class import*
from pdb_class import*

#======================================================================

class mopac_inp:
	"""Generator for MOPAC semi-empirical QM input files.
	
	Reads molecular structures and creates MOPAC input files with specified
	quantum mechanical parameters, charge state, and solvation model.
	Supports MOZYME approximation for large systems.
	
	Attributes:
		name (str): Input XYZ filename.
		inpnam (str): Output MOPAC input filename.
		charge (str): Charge specification for MOPAC.
		multi (int): Multiplicity (1=singlet, 2=doublet, etc.).
		mult (str): Multiplicity label for MOPAC format.
		solvent (bool): Include implicit solvation (eps=78.4).
		hamilt (str): QM method/hamiltonian (e.g., 'PM7', 'PM6').
		xyz (xyz_parser): Parser object with molecular geometry.
		mozyme (str): MOZYME keyword for large system approximation.
		mgf (str): Graph format output flag.
	"""

	def __init__(self  ,
				xyzfile,
				charge ,
				multi  ,
				inpnam ,
				mozyme , 
				mgf    , 
				method):
		"""Initialize MOPAC input generator from XYZ structure file.
		
		Args:
			xyzfile (str): Path to input XYZ coordinate file.
			charge (int): Formal charge on the system.
			multi (int): Spin multiplicity (1=singlet, 2=doublet).
			inpnam (str): Output MOPAC filename.
			mozyme (bool): Enable MOZYME large-system approximation.
			mgf (bool): Enable molecular graphics format output.
			method (str): Quantum mechanical method (e.g., 'PM7', 'PM6', 'RM1').
		"""

		self.name    = xyzfile
		self.inpnam  = inpnam
		self.charge  = ""
		self.multi   = multi
		self.mult    = "Singlet"
		self.solvent = True
		self.hamilt  = method
		self.xyz     = None
		self.mozyme  = ""
		self.mgf     = ""

		if not self.multi == 1:
			self.mult = "Doublet"

		self.charge = "charge=" +str(charge)
		if mozyme:
			self.mozyme = "mozyme"
			self.charge = ""
		self.xyz = xyz_parser(self.name)
		self.xyz.parse_xyz()
		if mgf:
			self.mgf = "graphf"
			
	def write_mop(self):
		"""Write MOPAC input file with atomic coordinates and parameters.
		
		Generates a MOPAC-formatted input file with method specification,
		charge, solvation, and atomic coordinates in proper format.
		"""

		mop_inp = open(self.inpnam,'w')
		mop_text = ''
		mop_text += '{0} 1SCF ALLVECS VECTOR aux {1} {2} {3} eps=78.4 large\n\n\n'.format(self.hamilt,self.mozyme,self.charge,self.mgf)

		for i in range(self.xyz.Natoms):
			mop_text +="{0}  {1}  1 {2}  1 {3} \n".format(self.xyz.AtomLabels[i],self.xyz.xCoord[i],self.xyz.yCoord[i],self.xyz.zCoord[i])

		mop_inp.write(mop_text)
		mop_inp.close()


def run_all(met):
	"""Batch process all PDB files and generate MOPAC inputs for multiple charge states.
	
	Finds all PDB files in current directory, converts them to XYZ format, and
	generates MOPAC input files for neutral, cation (charge=+1), and anion 
	(charge=-1) states. Creates shell script to run all MOPAC jobs.
	
	Args:
		met (str): Quantum mechanical method for all calculations (e.g., 'PM7').
	"""
	lists = glob.glob('*.pdb')
	for pdb in lists:
		a = protein(name=pdb,amber=True)
		a.pdb_parse(pdb)
		a.write_xyz()
	listxyz = glob.glob('*.xyz')
	for xyz in listxyz:
		a = mopac_inp(xyz,1,1,xyz[:-4]+"neutro.mop",met)
		a.write_mop()
		b = mopac_inp(xyz,2,2,xyz[:-4]+"cation.mop",met)
		b.write_mop()
		c = mopac_inp(xyz,0,2,xyz[:-4]+"anion.mop",met)
		c.write_mop()				
	script = open("run_mopac.sh",'w')
	script_text ="#!/bin/sh \n"
	listsmop = glob.glob('*.mop')
	for mop in listsmop:
		script_text += "/opt/mopac/MOPAC2016.exe {0}\n".format(mop)
	script.write(script_text)
	script.close()
