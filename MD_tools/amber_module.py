#/usr/bin/env python
# -*- coding: utf-8 -*-
"""AMBER force field preparation and QM/MM simulation module.

This module provides tools for preparing PDB structures for molecular dynamics
simulations using AMBER. It handles hydrogen addition, parameter generation,
energy minimization, and QM input file generation for quantum mechanical calculations.

Typical usage example:
    amber = amber_mod('structure.pdb', H_opt=True, lig='ligand.pdb')
    amber.tleap_call()
    amber.sandro()
    amber.run_sander()
    amber.get_pdb()
"""


import os 
from pdb_class import*
from gmx_module import*


class amber_mod:
	"""Prepares molecular structures for AMBER MD simulations.
	
	This class handles the full pipeline for preparing PDB structures including
	hydrogen optimization, parameter generation, minimization, and conversion
	to optimized coordinates.
	
	Attributes:
		name (str): Input PDB filename.
		lig (str): Ligand filename for parameter generation.
		nameAf (str): Intermediate PDB with hydrogens.
		clear (str): Cleaned PDB filename.
		H_opt (bool): Whether to optimize hydrogens.
		amb_reduce (str): Path to AMBER reduce executable.
		pdb4amber (str): Path to pdb4amber executable.
		tleap (str): Path to tleap executable.
		sander (str): Path to sander executable.
		ambpdb (str): Path to ambpdb executable.
	"""

	def __init__(self                                    ,
			     pdbfile                                 ,
			     H_opt=True                              ,
			     lig  ="none"                            ,
			     path="/home/igorchem/programs/amber18/bin"):
		"""Initialize AMBER preparation object and run initial hydrogen addition.
		
		Automatically runs the AMBER reduce tool to add hydrogens to the structure.
		
		Args:
			pdbfile (str): Path to input PDB file.
			H_opt (bool, optional): Optimize hydrogen positions. Defaults to True.
			lig (str, optional): Ligand filename for processing. Defaults to "none".
			path (str, optional): Base path to AMBER bin directory. 
				Defaults to "/home/igorchem/programs/amber18/bin".
		"""

		self.name = pdbfile
		self.lig  = lig
		self.nameAf = self.name[:-4] + "_h.pdb"
		self.clear = self.name[:-4] +"_clear_h.pdb"
		self.H_opt = H_opt		
		self.amb_reduce = path + "/reduce"
		self.pdb4amber = path + "/pdb4amber"
		self.tleap = path + "/tleap"
		self.sander = path + "/sander" 
		self.ambpdb = path + "/ambpdb"


		print("========== starting reduce ===========")

		print(self.amb_reduce +" " + self.name + " > " + self.nameAf)

		os.system(self.amb_reduce +" " + self.name + " > " + self.nameAf)
		
		print("========== end of reduce =========")

		#print(self.pdb4amber +" -i " + self.nameAf +" -o " + self.clear)
		
		#os.system(self.pdb4amber +" -i " + self.nameAf +" -o " + self.clear)

	
	def tleap_call(self):
		"""Generate AMBER topology and coordinate files using tleap.
		
		Creates tleap input script with FF14SB force field parameters and TIP3P
		solvation box, then runs tleap to generate prmtop and inpcrd files.
		"""

		tleap_in  = "tleap \n" 
		tleap_in += "source leaprc.protein.ff14SB \n"
		tleap_in += "prot = loadPdb " + self.nameAf + "\n"
		tleap_in += "source leaprc.water.tip3p \n"
		tleap_in += "solvatebox prot TIP3PBOX 10.0 \n"
		tleap_in += "saveamberparm prot prmtop inpcrd\n"
		tleap_in += "quit"

		tleap_file = open('tleap_in','w')
		tleap_file.write(tleap_in)
		tleap_file.close()

		os.system(self.tleap + " -f tleap_in" )

	def antechamber(self,charge=0):
		"""Generate AMBER parameters for ligand/cofactor using antechamber.
		
		Uses antechamber with BCC charge assignment and parmchk for force field
		parameters. Generates mol2 and frcmod files for ligand.
		
		Args:
			charge (int, optional): Formal charge on the ligand. Defaults to 0.
		"""
		string  = ('antechamber -i %s -fi %s -o %s -fo mol2 -c bcc -nc %i' %('pdb/'+self.lig, self.lig[-3:], 'mol2/'+self.lig[:-3]+'mol2', charge))
		os.system(string)    
		string  = ('parmchk -i %s -f mol2 -o %s' %('mol2/'+self.lig[:-3]+'mol2', 'mol2/'+self.lig[:-3]+'frcmod'))
		os.system(string)
        
        
	def sandro(self):
		"""Create AMBER sander input file for hydrogen optimization minimization.
		
		Generates minimization input with restraints on heavy atoms (restrained
		hydrogen optimization), 5000 maximum cycles with 500 initial steepest descent.
		"""

		if self.H_opt==True:
			mini_in =  "minimize hydrogens \n"
			mini_in += " &cntrl \n"
			mini_in += "  imin=1,maxcyc=5000,ncyc=500, \n"
			mini_in += "  ntb=1,ntmin=1,ntpr=100,ntr=1, \n"
			mini_in += "  restraintmask="+'"!@H="'+",igb=0, \n"
			mini_in += "  restraint_wt=500.0,cut=10.0 \n"
			mini_in += "/"

		mini_file = open("minimize.in",'w')
		mini_file.write(mini_in)
		mini_file.close()

	def run_sander(self):
		"""Execute AMBER sander with the prepared minimization input.
		
		Runs geometry optimization using prmtop and inpcrd files, produces
		minimized structure in minimize.rst.
		"""

		text_to_run = self.sander +" -O -i minimize.in -o minimize.out -p prmtop -c inpcrd -ref inpcrd -r minimize.rst -inf minimize.mdinfo"
		print(text_to_run)

		os.system(text_to_run)		
	
	def get_pdb(self):
		"""Convert optimized AMBER restraint file back to PDB format.
		
		Uses ambpdb to extract coordinates from minimize.rst, writes optimized
		PDB, then parses and re-writes with pdb_class for consistency.
		"""
	
		text_to_run = self.ambpdb +" -c inpcrd < minimize.rst > "+ self.name[:-4] +"_opt.pdb"
		
		os.system(text_to_run)

		pdb_opt = protein(self.name,amber=True)
		pdb_opt.pdb_parse(self.name[:-4]+"_opt.pdb")
		pdb_opt.write_pdb(self.name[:-4]+"_opt_final.pdb")

	def mopac_inp(self): 
		"""Generate MOPAC quantum chemistry input from optimized structure.
		
		Reads the optimized PDB, converts coordinates to XYZ format,
		creates MOPAC input file for QM calculations.
		"""

		opt_pdb = protein(self.name[:-4]+"_opt_final.pdb")
		opt_pdb.pdb_parse()
		opt_pdb.mopac_mop()
	
	def equilibration(self):
		"""Placeholder for equilibration protocol setup."""
		pass

	def production(self):
		"""Placeholder for production MD protocol setup."""
		pass

	def analysis_MD(self):
		"""Placeholder for MD trajectory analysis."""
		pass
'''	
a = amber_mod("1a2y_p1.pdb")
a.tleap_call()
a.sandro()
'''
