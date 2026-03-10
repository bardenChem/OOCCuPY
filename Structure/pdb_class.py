#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""PDB structure parsing and manipulation module.

This module provides classes for parsing and manipulating PDB (Protein Data Bank)
structure files. It includes functionality for reading atomic coordinates, defining
residues, calculating geometric properties, and performing various structural 
modifications such as pruning water molecules, removing ions, and splitting complexes.

Typical usage example:
    protein_obj = protein('structure.pdb')
    protein_obj.remove_waters()
    protein_obj.prune_ions()
    protein_obj.write_pdb('output.pdb')
"""

#=======================================================================

#load modules

import os
import math

#=======================================================================

# residues data dictionaries

RES_N0C ={'GLY':2,'ALA':3,'VAL':5,'PHE':9,'ILE':6,'LEU':6,'PRO':5,'MET':5,'ASP':4,'GLU':5,'LYS':6,'ARG':6,'SER':3,'THR':4,'TYR':9,'CYS':3,'ASN':4,'GLN':5,'HIS':6,'TRP':11}
RES_N0H ={'GLY':3,'ALA':5,'VAL':9,'PHE':9,'ILE':11,'LEU':11,'PRO':7,'MET':9,'ASP':4,'GLU':6,'LYS':13,'ARG':13,'SER':5,'THR':7,'TYR':9,'CYS':5,'ASN':6,'GLN':8,'HIS':8,'TRP':10}
RES_N0N ={'GLY':2,'ALA':3,'VAL':5,'PHE':9,'ILE':6,'LEU':6,'PRO':5,'MET':5,'ASP':4,'GLU':5,'LYS':6,'ARG':6,'SER':3,'THR':4,'TYR':9,'CYS':3,'ASN':4,'GLN':5,'HIS':6,'TRP':11}
RES_N0O ={'GLY':2,'ALA':3,'VAL':5,'PHE':9,'ILE':6,'LEU':6,'PRO':5,'MET':5,'ASP':4,'GLU':5,'LYS':6,'ARG':6,'SER':3,'THR':4,'TYR':9,'CYS':3,'ASN':4,'GLN':5,'HIS':6,'TRP':11}
RES_N0S ={'GLY':2,'ALA':3,'VAL':5,'PHE':9,'ILE':6,'LEU':6,'PRO':5,'MET':5,'ASP':4,'GLU':5,'LYS':6,'ARG':6,'SER':3,'THR':4,'TYR':9,'CYS':3,'ASN':4,'GLN':5,'HIS':6,'TRP':11}


ions = ["K+","Cl-","CL-","Na+","Mg+"]
#=======================================================================

# main classes coding

class pdb_atom:
	"""Represents a single atom from a PDB structure file.
	
	This class stores all relevant atomic properties read from a PDB file including
	coordinates, residue information, and chemical properties.
	
	Attributes:
		name (str): Atom name/identifier.
		Type (str): Atom type designation.
		element (str): Chemical element symbol.
		ptype (str): PDB atom type designation.
		xcoord (float): X-coordinate in Angstroms.
		ycoord (float): Y-coordinate in Angstroms.
		zcoord (float): Z-coordinate in Angstroms.
		num (int): Atom number/index.
		resNum (int): Residue number this atom belongs to.
		chain_t (str): Chain identifier.
		resTyp (str): Residue type/name.
		charge (float): Atomic charge.
		occ (float): Occupancy from PDB file.
		bfactor (float): B-factor/temperature factor.
	"""
	def __init__(self):
		self.name      = ''
		self.Type      = ''
		self.element   = ''
		self.ptype     = ''
		self.xcoord    = 0
		self.ycoord    = 0
		self.zcoord    = 0
		self.num       = 0
		self.resNum    = 0
		self.chain_t   = ''
		self.resTyp    = ''
		self.name      = ''
		self.charge    = 0
		self.occ       = 0
		self.bfactor   = 0

class residue:
	"""Represents a single amino acid residue in a protein structure.
	
	This class organizes atoms that belong to a single residue, including the
	backbone atoms and side chain atoms, along with properties like charge and
	hydrogen count.
	
	Attributes:
		name (str): Residue identifier (e.g., 'ALA15').
		typ (str): Residue type/name (e.g., 'ALA').
		num (int): Residue sequence number.
		alfaC (pdb_atom): Alpha carbon atom.
		Nitrogen (pdb_atom): Backbone nitrogen atom.
		carb (pdb_atom): Backbone carbon atom.
		oxygen (pdb_atom): Backbone oxygen atom.
		carbB (pdb_atom): Beta carbon atom (side chain start).
		r_atoms (list): All ordered atoms in residue.
		atomsNum (list): List of atom numbers belonging to this residue.
		hydrogen (int): Count of hydrogen atoms in residue.
		charge (float): Total charge of the residue.
		side_chain (list): Side chain atoms.
	"""

	def __init__(self):

		obj = pdb_atom

		self.name        = ''
		self.typ         = ''
		self.num         = ''
		self.alfaC       = obj
		self.Nitrogen    = obj
		self.carb		 = obj
		self.oxygen      = obj
		self.carbB       = obj
		self.r_atoms     = []
		self.atomsNum    = []
		self.hydrogen    = 0
		self.charge      = 0
		self.side_chain  = []

	def reorg(self):
		"""Reorganize atoms into backbone-first, then side-chain order.
		
		Populates r_atoms list with backbone atoms (N, CA, C, O, CB) followed
		by all side chain atoms in a consistent order suitable for analysis.
		"""

		self.r_atoms.append(self.Nitrogen)
		self.r_atoms.append(self.alfaC)
		self.r_atoms.append(self.carb)
		self.r_atoms.append(self.oxygen)
		self.r_atoms.append(self.carbB)

		for i in range(len(self.side_chain)):
			self.r_atoms.append(self.side_chain[i])
			#print(self.r_atoms[i].ptype)

class protein:
	"""Represents a protein structure parsed from a PDB (Protein Data Bank) file.
	
	This class handles reading PDB files, organizing atoms into residues, calculating
	geometric properties (center, bounding box), and provides methods for structure
	manipulation including water removal, ion pruning, charge assignment, and file output.
	
	Attributes:
		name (str): Path to the PDB file.
		chain (list): List of residue objects.
		resN (int): Number of residues in the chain.
		atoms (list): All pdb_atom objects in the structure.
		total_charge (float): Total charge of the protein.
		waters (list): Water molecule residues.
		up_vertice (list): Maximum coordinate box corner [x, y, z].
		down_vertice (list): Minimum coordinate box corner [x, y, z].
		protein_center (list): Geometric center of protein [x, y, z].
	"""

	def __init__(self,name,reorg=False):
		"""Initialize protein object and parse PDB file.
		
		Reads a PDB file and populates atomic and residue data structures.
		Calculates bounding box and geometric center of the structure.
		
		Args:
			name (str): Path to the PDB file to parse.
			reorg (bool, optional): Whether to reorganize atoms. Defaults to False.
		"""
		self.name           = name
		self.chain          = []
		self.resN           = 0
		self.atoms          = []
		self.total_charge   = 0
		self.waters         = []
		self.up_vertice     = [0,0,0]
		self.down_vertice   = [0,0,0]
		self.protein_center = [0,0,0]
		
		#---------------------------------------------------------------#
		# pdb parser 
		#---------------------------------------------------------------#
		
		i = 1
		pdb_file = open(self.name,'r')
		for line in pdb_file:
			if line[0:4]=="ATOM" or line[0:6]=="HETATM":
				a = pdb_atom()
				a.num     = i
				a.ptype   = line[12:16]
				a.resTyp  = line[17:20]
				a.chain_t = line[21:22]
				a.resNum  = int(line[22:26])
				a.xcoord  = float(line[30:38])
				a.ycoord  = float(line[38:46])
				a.zcoord  = float(line[46:54])
				a.occ     = float(line[56:60])
				a.bfactor = float(line[61:66])
				a.name    = a.Type + str(a.num)
				if "Na+" in a.ptype or "Cl-" in a.ptype or "Zn" in a.ptype:
					a.element = a.ptype[0:3]
				else:
					a.element = a.ptype[0:2]					
				if a.element[0] =="1" or a.element[0]=="2" or a.element[0]=="3":
					a.element = "H"
				elif a.element == "He":
					a.element = "H"	
				elif a.ptype == "OXT":
					a.ptype = "O"					
				self.atoms.append(a)				
				i+=1
		pdb_file.close()													
				
		#---------------------------------------------------------------
		# residue definition 
		# --------------------------------------------------------------
		
		r = 1
		resnum_i = self.atoms[0].resNum
		print(len(self.atoms),resnum_i)
		for i in range(len(self.atoms)):
			if  i==0 or not self.atoms[i].resNum == resnum_i:
				a      = residue()
				a.name = self.atoms[i].resTyp + str(r)
				a.typ  = self.atoms[i].resTyp
				a.num  = r
				a.atomsNum.append(self.atoms[i].num)
				self.chain.append(a)
				r+=1
				resnum_i = self.atoms[i].resNum
				self.atoms[i].resNum = r-1
			else:
				self.atoms[i].resNum = r-1

		
		'''
		for i in range(len(self.chain)):
			for j in range(len(self.atoms)):
				if self.chain[i].num == self.atoms[j].resNum:
					self.chain[i].atomsNum.append(self.atoms[j].num)
					if self.atoms[j].ptype == 'N':
						self.chain[i].Nitrogen = self.atoms[j]
					elif self.atoms[j].ptype == 'CA':
						self.chain[i].alfaC = self.atoms[j]
					elif self.atoms[j].ptype =='C':
						self.chain[i].carb = self.atoms[j]
					elif not self.chain[i].typ =='GLY' and self.atoms[j].ptype == 'CB':
						self.chain[i].carbB = self.atoms[j]
					elif self.chain[i].typ =='GLY' and self.atoms[j].ptype == 'H':
						self.chain[i].carbB = self.atoms[j]
					elif self.atoms[j].ptype == 'O' or self.atoms[j].ptype == 'OC1':
						self.chain[i].oxygen = self.atoms[j]
					elif self.atoms[j].element =='H':
						self.chain[i].hydrogen +=1
						self.chain[i].side_chain.append(self.atoms[j])
					else:
						self.chain[i].side_chain.append(self.atoms[j])

		if reorg == True:
			for i in range(len(self.chain)):
				self.chain[i].reorg()

		for i in range(len(self.chain)):
			for j in range(len(self.chain[i].r_atoms)):
				print(i,j,self.chain[i].typ,i+j,self.chain[i].r_atoms[j].ptype)

			cnt = 0
			for i in range(len(self.chain)):
				for j in range(len(self.chain[i].r_atoms)):
					self.atoms[cnt] = self.chain[i].r_atoms[j]
					cnt +=1
		'''

		self.resN = len(self.chain)
				
		#--------------------------
		#define geometric properties of the protein
		#--------------------------
		xc =[]
		yc =[]
		zc =[]
		for i in range(len(self.atoms)):
			xc.append(self.atoms[i].xcoord) 
			yc.append(self.atoms[i].ycoord) 
			zc.append(self.atoms[i].zcoord)
		
		self.up_vertice[0]   = max(xc)
		self.up_vertice[1]   = max(yc)
		self.up_vertice[2]   = max(zc)
		self.down_vertice[0] = min(xc)
		self.down_vertice[1] = min(yc)
		self.down_vertice[2] = min(zc)
		
		self.protein_center[0] = (self.up_vertice[0] + self.down_vertice[0])/2
		self.protein_center[1] = (self.up_vertice[1] + self.down_vertice[1])/2
		self.protein_center[2] = (self.up_vertice[2] + self.down_vertice[2])/2
				
		#print properties  
		print("PDB file: "+self.name)
		print("Number of atoms: "+str(len(self.atoms)))
		print("Number of residues: "+str(len(self.chain)))
		print(self.protein_center)
		print(self.up_vertice)
		print(self.down_vertice)
		
	#-------------------------------------------------------------------	
	
	def remove_atom(self,i):
		"""Remove a single atom by index.
		
		Args:
			i (int): Index of the atom to remove from the atoms list.
		"""
		del self.atoms[i]		
		
	#-------------------------------------------------------------------
	
	def remove_residue(self,i):
		"""Remove a residue and all its atoms by residue index.
		
		Args:
			i (int): Index of the residue to remove from the chain list.
		"""
		
		k=0
		while not k == len(self.atoms):
			if self.atoms[k].resNum == i:
				del self.atoms[k]
			k +=1
		del self.chain[i]		
		
	#-------------------------------------------------------------------

	def prune_pdb(self):
		"""Remove atoms with alternate locations (B records) from structure.
		
		Scans atoms for residue types starting with 'B' (indicating alternate
		conformations) and removes them.
		"""
		a = []
		for i in range(len(self.atoms)):
			if self.atoms[i].resTyp[0]=="B":
				a.append(i)
				
		for i in sorted(a,reverse=True):
			del self.atoms[i]
	
	#-------------------------------------------------------------------
	
	def prune_water(self,radius,res):
		"""Remove water molecules beyond a specified distance threshold.
		
		Identifies water molecules (HOH/WAT/SOL) that are beyond the specified
		radius from a reference point (protein center or specific residue) and
		removes them along with their atoms.
		
		Args:
			radius (float): Distance cutoff in Angstroms.
			res (int): Residue index for reference point (0 = protein center).
		"""
		reference = [self.protein_center[0],self.protein_center[1],self.protein_center[2]]
		if res > 0:
			reference[0] = self.atoms[self.chain[res-1].atomsNum[0]].xcoord
			reference[1] = self.atoms[self.chain[res-1].atomsNum[0]].ycoord
			reference[2] = self.atoms[self.chain[res-1].atomsNum[0]].zcoord
		rmv_list = []
		rmv_list_a = []
		dist_tmp = 0
		for i in range(len(self.chain)):
			if self.chain[i].typ == "HOH" or  self.chain[i].typ == "WAT" or self.chain[i].typ == "SOL":
				oxygen = self.atoms[self.chain[i].atomsNum[0]]
				xc = (oxygen.xcoord - reference[0])**2
				yc = (oxygen.ycoord - reference[1])**2
				zc = (oxygen.zcoord - reference[2])**2
				dist_tmp = math.sqrt(xc+yc+zc)
				if dist_tmp > radius:
					rmv_list.append(i-1)
					rmv_list_a.append(self.chain[i].atomsNum[0]-1)
					rmv_list_a.append(self.chain[i].atomsNum[0])
					rmv_list_a.append(self.chain[i].atomsNum[0]+1)
				#else:
				#	print(dist_tmp)
		
		for i in sorted(rmv_list,reverse=True):
			del self.chain[i]
			
		for i in sorted(rmv_list_a,reverse=True):
			del self.atoms[i]
			
			
	#-------------------------------------------------------------------
	
	def prune_ions(self):
		"""Remove all ionic species from the structure.
		
		Identifies and removes all atoms with residue types matching known ions
		(K+, Cl-, Na+, Mg+) from the structure.
		"""
		ions_rm    = []
		for i in range(len(self.atoms)):
			if self.atoms[i].resTyp in ions:
				ions_rm.append(i)
		
		for i in sorted(ions_rm,reverse=True):
			del self.atoms[i]
							
	#-------------------------------------------------------------------
						
	def split_complex(self,lign):
		"""Extract and save a ligand/component as a separate PDB file.
		
		Separates a specific residue type from the protein structure and writes
		it to a new PDB file with '_lig.pdb' suffix.
		
		Args:
			lign (str): Three-letter residue code of ligand to extract (e.g., 'ATP').
		"""
		lig = []
		atoms_swap = []
		a   = []
		for i in range(len(self.atoms)):
			if self.atoms[i].resTyp==lign:
				lig.append(self.atoms[i])
				a.append(i)
				
		for i in sorted(a,reverse=True):
			del self.atoms[i]
				
		input_text ="HEADER {0} pdb file\n".format(self.name)

		i=1
		for atom in lig:
			input_text += "ATOM {0:6} {1:4} {2:2} {3:<1} {4:<7} {5:7.3f} {6:7.3f} {7:7.3f} {8:>5.2f} {9:>4.2f} \n".format(i,atom.ptype,atom.resTyp,atom.chain_t,atom.resNum,atom.xcoord,atom.ycoord,atom.zcoord,atom.occ,atom.bfactor)
			i+=1

		pdb = open(self.name[:-4]+"_lig.pdb",'w')
		pdb.write(input_text)
		pdb.close()

	def remove_waters(self):
		"""Remove all water molecules from the structure.
		
		Identifies water residues (WAT, HOH, SOL) and removes them from both
		the atoms and chain lists.
		"""
		a = []
		for i in range(len(self.atoms)):
			if self.atoms[i].resTyp=="WAT" or self.atoms[i].resTyp=="HOH" or self.atoms[i].resTyp=="SOL":
				a.append(i)
				
		for i in sorted(a,reverse=True):
			del self.atoms[i]
				


	def charge_res(self):
		"""Assign formal charges to residues based on hydrogen count.
		
		Determines protonation state and formal charge of ionizable residues
		(ASP, GLU, LYS, ARG, HIS, CYS) based on their hydrogen atom count.
		Accounts for terminal residue special cases.
		"""

		for i in range(len(self.chain)):

			# for asp and glu
			if self.chain[i].typ == 'ASP' or self.chain[i].typ == 'GLU':
				if self.chain[i].num == self.chain[0].num:
					if self.chain[i].hydrogen   == RES_N0H[self.chain[i].typ].typ +1:
						self.chain[i].charge = -1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ].typ +2:
						self.charge = 0
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ].typ +3:
						self.chain[i].charge = 1
					else:
						continue
				elif self.chain[i].num == self.chain[-1].num:
					if self.chain[i].hydrogen   == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge = -2
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ] +1:
						self.chain[i].charge = -1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ] +2:
						self.chain[i].charge = 0
					else:
						continue
				else:
					if self.chain[i].hydrogen   == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge = -1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+1:
						self.chain[i].charge = 0
					else:
						continue
			#for lys and arg
			elif self.chain[i].typ == 'LYS' or self.chain[i].typ == 'ARG':
				if self.chain[i].num == self.chain[0].num:
					if self.chain[i].hydrogen   == RES_N0H[self.chain[i].typ]+1:
						self.chain[i].charge = 1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+2:
						self.chain[i].charge = 2
					else:
						continue
				elif self.chain[i].num == self.chain[-1].num:
					if self.chain[i].hydrogen   == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge = 0
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ] +1:
						self.chain[i].charge = 1
					else:
						continue
				else:
					if self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge = 1
					else:
						continue
			#for his
			elif self.chain[i].typ == 'HIS':
				if self.chain[i].num == self.chain[0].num:
					if self.chain[i].hydrogen   == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge=0
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+1:
						self.chain[i].charge =1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ] +2:
						self.chain[i].charge = 2
					else:
						continue
				elif self.chain[i].num == self.chain[-1].num:
					if self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge=0
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+1:
						self.chain[i].charge=1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]-1:
						self.chain[i].charge=-1
					else:
						continue
				else:
					if self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge = 1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]-1:
						self.chain[i].charge=0
					else:
						continue
			#for cys
			elif self.chain[i].typ == 'CYS':
				if self.chain[i].num == self.chain[0].num:
					if self.chain[i].hydrogen   == RES_N0H[self.chain[i].typ]+1:
						self.chain[i].charge = 0
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+2:
						self.chain[i].charge = 1
					else:
						continue
				elif self.chain[i].num == self.chain[-1].num:
					if self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge=-1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+1:
						self.chain[i].charge=0
					else:
						continue
				else:
					if self.chain[i].hydrogen == RES_N0H[self.chain[i].typ] or self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]-1:
						self.chain[i].charge=0
					else:
						continue
			#for any residue
			else:
				if self.chain[i].num == self.chain[0].num:
					if self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+2:
						self.chain[i].charge=1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+1:
						self.chain[i].charge=0
					else:
						continue
				elif self.chain[i].num == self.chain[-1].num:
					if self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge=-1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+1:
						self.chain[i].charge=0
					else:
						continue
				else:
					self.chain[i].charge = 0


	def write_xyz(self):
		"""Write structure to XYZ format file.
		
		Exports atomic coordinates in extended XYZ format (element, x, y, z).
		Output file is named with .xyz extension based on input pdb filename.
		"""

		input_text = '{0} \n \n'.format(len(self.atoms))
		xyz=open(self.name[:-4] +'.xyz','w')
		
		for atom in self.atoms:
			if atom.element[0] == "H":				
				input_text += '{0} {1} {2} {3} \n'.format(atom.element[0],atom.xcoord,atom.ycoord,atom.zcoord)
			else:
				input_text += '{0} {1} {2} {3} \n'.format(atom.element,atom.xcoord,atom.ycoord,atom.zcoord)


		xyz.write(input_text)
		xyz.close()

	def write_pdb(self,filename):
		"""Write structure to PDB format file.
		
		Exports the protein structure to a standard PDB format file with proper
		formatting for atom records including coordinates, occupancy, and B-factors.
		
		Args:
			filename (str): Path to the output PDB file.
		"""

		input_text ="HEADER {0} pdb file\n".format(self.name)

		i=1
		for atom in self.atoms:
			input_text += "ATOM {0:6} {1:4} {2:2} {3:<1} {4:<7} {5:7.3f} {6:7.3f} {7:7.3f} {8:>5.2f} {9:>4.2f} \n".format(i,atom.ptype,atom.resTyp,atom.chain_t,atom.resNum,atom.xcoord,atom.ycoord,atom.zcoord,atom.occ,atom.bfactor)
			i+=1

		pdb = open(filename,'w')
		pdb.write(input_text)
		pdb.close()

	def mopac_mop(self           ,
				  mode = "SP"    ,
				  mozyme = True  ,
				  solvent =True  ,
				  method ="PM7"  ):
		"""Generate MOPAC input file for semi-empirical quantum chemistry.
		
		Creates a MOPAC .mop input file for quantum mechanical calculations with
		specified method and solvation model. Supports MOZYME approximation for
		large systems.
		
		Args:
			mode (str, optional): Calculation mode. Defaults to "SP" (single point).
			mozyme (bool, optional): Enable MOZYME approximation. Defaults to True.
			solvent (bool, optional): Include solvation (eps=78.4). Defaults to True.
			method (str, optional): QM method (e.g., 'PM7', 'PM6'). Defaults to "PM7".
		"""

		sol  = ""
		mozy = ""
		if solvent:
			sol = "EPS=78.4 RSOLV1.3"
		
		if mozyme:
			mozy = "mozyme"


		input_text = "{0} 1SCF large aux allvecs {1} {2}\n\n".format(method,mozy,sol,cutoff)
		chain = ""

		i=1
		for atom in self.atoms:
			input_text += "{0} {1}  1  {2} 1 {3} \n".format(atom.element,atom.xcoord,atom.ycoord,atom.zcoord)
			i+=1


		input_file =open(self.name[:-4] + ".mop",'w')
		input_file.write(input_text)
		input_file.close()


Objeto_exemplo = protein("nome_do_arquivo.pdb") #criei o objeto
Objeto_exemplo.write_xyz() #usando método 
