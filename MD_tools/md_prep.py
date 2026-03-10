#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""MD Preparation Module for OOCCuPY.

Provides utilities for preparing molecular dynamics simulations using AMBER
and GROMACS. This module handles ligand and protein preparation, parametrization
with AMBER, and structure optimization.

The module manages the complete MD preparation workflow:
- Extraction and parametrization of ligands/substrates
- Preparation of protein/enzyme structures
- Creation of solvated protein-ligand complexes
- Integration with AMBER and GROMACS for simulation setup

Attributes:
    cofac_list (list): List of known cofactors (ATP, NADPH, NAD+, etc.)
        with pre-computed parameters available in the cofac/ directory.
    atp_list (list): List of ATP ribose atom names for special handling.
    path_amber (str): Path to AMBER binary executables.
    path_cofac (str): Path to pre-computed cofactor parameter files.
"""

from pdb_class import *
from gmx_module import *
import parmed as pmd
import glob
import sys, os


# Known cofactors with pre-computed parameters
cofac_list = ["ATP", "atp", "NADPH", "NADP+", "NADH", "NAD+", "ADP"]
# ATP ribose atom names for special processing
atp_list = ["O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]

# AMBER and GROMACS binary paths
path_amber = "/home/igorchem/programs/amber18/bin"
Reduce = path_amber + "/reduce"  # Hydrogen addition/removal tool
pdb4 = path_amber + "/pdb4amber"  # PDB format converter
antech = path_amber + "/antechamber"  # Ligand parametrization
parmchk = path_amber + "/parmchk2"  # Parameter checking
tleap = path_amber + "/tleap"  # AMBER simulation preparation
pymol = "/home/igorchem/programs/bin/pymol"  # Molecular visualization
path_cofac = "/home/igorchem/OOCCuPY/cofac/"  # Pre-computed cofactor parameters


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def my_replace(fl, old, new):
    """Replace text in a file.

    Performs a simple find-and-replace operation on a text file.

    Args:
        fl (str): File path to edit.
        old (str): Text to find and replace.
        new (str): Replacement text.

    Note:
        This function should be relocated to a dedicated utilities library.
    """
    with open(fl, 'r+') as f:
        s = f.read()
        s = s.replace(old, new)
        f.seek(0)
        f.write(s)
			
#***********************************************************************
def fix_cofac_atoms(lig):
    """Fix atom types for cofactor molecules (ATP, NADPH, etc.).
    
    Processes cofactor PDB files to fix atom type nomenclature and rename
    residues appropriately. Replaces prime notation in atom names (e.g., O5')
    with asterisks (e.g., O5*) for AMBER compatibility.
    
    Args:
        lig (str): Ligand name/identifier (e.g., 'ATP').
    
    Returns:
        str: Modified ligand name for the processed PDB file.
    
    Raises:
        FileNotFoundError: If expected cofactor parameter files are missing.
    """
	new_atoms   = []
	result_text = ""
	LIG         = lig
	
	cofac = protein(lig)
	cofac.pdb_parse()
	
	
	if lig[:3]  == "ATP":
		LIG  = "atp"
		os.system("cp " + path_cofac +"atp* .")
	pdb_res  = open(LIG+".pdb",'w')
	
	for atom in cofac.atoms:
		atom.resType = "atp"
		if atom.ptype[-1:] == "'":
			atom.ptype   = atom.ptype[:-1]+"*"		
		new_atoms.append(atom)
	
	i=0
	for atom in new_atoms:
		result_text+= "ATOM {0:6} {1:4} {2:2} {3:<1} {4:<7} {5:7.3f} {6:7.3f} {7:7.3f} {8:>5.2f} {9:>4.2f} \n".format(i,atom.ptype,LIG,atom.chain_t,atom.resNum,atom.xcoord,atom.ycoord,atom.zcoord,atom.occ,atom.bfactor)
		i+=1
	pdb_res.write(result_text)
	pdb_res.close()

	return(LIG)
	
#***********************************************************************
def pdb_cat(pdb1, pdb2, flname):
    """Concatenate two PDB files into a single structure.
    
    Combines two protein/ligand structures into one PDB file, useful for
    creating protein-ligand complexes.
    
    Args:
        pdb1 (str): Path to first PDB file (usually protein).
        pdb2 (str): Path to second PDB file (usually ligand).
        flname (str): Output file name for concatenated structure.
    
    Note:
        This function should be refactored into the pdb_class module.
    """
	
	pdba = protein(pdb1) 
	pdbb = protein(pdb2) 
	pdba.pdb_parse()
	pdbb.pdb_parse()	

	result_text = "HEADER complex of {} {}".format(pdb1,pdb2)
	
	pdb_res = open(flname,"w")
	
	i = 0
	for atom in pdba.atoms:
		result_text+= "ATOM {0:6} {1:4} {2:2} {3:<1} {4:<7} {5:7.3f} {6:7.3f} {7:7.3f} {8:>5.2f} {9:>4.2f} \n".format(i,atom.ptype,atom.resTyp,atom.chain_t,atom.resNum,atom.xcoord,atom.ycoord,atom.zcoord,atom.occ,atom.bfactor)
		i+=1
	for atom in pdbb.atoms:
		result_text+= "ATOM {0:6} {1:4} {2:2} {3:<1} {4:<7} {5:7.3f} {6:7.3f} {7:7.3f} {8:>5.2f} {9:>4.2f} \n".format(i,atom.ptype,atom.resTyp,atom.chain_t,atom.resNum,atom.xcoord,atom.ycoord,atom.zcoord,atom.occ,atom.bfactor)
		i+=1	
		
	pdb_res.write(result_text)
	pdb_res.close
	

#=======================================================================
			
class md_prep:
    """Prepare molecular systems for AMBER/GROMACS molecular dynamics simulations.
    
    This class handles the complete workflow for preparing MD simulations:
    1. Extraction and parametrization of ligands
    2. Preparation of protein/enzyme structures
    3. Creation of solvated complexes
    4. Structure optimization and equilibration
    
    Attributes:
        pdb (str): Input PDB file path.
        current_pdb (str): Path to current PDB file being processed.
        lig (list): List of ligand/substrate names.
        net_charge (int): Net charge of the system.
        lig_charge (list): Charge of each ligand.
        num_lig (int): Number of ligands/substrates in system.
    """
    #-------------------------------------------
    
    def __init__(self, pdb):
        """Initialize the md_prep object.
        
        Args:
            pdb (str): Path to input PDB file containing protein structure.
        """
		self.pdb         = pdb
		self.current_pdb = pdb 
		self.lig         = []
		self.net_charge  = 0
		self.lig_charge  = []
		self.num_lig     = 0 
		
	#-------------------------------------------
	
    def prepare_lig(self, nlig, lign, chg, mult, rwat=True, lig_hy="False"):
        """Prepare and parametrize ligand(s) and extract from protein structure.
        
        Processes ligands by:
        1. Extracting from PDB structure
        2. Adding/removing hydrogens as needed
        3. Parametrizing with ANTECHAMBER or loading pre-computed parameters
        4. Generating AMBER library and force field files
        
        Args:
            nlig (int): Number of ligands to process.
            lign (list): List of ligand residue names.
            chg (list): List of ligand charges (one per ligand).
            mult (list): List of ligand multiplicities (one per ligand).
            rwat (bool, optional): Whether to remove water molecules. Defaults to True.
            lig_hy (str, optional): Whether to add hydrogens to ligands ('T'/'F'). 
                Defaults to "False".
        
        Note:
            Known cofactors (ATP, NADPH, etc.) use pre-computed parameters
            from the cofac/ directory. Unknown ligands are parametrized
            with ANTECHAMBER and charge obtained from BCC method.
        """
		
		lig_h = False
		if lig_hy == "T":
			lig_h = True
		pdb              = protein(self.pdb)		
		self.num_lig     = int(nlig)
		self.current_pdb = self.pdb
		
		print("Paramters parsed:\n")
		print("pdb file: " + self.pdb)
		print("Num of ligands: "+str(self.num_lig))
		
		for i in range(len(lign)):
			print("Lig #" + str(i)+": " + lign[i]) 
			print("Lig charge #" + str(i)+": " + str(chg[i])) 
			print("Lig multiplicity #" + str(i)+": " + str(mult[i])) 
			
		#pdb.pdb_parse()
		
		#if rwat:
		#	pdb.remove_waters()
		
		
		for j in range(self.num_lig):
			lig = []
			self.lig_charge.append(int(chg[j]))
			self.lig.append(lign[j])
			
			'''		
			for atom in pdb.atoms:
				if atom.resTyp == lign[j]:
					lig.append(atom)
			
			text_lig = "HEADER LIG\n"
			lig_pdb = open(lign[j]+".pdb",'w')
			i=0
			
			for atom in lig:
				text_lig += "ATOM {0:6} {1:4} {2:2} {3:<1} {4:<7} {5:7.3f} {6:7.3f} {7:7.3f} {8:>5.2f} {9:>4.2f} \n".format(i,atom.ptype,atom.resTyp,atom.chain_t,atom.resNum,atom.xcoord,atom.ycoord,atom.zcoord,atom.occ,atom.bfactor)
				i+=1			
			lig_pdb.write(text_lig)
			lig_pdb.close()
			'''
			
			
			print("=======================================================")
			print("Removing the ligand atoms from the provided PDB")
			if self.num_lig == 1:
				print("grep "+self.lig[j]+" "+ self.pdb+" > "+self.lig[j]+".pdb ")				
				os.system("grep "+self.lig[j]+" "+ self.pdb+" > "+self.lig[j]+".pdb ")
				print("grep -v "+self.lig[j]+" "+self.pdb+" > "+self.pdb[:-4]+"_w.pdb ")				
				os.system("grep -v "+self.lig[j]+" "+ self.pdb+" > "+self.pdb[:-4]+"_w.pdb ")
				self.current_pdb  = self.pdb[:-4]+"_w.pdb"
			elif self.num_lig > 1:
				print("grep -v "+self.lig[j]+" "+self.current_pdb+" > "+self.pdb[:-4]+"_w.pdb ")
				os.system("grep -v "+self.lig[j]+" "+ self.current_pdb+" > "+self.pdb[:-4]+"_"+str(j)+"_.pdb")
				self.current_pdb = self.pdb[:-4]+"_"+str(j)+"_.pdb"
			
			print("=======================================================")
			if lig_h:
				print("Removing and adding hydrogens in the ligand pdb")
				print(Reduce +"-Trim "+self.lig[j]+".pdb > " + self.lig[j]+"_h.pdb")
				os.system(Reduce+"-Trim "+self.lig[j]+".pdb > " + self.lig[j]+"_h.pdb")
				if not self.lig[j] in cofac_list:
					print(Reduce +self.lig[j]+"_h.pdb > " + self.lig[j]+".pdb")
					os.system(Reduce +self.lig[j]+"_h.pdb > " + self.lig[j]+".pdb")
				else:
					if self.lig[j] in cofac_list:
						os.system("cp " +self.lig[j]+"_h.pdb "+ self.lig[j]+".pdb")
			else:
				if self.lig[j] in cofac_list:
					print("Removing and adding hydrogens in the ligand pdb")
					print(Reduce +"-Trim "+self.lig[j]+".pdb > " + self.lig[j]+"_h.pdb")
					os.system(Reduce+"-Trim "+self.lig[j]+".pdb > " + self.lig[j]+"_h.pdb")	
					os.system("cp " +self.lig[j]+"_h.pdb "+ self.lig[j]+".pdb")
					
			os.system("cat < "+ self.lig[j]+".pdb")	
			input()
			if self.lig[j] in cofac_list:
				os.system( "cp " + path_cofac + "*lib "+	os.getcwd() )		
				os.system( "cp " + path_cofac + "*frcmod "+	os.getcwd() )
				os.system( "cp " + path_cofac + "*prep   "+	os.getcwd() )
				
				print("=======================================================")
				print("Ligand parameters will be loaded instead of created with ANTECHAMBER")
				self.lig[j] = fix_cofac_atoms(lign[j]+".pdb")	
					
				fl = os.listdir('.')
				print(fl)
				if self.lig[j]+".frcmod" in fl:
					print("FRCMOD OK...")
				else:				
					print("FRCMOD not found")				
					sys.exit()	
				if self.lig[j] +".lib" in fl:		
					print("LIB OK...")
					self.lig[j] = self.lig[j]+".pdb"
				else:
					print("=======================================================")	
					print("Creating tleap input to save ligand library")
					tleap_in =  "source leaprc.gaff2\n"
					tleap_in += "loadamberparams " +self.lig[j][:-4]+ ".frcmod\n"
					tleap_in +=  self.lig[j][:-4]+" = loadPdb "+self.lig[j]+"\n"
					tleap_in += "check "+self.lig[j][:-4]+"\n"					
					tleap_in += "saveoff " +self.lig[j][:-4]+" "+self.lig[j][:-4]+".lib \n"		
					tleap_in += "quit"
					tleap_file = open("tleap_in"+"_"+self.lig[j][:-4],'w')
					tleap_file.write(tleap_in)
					tleap_file.close()
					print("=======================================================")
					print("Run tleap and save the library with parameter ligands.")
					print(tleap + " -f tleap_in")
					os.system(tleap + " -f tleap_in"+"_"+self.lig[j][:-4])
				
			else:
				print("=======================================================")
				print("Ligand parameters will be created with ANTECHAMBER.")		
			
				par = False
				fl = os.listdir('.')
				if self.lig[j]+".frcmod" in fl:
					print("Found parameters for " + self.lig[j])
					print("FRCMOD file found for this ligand, antechamber parametrization will be skipped!") 
					self.lig[j] = self.lig[j]+".pdb"
					par = True
										
				if not par:
					print("===================================================")
					print("Run ANTECHAMBER:")
					print(antech+" -i "+self.lig[j]+".pdb -fi pdb -o "+self.lig[j]+".mol2 -fo mol2 -c bcc -nc "+chg[j]+" -m "+mult[j])
					os.system(antech + " -i " + self.lig[j]+".pdb -fi pdb -o " + self.lig[j]+".mol2  -fo mol2 -c bcc -nc "+chg[j]+" -m "+mult[j] )
					os.system("rm ANTECHAMBER*")
					print("===================================================")
					print("Run Pamchek and generate frcmod")				
					print(parmchk+" -i "+self.lig[j]+".mol2 -f mol2 -o " +self.lig[j]+".frcmod")
					print(parmchk+" -i "+self.lig[j]+".mol2 -f mol2 -o " +self.lig[j]+".frcmod")
					os.system(parmchk + " -i "+ self.lig[j]+".mol2 -f mol2 -o " + self.lig[j]+".frcmod")
					input()			
					print("=======================================================")	
					print("Creating tleap input to save ligand library")
					tleap_in = "source leaprc.gaff2 \n"
					tleap_in += "loadamberparams " +self.lig[j]+ ".frcmod\n"
					tleap_in +=  self.lig[j]+" = loadmol2 "+self.lig[j]+".mol2\n"
					tleap_in += "check "+self.lig[j]+"\n"					
					tleap_in += "saveoff " +self.lig[j]+" "+self.lig[j]+".lib \n"		
					tleap_in += "quit"
		
					tleap_file = open("tleap_in"+"_"+self.lig[j],'w')
					tleap_file.write(tleap_in)
					tleap_file.close()
					print("=======================================================")
					print("Run tleap and save the library with parameter ligands.")
					print(tleap + " -f tleap_in")
					os.system(tleap + " -f tleap_in"+"_"+self.lig[j])
					self.lig[j] = self.lig[j]+".pdb"
			
		os.system("sed 's/OXT/O  /' "+self.current_pdb+" > "+self.pdb[:-4]+"_wl.pdb ")
		self.current_pdb = self.pdb[:-4] +"_wl.pdb"

    def build_complex(self, addH=True):
        """Build complete solvated protein-ligand complex for simulation.
        
        Performs the following steps:
        1. Adds hydrogens to protein structure using Reduce
        2. Converts PDB format for AMBER compatibility
        3. Concatenates protein and ligand(s) into complex
        4. Generates AMBER topology and coordinate files
        5. Converts to GROMACS format
        
        Args:
            addH (bool, optional): Whether to add hydrogens to protein. 
                Defaults to True.
        
        Generates:
            prmtop: AMBER topology file
            inpcrd: AMBER coordinate file
            top: GROMACS topology file
            gro: GROMACS coordinate file
        """
        
        print("=======================================================")
		print("Preparing Receptor/enzyme!")
		print(Reduce +" -Trim "+self.current_pdb+  " > " +self.current_pdb[:-4]+"_p.pdb") 
		os.system(Reduce+" -Trim "+self.current_pdb+" > "+self.current_pdb[:-4]+"_p.pdb")
		print("=======================================================")
		print(pdb4 +  self.current_pdb[:-4]+"_p.pdb"+ " > " +self.current_pdb[:-4]+"_c.pdb")
		os.system(pdb4 +  self.current_pdb[:-4]+"_p.pdb"+ " > " +self.current_pdb[:-4]+"_c.pdb")
		#os.system(pdb4 +  self.current_pdb+" > " +self.current_pdb[:-4]+"_c.pdb")
		
		self.current_pdb = self.current_pdb[:-4]+"_c.pdb"
		input()
		print("=======================================================")
		print("Concatenating Receptor/Enzyme with ligand/substrate")
		for i in range(self.num_lig):
			print(self.lig[i])
			pdb_cat(self.current_pdb,self.lig[i],self.current_pdb)	
						
		os.rename(self.current_pdb,self.pdb[:-4]+"_comp.pdb")
		self.current_pdb = self.pdb[:-4]+"_comp.pdb"

		print("=======================================================")	
		#Creating tleap input to save ligand library
		tleap_in =  "source oldff/leaprc.ff99SB \n"		
		tleap_in += "source leaprc.gaff2 \n"		
		tleap_in += "source leaprc.water.tip3p \n"
		for i in range(self.num_lig):
			tleap_in += "loadamberparams " + self.lig[i][:-4] + ".frcmod\n"
			tleap_in += "loadoff " + self.lig[i][:-4] + ".lib\n"
		tleap_in += "complex = loadPdb " + self.current_pdb+ " \n"
		tleap_in += "solvatebox complex TIP3PBOX 12.0 \n"
		tleap_in += "check complex \n"
		tleap_in += "addions2 complex Na+ 0\n"
		tleap_in += "addions2 complex Cl- 0\n"
		tleap_in += "savePdb complex "+self.current_pdb+"\n"
		tleap_in += "saveamberparm complex " +self.pdb[:-4]+".prmtop "+ self.pdb[:-4] +".inpcrd\n"
		tleap_in += "quit"

		tleap_file = open('tleap_in','w')
		tleap_file.write(tleap_in)
		tleap_file.close()
		
		print(tleap + " -f tleap_in")
		os.system(tleap + " -f tleap_in")
		
		#---------------------------------------------------------------
		
		fl = os.listdir('.')
		gromp = False
		print(self.pdb[:-4]+".top")
		if self.pdb[:-4]+".top" in fl:
			print("Found gromacs parameters for " + self.pdb)
			gromp = True
		elif not self.pdb[:-4]+".prmtop" in fl:
			print("prmtop not in folder for gromacs conversion")
			sys.exit()
		else:
			print("===================================================")
			print("Saving topologies for gromacs")
			ap = pmd.load_file(self.pdb[:-4]+".prmtop",self.pdb[:-4] +".inpcrd")
			ap.save(self.pdb[:-4]+".top")
			ap.save(self.pdb[:-4]+".gro")

    def min_gromacs(self):
        """Perform energy minimization and equilibration with GROMACS.
        
        Executes the following GROMACS workflow:
        1. Structure minimization (EM)
        2. NVT equilibration (constant volume)
        3. NPT equilibration (constant pressure)
        4. MD production setup
        
        Requires pre-generated topology and coordinate files in GROMACS format.
        Generates .mdp input files (em.mdp, nvt.mdp, npt.mdp, md.mdp) and
        outputs final structures ready for production MD.
        """
        
        gromacs_inp()
		os.system("sed 's/WAT/SOL/' "+self.current_pdb+" > "+self.current_pdb[:-4]+"_t.pdb")
		os.rename(self.current_pdb[:-4]+"_t.pdb",self.current_pdb)
		os.system("sed 's/WAT/SOL/' "+self.pdb[:-4]+".top > "+self.pdb[:-4]+"_t.top")
		os.rename(self.pdb[:-4]+"_t.top",self.pdb[:-4]+".top")	
		os.system("sed 's/WAT/SOL/' "+self.pdb[:-4]+".gro > "+self.pdb[:-4]+"_t.gro")
		os.rename(self.pdb[:-4]+"_t.gro",self.pdb[:-4]+".gro")
		
		print("=======================================================")
		print("Preparing the gromacs input files for structure minimization.")
		text_to_run = "/usr/bin/gmx_d" +" grompp -f em.mdp -c "+self.pdb[:-4]+ " -p "+ self.pdb[:-4] +".top -o em.tpr -maxwarn 50"
		os.system(text_to_run)
		
		print("=======================================================")
		print("Running minimization in gromacs.")
		text_to_run = "/usr/bin/gmx_d" + " mdrun -v -deffnm em"
		os.system(text_to_run)	
		
		print("=======================================================")
		print("Writting minimized structure pdb")
		print("/usr/bin/gmx_d editconf -f em.gro -o "+self.pdb[:-4]+"_min.pdb")
		os.system("gmx editconf -f em.gro -o "+self.pdb[:-4]+"_min.pdb")	
		
		print("=======================================================")
		print("Preparing gromacs NVT equilibration")
		print("gmx grompp -f nvt.mdp -c em.gro -p "+self.pdb[:-4]+".top -o nvt.tpr -maxwarn 50")
		os.system("gmx grompp -f nvt.mdp -c em.gro -p "+self.pdb[:-4]+".top -o nvt.tpr -maxwarn 50")
		
		print("=======================================================")
		print("Running nvt equilibration in gromacs.")
		text_to_run = "/usr/bin/gmx_d" + " mdrun -v -deffnm nvt"
		os.system(text_to_run)	
		
		print("=======================================================")
		print("Preparing gromacs NPT equilibration")
		print("gmx grompp -f npt.mdp -c nvt.gro -p "+self.pdb[:-4]+".top -o npt.tpr -maxwarn 50")
		os.system("gmx grompp -f npt.mdp -c nvt.gro -p "+self.pdb[:-4]+".top -o npt.tpr -maxwarn 50")
		
		print("=======================================================")
		print("Running npt equilibration in gromacs.")
		text_to_run = "/usr/bin/gmx_d" + " mdrun -v -deffnm npt"
		os.system(text_to_run)	
		
		print("=======================================================")
		print("Preparing gromacs NPT equilibration")
		print("gmx grompp -f md.mdp -c npt.gro -p "+self.pdb[:-4]+".top -o md.tpr -maxwarn 50")
		os.system("gmx grompp -f md.mdp -c npt.gro -p "+self.pdb[:-4]+".top -o md.tpr -maxwarn 50")
		
		
    def production(self):
        """Start production molecular dynamics simulation.
        
        Note:
            This method is not yet implemented.
        """
        pass 
		
    def organize_files(self):
        """Organize output files into directory structure.
        
        Note:
            This method is not yet implemented.
        """
        pass
		
    def process_traj(self):
        """Process molecular dynamics trajectory files.
        
        Note:
            This method is not yet implemented.
        """
        pass
		
		

