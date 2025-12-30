#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = MopacQCMMinput.py


from .commonFunctions import *
from pMolecule import *
from pBabel import ImportCoordinates3
from pCore import Unpickle

#************************************************************
class MopacQCMMinput:
	'''
	Class to set methods to creat inputs for run QC/MM in mopac
	'''
	#---------------------------------------------------------
	def __init__(self,_parameters):
		'''
		'''
		self.system 		= _parameters["active_system"]
		self.baseName       = _parameters["basename"]
		self.QCatoms  		= []		
		self.gradVectors 	= []
		self.molinFile      = None 
		self.inputFile 		= None
		self.atomsDict		= {}
		self.pars           = None
		self.charges        = []

		self.Check_Pars(_parameters)

		if hasattr(self.system,"mmState"): 
			self.charges = self.system.mmState.charges
		else: self.charges[len(self.system.atoms)]

		if self.pars["cood_name"]:
			self.system.coordinates3 = Unpickle(self.pars["cood_name"])
		
		for i in self.system.atoms.items:
			symbol = GetAtomicSymbol( i.atomicNumber )
			index  = i.index
			x      = self.system.coordinates3[i.index, 0]
			y      = self.system.coordinates3[i.index, 1]
			z      = self.system.coordinates3[i.index, 2]
			self.atomsDict[index] = [ symbol,x,y,z,self.charges[index] ]
		self.QCatoms = list(self.system.qcState.qcAtoms)
		self.BAatoms = list(self.system.qcState.boundaryAtoms)

	#------------------------------------------------------------------
	def Check_Pars(self,_parameters):
		'''
		'''
		self.pars = { 
			"cood_name": None ,
			"Hamiltonian":"am1",
			"QCcharge":0       ,
			"multiplicity":1   ,
			"keywords":[]      , 
		}
		for key in _parameters.keys(): self.pars[key] = _parameters[key]

	#==================================================================
	def CalculateGradVectors(self):
		'''
		Calculate the grad vectors for the mol.in
		'''
		PHI      = 0.0 
		distance = 0.0
		indx     = 0
		if hasattr(self.system,"mmState"): 		
			#----------------------------------		
			for j in range( len(self.QCatoms) ):
				indx = 0 
				for i in self.system.atoms.items:						 
					distance = self.system.coordinates3.Distance( indx, self.QCatoms[j] )
					if not indx == self.QCatoms[j]:
						PHI 	+= self.charges[indx]/ distance
						indx    += 1
			
				PHI *= 332
				self.gradVectors.append(PHI)
				PHI=0		
		
	#===================================================================
	def write_input(self,_crd_name):
		'''
		Write the input files and grad vectors file
		'''
		#-----------------------------------------------------
		MULT        = "singlet"
		if 	 self.pars["multiplicity"] == 2: MULT = "doublet"
		elif self.pars["multiplicity"] == 3: MULT = "triplet"
		elif self.pars["multiplicity"] == 4: MULT = "quartet"
		elif self.pars["multiplicity"] == 5: MULT = "quintet"
		#-----------------------------------------------------	
		# creating files objects and setting paths	
		mol_file_name = os.path.join( self.baseName,"mol.in")
		self.mop_file_name = os.path.join(self.baseName, _crd_name[:-4] + "_" + self.pars["Hamiltonian"]+ ".mop" )
		mol_file  = open( mol_file_name, "w" )
		mop_file  = open( self.mop_file_name, "w" )
		# Getting info to write PDB file 
		sequence = getattr( self.system, "sequence", None )
		if sequence is not None: pdb_file  = open( self.mop_file_name[:-4]+".pdb","w")
		# Writting the files
		molInText = "\n{} 0\n".format( len(self.QCatoms) )
		mop_text  = self.pars["Hamiltonian"] + " 1SCF charge={} {} ".format(self.pars["QCcharge"],MULT)
		pdb_text  = ""
		for _key in self.pars["keywords"]:
			mop_text += _key
			mop_text += " "
		#-------------------
		cnt = 1
		mop_text+="\n\n\n"

		a1    = ""
		a2    = ""
		A1res = "UKN" 
		A2res = "UKN"
		for i in self.QCatoms:
			if sequence is not None:
				a1    = self.system.atoms.items[i]
				A1res = a1.parent.label.split(".")
				if i in self.BAatoms:
					a2    = self.system.atoms.items[ self.QCatoms[cnt+1] ]
					A2res = a2.parent.label.split(".")	
					pdb_text +=  "ATOM {0:6} {1:4} {2:2} {3:<1} {4:<7} {5:7.3f} {6:7.3f} {7:7.3f} {8:>5.2f} {9:>4.2f} \n".format(cnt,"H",A2res[0],"A",A2res[1],self.atomsDict[i][1],self.atomsDict[i][2],self.atomsDict[i][3],1.00,0.00)
				else:
					pdb_text +=  "ATOM {0:6} {1:4} {2:2} {3:<1} {4:<7} {5:7.3f} {6:7.3f} {7:7.3f} {8:>5.2f} {9:>4.2f} \n".format(cnt,a1.label,A1res[0],"A",A1res[1],self.atomsDict[i][1],self.atomsDict[i][2],self.atomsDict[i][3],1.00,0.00)
			if i in self.BAatoms:								
				mop_text += "{} {} 1 {} 1 {} 1\n".format("H",self.atomsDict[i][1],self.atomsDict[i][2],self.atomsDict[i][3])
			else:
				mop_text += "{} {} 1 {} 1 {} 1\n".format(self.atomsDict[i][0],self.atomsDict[i][1],self.atomsDict[i][2],self.atomsDict[i][3])
			cnt+=1
		idx = 0 
		for i in self.QCatoms:
			molInText += "{} {} {} {} {}\n".format(self.atomsDict[i][0],self.atomsDict[i][1],self.atomsDict[i][2],self.atomsDict[i][3],self.gradVectors[idx])
			idx += 1
		# writing to objects
		mop_file.write(mop_text)
		mop_file.close()
		mol_file.write(molInText)
		mol_file.close()
		if sequence is not None:
			pdb_file.write(pdb_text)
			pdb_file.close()
	#--------------------------------------------------------
	def Execute(self, mopac_path="/opt/mopac/bin/mopac"):
		'''
		'''

		command = mopac_path + " " + self.mop_file_name 
		if os.path.exists(mopac_path): 
			os.system(command)
		else: 
			mopac_path =  "/opt/apps/mopac/2016/bin/mopac "
			command = mopac_path + self.mop_file_name
			os.system(command)
	#--------------------------------------------------------
	def GetEnergy(self):
		'''
		Read ARC file to get and return the total energy
		'''
		arcfile = open(self.mop_file_name[:-4]+".arc","r")
		energy = 0.0 		
		for line in arcfile:
			line2 = line.split()
			if len(line2) == 9:
				if line2[0] == "HEAT" and line2[2] == "FORMATION":
					energy = 4.184*float(line2[4])
					break

		arcfile.close()
		return(energy)

#======================================================================================================#
#======================================END OF FILE=====================================================#
#======================================================================================================#