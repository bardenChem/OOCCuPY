#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ReactionCoordinate module - Definition and management of reaction coordinates.

This module provides the ReactionCoordinate class for defining various types
of reaction coordinates (distance, angle, dihedral, etc.) used in scanned
simulations, umbrella sampling, and constraint-based calculations.

Authors: Igor Barden Grillo, contributors
"""

#FILE = ReactionCoordinate.py

#=============================================================================
#from tomlkit import key

from .commonFunctions import *
from pMolecule import *
#*****************************************************************************
class ReactionCoordinate:
	"""Define and manage reaction coordinates for constrained simulations.
	
	Handles definition of various reaction coordinate types (distance, angle,
	dihedral, multiple distance, etc.) with optional mass-weighted constraints.
	Supports automatic labeling from system structure and parameter definition
	for scanned and restricted simulations.
	
	Attributes:
		atomsSel (list): Original atom selection.
		atoms (list): Processed atom indices.
		nAtoms (int): Number of atoms in coordinate definition.
		Type (str): Coordinate type (Distance, Angle, Dihedral, etc.)
		increment (float): Scan increment size.
		minimumD (float): Initial/minimum value of coordinate.
		label (str): Human-readable coordinate label.
		label2 (str): Alternative coordinate label.
		weight13, weight31 (float): Weights for multiple distance coordinate.
	"""
	def __init__(self,_atoms,_massConstraint,_type="Distance"):
		"""Initialize ReactionCoordinate object.
		
		Args:
			_atoms (list): Atom indices or patterns defining coordinate.
			_massConstraint (bool): Apply mass-weighted constraint.
			_type (str): Coordinate type. Options: Distance, Angle, Dihedral,
				multipleDistance, Tether. Default: "Distance"
		"""
		self.atomsSel	    = _atoms
		self.atoms          = []
		self.nAtoms 		= len(_atoms)
		self.massConstraint = _massConstraint
		self.Type 			= _type
		self.weight13 		=  1.0
		self.weight31 		= -1.0
		self.period 		= 360.0
		self.increment      = 0.0
		self.minimumD  		= 0.0
		self.maximumD		= 0.0
		self.label 			= "Reaction Coordinate"
		self.label2         = "ReactionCoordinate"
		self.reaction_type 	= "unknown"
		self.nsteps		    = 0
		self.AB             = 0.0
		self.BC			    = 0.0

		for atom in self.atomsSel:
			try: self.atoms.append( atom )
			except:self.atoms.append( atom[0] )


		if self.Type == "Distance" or self.Type == "distance":
			if self.nAtoms == 3:
				self.Type = "multipleDistance"
	#==========================================================================================================
	def GetRCLabel(self,_molecule):
		"""Generate human-readable labels for the reaction coordinate.
		
		Extracts atom names and residue information from the molecular system
		to create informative labels for plots and output files.
		
		Args:
			_molecule (System): Molecular system for label extraction.
		"""
		sequence = getattr( _molecule, "sequence", None )		
		if self.Type == "multipleDistance":
			A1 = _molecule.atoms.items[ self.atoms[0] ]
			A2 = _molecule.atoms.items[ self.atoms[1] ]
			A3 = _molecule.atoms.items[ self.atoms[2] ]			
			if not sequence == None:
				A1res = A1.parent.label.split(".")
				A2res = A2.parent.label.split(".")
				A3res = A3.parent.label.split(".")
				#resName1 
				self.label =  A1.label + "(" + A1res[0] + A1res[1] + ")-"
				self.label += A2.label + "(" + A2res[0] + A2res[1] + ")--"
				self.label += A3.label + "(" + A3res[0] + A3res[1] + ") $\AA$"
			else: 
				self.label  = A1.label + "-" + A2.label +"-"+ A3.label
				self.label2 = A1.label + "-" + A2.label               
		elif self.Type == "Distance" or self.Type == "distance":			
			A1 = _molecule.atoms.items[ self.atoms[0] ]
			A2 = _molecule.atoms.items[ self.atoms[1] ]
			if not sequence == None:
				A1res = A1.parent.label.split(".")
				A2res = A2.parent.label.split(".")
				self.label =  A1.label + "(" + A1res[0] + A1res[1] + ")--"
				self.label += A2.label + "(" + A2res[0] + A2res[1] + ") $\AA$"	
				self.label2 =A1.label + "(" + A1res[0] + A1res[1] + ")--"
				self.label2+=A2.label + "(" + A2res[0] + A2res[1] + ")"
			else: 
				self.label  = A1.label + "-" + A2.label			
				self.label2 = A1.label + "-" + A2.label 
		#.--------------------------
		elif self.Type == "Dihedral":
			A1 = _molecule.atoms.items[ self.atoms[0] ]
			A2 = _molecule.atoms.items[ self.atoms[1] ]
			A3 = _molecule.atoms.items[ self.atoms[2] ]
			A4 = _molecule.atoms.items[ self.atoms[3] ]
			if not sequence == None:
				A1res = A1.parent.label.split(".")
				A2res = A2.parent.label.split(".")
				A3res = A3.parent.label.split(".")
				A4res = A4.parent.label.split(".")
				self.label =  A1.label + "(" + A1res[0] + A1res[1] + ")-"
				self.label += A2.label + "(" + A2res[0] + A2res[1] + ")-"
				self.label += A3.label + "(" + A3res[0] + A3res[1] + ")-"
				self.label += A4.label + "(" + A4res[0] + A4res[1] + ") $\AA$"
			else: self.label =  A1.label + "-" + A2.label +"-" + A3.label +"-"+A4.label + "$\AA$"
	#==================================================================================================
	def SetInformation(self,_molecule,_dincre,_dminimum=None,_sigma_pk1_pk3=None,_sigma_pk3_pk1=None):
		"""Define coordinate values and scanning parameters.
		
		Args:
			_molecule (System): Molecular system for coordinate calculation.
			_dincre (float): Scan increment size.
			_dminimum (float): Initial coordinate value. If None, calculated.
			_sigma_pk1_pk3 (float): Weight for first atom pair (multiple distance).
			_sigma_pk3_pk1 (float): Weight for second atom pair (multiple distance).
		"""	
		self.increment = _dincre		
		set_pars = True
		if not _dminimum == None:  
			self.minimumD = _dminimum
			set_pars      = False
		if not _sigma_pk1_pk3 == None: 	self.weight13 = _sigma_pk1_pk3
		if not _sigma_pk3_pk1 == None:	self.weight31 = _sigma_pk3_pk1		
		if set_pars: 
			
			# ============================================================
			# TWO-CENTER REACTIONS (Association/Dissociation)
        	# ============================================================
			if self.Type == "Distance" or self.Type == "distance":
				dist_initial = _molecule.coordinates3.Distance(self.atoms[0], self.atoms[1])
            
            	# Detect if we're starting with bonded or non-bonded atoms
            	# Typical bond lengths: H-H ~0.74 Å, C-H ~1.09 Å, O-H ~0.96 Å, C-C ~1.54 Å
				is_bonded = dist_initial < 2.0 

				if is_bonded:
					self.reaction_type = "dissociation"
					self.minimumD = dist_initial
					self.maximumD = dist_initial + 2.5 # Set maximum distance for dissociation
				else:
					self.reaction_type = "association"
					self.minimumD = dist_initial
					self.maximumD = dist_initial - 1.0 # Set maximum distance for association

				if self.maximumD <= self.minimumD:

					# Fallback logic
					self.maximumD = self.minimumD + 2.0 if is_bonded else max(self.minimumD - 1.0, 0.5)
			
			elif self.Type == "multipleDistance":
				if self.massConstraint:
					# Mass-weighted coordinates (for H-transfer)
					atomic_n1 = _molecule.atoms.items[self.atoms[0]].atomicNumber
					atomic_n3 = _molecule.atoms.items[self.atoms[2]].atomicNumber
					mass_a1 = GetAtomicMass(atomic_n1)
					mass_a3 = GetAtomicMass(atomic_n3)
					self.weight13 = mass_a1/(mass_a1+mass_a3)
					self.weight31 = mass_a3/(mass_a1+mass_a3)

					self.AB = _molecule.coordinates3.Distance(self.atoms[0], self.atoms[1])
					self.BC = _molecule.coordinates3.Distance(self.atoms[1], self.atoms[2])					
					# Transfer coordinate: weighted difference
					self.minimumD = (self.weight13 * self.AB) - (self.weight31 * self.BC )
					print(f"DEBUG: d(A-B) = {self.AB:.3f}, d(B-C) = {self.BC:.3f}")
					print(f"DEBUG: weight13 = {self.weight13:.3f}, weight31 = {self.weight31:.3f}")
					print(f"DEBUG: Calculated RC = {self.minimumD:.3f}")
					print(f"DEBUG: Expected RC = {self.AB - self.BC:.3f}")
					# For transfer, we scan from initial to final state
					self.reaction_type = "transfer"
					self.SetTransferLimits()
				else:
					# Simple coordinate: d(A-B) - d(B-C)
					self.AB = _molecule.coordinates3.Distance(self.atoms[0], self.atoms[1])
					self.BC = _molecule.coordinates3.Distance(self.atoms[1], self.atoms[2])
					self.minimumD = self.AB - self.BC					              
					# Transfer limits
					self.reaction_type = "transfer"
					self.SetTransferLimits()
			
			elif self.Type == "Dihedral": 
				self.minimumD = _molecule.coordinates3.Dihedral(self.atoms[0],self.atoms[1],self.atoms[2],self.atoms[3])
	
		self.DefineSteps()
	#--------------------------------------------------------------------------------------
	def GetDissociationLimit(self, bond_length):
		"""Set maximum distance for bond breaking based on bond type."""
    
   		# Different bonds break at different distances
		# Typically 1.5-2.5 times the equilibrium bond length
    
		if bond_length < 1.0:
		# Very short bond (H-H, etc.)
			return bond_length + 2.0
		elif bond_length < 1.3:
        # Typical C-H, N-H, O-H
			return bond_length + 2.0
		elif bond_length < 1.6:
		# C-C, C-N, C-O single bonds
			return bond_length + 2.2
		elif bond_length < 2.0:
			# Longer bonds (metal-ligand, etc.)
			return bond_length + 2.5
		else:
			# Unusually long - fallback
			return bond_length * 2.0  # Use multiplier to avoid "exploding"
		

	def GetAssociationLimit(self, current_distance):
		"""Set maximum approach distance for bond formation."""    
	
		# For association, we scan from current distance DOWN to bond length
		# The "maximum" in your coordinate might actually be the STARTING point
		# while "minimum" is the target bond length
    
    	# Typical bond lengths for common pairs (in Å)
		typical_bonds = {
			('H', 'H'): 0.74, ('H', 'C'): 1.09, ('H', 'N'): 1.01, ('H', 'O'): 0.96,
			('C', 'C'): 1.54, ('C', 'N'): 1.47, ('C', 'O'): 1.43,
			('N', 'N'): 1.45, ('N', 'O'): 1.40, ('O', 'O'): 1.48
		}
    
    	# Estimate target bond length
		target_bond_length = self.EstimateBondLength(current_distance, typical_bonds)
    
		if current_distance > target_bond_length:
        	# We're scanning towards shorter distances
        	# Current distance is maximum, target bond length is minimum
			return target_bond_length
		else:
        	# Already closer than typical bond - unusual
        	# Scan slightly closer to find optimal bond
			return max(current_distance - 1.0, 0.5)
		
	#--------------------------------------------------------------------------------------
	def SetTransferLimits(self):
		"""Set appropriate limits for 3-center transfer reactions."""
    
		# Correct swapped coordinate
		final_RC = self.weight13 * self.BC + self.weight31 * self.AB
    
		scan_range = abs(final_RC - self.minimumD)
		buffer = max(0.2, scan_range * 0.05)   # at least 0.2 Å buffer
    
		if final_RC > self.minimumD:
			self.minimumD = self.minimumD - buffer
			self.maximumD = final_RC + buffer	
		else:
			self.minimumD = final_RC - buffer
			self.maximumD = self.minimumD + buffer
    
    # Ensure we cross zero (important for TS)
		if self.minimumD > 0:
			self.minimumD = -abs(self.minimumD)
		if self.maximumD < 0:
			self.maximumD = abs(self.maximumD)
    
    # Optional: clamp only if extremely large (> 5 Å total)
		max_range = 5.0
		if abs(self.maximumD - self.minimumD) > max_range:
			print(f"  Large scan range ({abs(self.maximumD - self.minimumD):.2f}) for transfer")
			print(f"  Clamping to ±{max_range/2:.1f}")
			self.minimumD = max(self.minimumD, -max_range/2)
			self.maximumD = min(self.maximumD, max_range/2)
#---------------------------------------------------------------------------------------
	def EstimateBondLength(self, current_distance, typical_bonds_dict):
		"""Estimate expected bond length from atomic numbers."""
    
		# Try to get atom types from your molecule object
		# This assumes you have access to atomic numbers
		try:
			atom1 = self.atoms[0]
			atom2 = self.atoms[1]
			# You'll need to get element symbols from your molecule object
			# This is pseudocode - adapt to your actual data structure
			elem1 = self.GetElementSymbol(atom1)
			elem2 = self.GetElementSymbol(atom2)
        
			key = tuple(sorted([elem1, elem2]))
			if key in typical_bonds_dict:
				return typical_bonds_dict[key]
		except:
			pass
    
    	# Fallback based on current distance
		if current_distance < 1.5:
			return current_distance  # Assume we're scanning from near equilibrium
		else:
			return 1.5  # Generic bond length
		
	def DefineSteps(self):
		"""Calculate number of steps for the scan."""
    
		print("Defining steps for the scan...")
    
		# Calculate range
		if self.reaction_type == "association":
		# For association, we're going from large to small distance
			coordinate_range = abs(self.maximumD - self.minimumD)
		elif self.reaction_type == "dissociation":
			coordinate_range = abs(self.minimumD - self.maximumD)
		elif self.reaction_type == "transfer":
			coordinate_range = abs(self.maximumD - self.minimumD)
       

		# Calculate steps with safety check
		if self.increment > 0:
			nsteps = int(abs(coordinate_range / self.increment)) + 1
            
			# Safety: limit maximum steps to reasonable number
			max_steps = 36  # Prevent runaway calculations
			if nsteps > max_steps:
				print(f"Warning: {nsteps} steps exceeds maximum ({max_steps})")
				print(f"  Adjusting increment from {self.increment:.3f}")
				self.increment = coordinate_range / (max_steps - 1)
				nsteps = max_steps
				print(f"  New increment: {self.increment:.3f}")            
			self.nsteps = nsteps
		else:
			raise ValueError(f"Increment must be positive, got {self.increment}")
		
	#==================================================================================================
	def Print(self):
		"""Print reaction coordinate information to console."""
		print( "Printing reaction coordinate information:")
		print( "\tAtoms Indices: {}".format(self.atoms) )
		print( "\tLabel: {}".format(self.label) )
		print( "\tType: {}".format(self.reaction_type) )
		print( "\tWeight N1:{} ".format(self.weight13) )
		print( "\tWeight N2:{} ".format(self.weight31) )
		print( "\tIncrement:{} ".format(self.increment) )
		print( "\tInitial distance:{}".format(self.minimumD) )		
		print( "\tMaximum distance:{}".format(self.maximumD) )	
		if self.reaction_type == "transfer":
			print( "\tScan range: {:.2f} to {:.2f} ({} steps)".format(self.minimumD, self.maximumD, self.nsteps) )
			print( "\tInitial AB: {:.3f}, BC: {:.3f}".format(self.AB, self.BC) )
		elif self.reaction_type == "association":
			print( "\tScanning from {:.2f} to {:.2f} ({} steps)".format(self.maximumD, self.minimumD, self.nsteps) )	
		elif self.reaction_type == "dissociation":
			print( "\tScanning from {:.2f} to {:.2f} ({} steps)".format(self.minimumD, self.maximumD, self.nsteps) )


#==================================================================================
#=============================END OF FILE==========================================
#==================================================================================
