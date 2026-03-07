#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SimulationSystem module - Molecular system wrapper and orchestration.

This module provides the SimulationSystem class which wraps the pDynamo System class
with additional functionality for system setup, QM/MM region definition, reaction
coordinate management, and file I/O. Supports multiple input formats including
pDynamo PKL, AMBER, GROMACS, and coordinate files.

Authors: Igor Barden Grillo, contributors
"""

#FILE = CoreInterface.py


import os, glob, sys
#--------------------------------------------------------------
#Loading own libraries
from .commonFunctions import *
from .Simulation import Simulation
from .Analysis import *
from .QuantumMethods import *
#--------------------------------------------------------------
#loading pDynamo Libraries
from pBabel                    import *                                     
from pCore                     import *                                     
from pMolecule                 import *                              
from pMolecule.MMModel         import *
from pMolecule.NBModel         import *                                     
from pMolecule.QCModel         import *
from pSimulation               import PruneByAtom

from pBabel            import ExportSystem                    , \
                              GromacsDefinitionsFileReader    , \
                              GromacsParameterFileReader      , \
                              ImportCoordinates3              , \
                              ImportSystem
from .ReactionCoordinate import ReactionCoordinate
#==========================================================================

#**************************************************************************
class SimulationSystem:
    """Wrapper for pDynamo System class with enhanced functionality.
    
    Manages molecular system setup, coordinate import from various formats,
    QM/MM region definition, quantum method assignment, reaction coordinate
    management, and spherical pruning. Provides convenient class methods for
    initialization from multiple file formats.
    
    Attributes:
        baseName (str): Base name derived from input file.
        label (str): System label/identifier.
        system (System): The underlying pDynamo System object.
        Hybrid (bool): Whether system is QM/MM hybrid.
        quantumRegion (list): Atom indices in quantum region.
        protein (bool): Whether system is a protein.
        reactionCoordinates (list): List of ReactionCoordinate objects.
        refEnergy (float): Reference energy for comparisons.
        rcs (int): Number of reaction coordinates defined.
    """  
    #.-------------------------------------------------------------------------
    def __init__(self,_label="No specified"):
        """Initialize SimulationSystem instance.
        
        Args:
            _label (str): System label/identifier. Default: "No specified"
        """        
        self.baseName            = None
        self.label               = _label
        self.system              = None # Instance of the System pDnamo class
        self.Hybrid              = None 
        self.quantumRegion       = []
        self.protein             = False 
        self.reactionCoordinates = [] 
        self.refEnergy           = 0.0 
        self.rcs                 = 0 

    #===================================================================================
    @classmethod
    def From_PKL(selfClass,_pklPath,_FolderName=None):
        """Initialize SimulationSystem from a pDynamo PKL file.
        
        Args:
            _pklPath (str): Path to pDynamo PKL system file.
            _FolderName (str): Optional folder name (unused).
            
        Returns:
            SimulationSystem: Initialized instance with system loaded from PKL.
        """
        self                = selfClass()
        self.system         = ImportSystem(_pklPath)
        _name               = os.path.basename(_pklPath)
        self.baseName       = _name[:-4]
        if not self.system.nbModel:      
            if self.system.symmetryParameters is not None: 
                self.system.DefineNBModel( NBModelCutOff.WithDefaults( ) )
            else:
             self.system.DefineNBModel( NBModelFull.WithDefaults( ) )

        # test if is quantum, if hass mmModel and NbModel
        return(self)
    #=================================================================================== 
    @classmethod
    def From_AMBER(selfClass,_topologyFile,_coordinateFile):
        """Initialize SimulationSystem from AMBER topology and coordinate files.
        
        Args:
            _topologyFile (str): Path to AMBER topology (.prmtop) file.
            _coordinateFile (str): Path to coordinate (.inpcrd) file.
            
        Returns:
            SimulationSystem: Initialized instance with AMBER system imported.
        """
        self = selfClass()
        self.NBmodel = NBModelCutOff.WithDefaults()      
        self.system               = ImportSystem(_topologyFile)
        self.system.coordinates3  = ImportCoordinates3(_coordinateFile)
        self.system.DefineNBModel( NBModelCutOff.WithDefaults () )
        self.NBmodel = self.system.nbModel
        _name                     = os.path.basename(_topologyFile)
        self.baseName             = _name[:-4]        
        return(self)     
    #===================================================================================
    @classmethod
    def From_Gromacs(selfClass,_topologyFile,_coordinateFile):
        """Initialize SimulationSystem from GROMACS topology and coordinate files.
        
        Args:
            _topologyFile (str): Path to GROMACS topology (.top) file.
            _coordinateFile (str): Path to coordinate (.gro) file.
            
        Returns:
            SimulationSystem: Initialized instance with GROMACS system imported.
        """
        self = selfClass()
        parameters   = GromacsParameterFileReader.PathToParameters ( _topologyFile )
        print(parameters)
        input()
        self.system  = GromacsDefinitionsFileReader.PathToSystem   ( _topologyFile, parameters = parameters )
        self.system.coordinates3 = ImportCoordinates3              ( _coordinateFile )
        if self.system.symmetryParameters is not None: self.system.DefineNBModel ( NBModelCutOff.WithDefaults ( ) )
        else                                         : self.system.DefineNBModel ( NBModelFull.WithDefaults   ( ) )
        self.baseName = os.path.basename(_coordinateFile[:-4])
        self.NBmodel  = self.system.nbModel
        self.system.nbModel.Summary()
        return(self)        
    #===================================================================================
    @classmethod
    def From_Coordinates(selfClass,_coordinateFile):
        """Initialize SimulationSystem from a coordinate file.
        
        Args:
            _coordinateFile (str): Path to coordinate file (PDB, MOL2, etc.).
            
        Returns:
            SimulationSystem: Initialized instance with system from coordinates.
        """
        self = selfClass()        
        self.system   = ImportSystem(_coordinateFile)
        _name         = os.path.basename(_coordinateFile)
        self.baseName = _name[:-4]
        return(self) 
    #===================================================================================
    @classmethod
    def Protein_From_Coordinates(selfClass,_coordinateFile,_modelNumber=1):
        """Initialize protein system from coordinate file with OPLS force field.
        
        Args:
            _coordinateFile (str): Path to protein coordinate file.
            _modelNumber (int): Model number to use if multiple models. Default: 1
            
        Returns:
            SimulationSystem: Initialized protein instance with OPLS parameters.
        """
        self = selfClass()
        self.system      = ImportSystem(_coordinateFile, modelNumber = _modelNumber, useComponentLibrary = True)
        self.system.DefineMMModel( MMModelOPLS.WithParameterSet("protein") )
        self.system.DefineNBModel( NBModelCutOff.WithDefaults() )                   
        _name        = os.path.basename(_coordinateFile)
        self.baseName= _name[:-4]
        self.protein = True
        return(self) 
    #====================================================================================
    def Check(self):
        """Calculate and return single-point energy of the system.
        
        Returns:
            float: Total system energy.
        """
        self.refEnergy = self.system.Energy(doGradients = False)
        self.system.Summary()
        return self.refEnergy    
    #====================================================================================
    def Spherical_Pruning(self,_centerAtom,_radius,_DEBUG=False):
        """Perform spherical pruning around a center atom.
        
        Removes atoms outside specified radius from a reference atom,
        keeping only a spherical region of interest.
        
        Args:
            _centerAtom (str): Atom pattern for center atom.
            _radius (float): Sphere radius in Angstroms.
            _DEBUG (bool): Print debug information. Default: False
        """
        #---------------------------------------------------
        oldSystem = Clone(self.system)
        #---------------------------------------------------
        atomref      = AtomSelection.FromAtomPattern( oldSystem, _centerAtom )
        core         = AtomSelection.Within(oldSystem,atomref,_radius)
        core2        = AtomSelection.ByComponent(oldSystem,core)
        #---------------------------------------------------
        newLabel    = self.label + "_pruned"
        NBmodel     = self.system.nbModel
        self.system = None
        self.system = PruneByAtom( oldSystem,Selection(core2) )
        self.system.DefineNBModel( NBmodel )      
        self.label  = newLabel
        

    #======================================================================================
    def Setting_Free_Atoms(self,_centerAtom,_radius,_DEBUG=False):
        """Set mobile atoms within specified radius from center atom.
        
        Atoms within radius are marked as mobile; atoms outside are frozen.
        
        Args:
            _centerAtom (str): Atom pattern for center atom.
            _radius (float): Radius in Angstroms.
            _DEBUG (bool): Print debug information. Default: False
        """
        #-----------------------------------------------------
        atomref = AtomSelection.FromAtomPattern(self.system, _centerAtom)
        core    = AtomSelection.Within(self.system,atomref,_radius)        
        mobile  = AtomSelection.ByComponent(self.system,core)  
        
        #-----------------------------------------------------
        newLabel= self.system.label + "_fixed"       
        #------------------------------------------------------        
        self.system.freeAtoms = mobile       
        self.system.label     = newLabel  
    #=========================================================================
    def Set_QC_Method(self,_parameters,_DEBUG=False):
        """Define quantum chemistry method for QM/MM simulations.
        
        Args:
            _parameters (dict): QM method parameters including method_class,
                Hamiltonian, functional, basis, etc.
            _DEBUG (bool): Export QC system to PDB. Default: False
        """
        if len(self.quantumRegion) > 0: _parameters["region"] = self.quantumRegion
        _parameters["active_system"] = self.system 
        qs = QuantumMethods(_parameters)
        qs.Set_QC_System()
        if not "method_class" in _parameters: _parameters["method_class"] = "SMO"
        if _DEBUG: qs.Export_QC_System()
        newLabel = self.system.label + "QC_system_"
        if "Hamiltonian" in _parameters: newLabel += _parameters["Hamiltonian"] 
        if "functional"  in _parameters: newLabel += _parameters["functional"] 
        self.system.label += newLabel
        self.system = qs.system
        
    #=========================================================================
    def Set_QCMM_Region(self,_pat_list,_centerAtom=None,_radius=None,_DEBUG=False):
        """Define QM region for QM/MM simulations.
        
        Args:
            _pat_list (list): List of atom patterns to include in QM region.
            _centerAtom (str): Alternative: center atom for spherical selection.
            _radius (float): Radius for spherical selection around center atom.
            _DEBUG (bool): Print debug information. Default: False
        """
        if len(_pat_list) > 0:
            for pat in _pat_list:
                _sel =  AtomSelection.FromAtomPattern(self.system,pat)
                self.quantumRegion += _sel

        if _centerAtom:
            atomRef = AtomSelection.FromAtomPattern(self.system,_centerAtom)
            core    = AtomSelection.Within(self.system,atomRef,_radius)
            self.quantumRegion = AtomSelection.ByComponent(self.system,core) 

        self.quantumRegion = list(self.quantumRegion)
       
    #=========================================================================
    def Set_Reaction_crd(self,atoms_rc,_type,_mass_c):
        """Define a reaction coordinate for the system.
        
        Args:
            atoms_rc (list): Atom indices or patterns defining the coordinate.
            _type (str): Type of coordinate (Distance, Angle, Dihedral, etc.)
            _mass_c (bool): Whether to apply mass-weighted constraint.
        """
        _atom_pat = []
        pat = -1
        for atom in atoms_rc:
            try:
                pat= int(atom)
            except:
                pat=AtomSelection.FromAtomPattern(self.system, atom)[0]
            if (pat ==-1): 
                print("Error in selection of atom pattern!!")
            _atom_pat.append(pat )
        
        
        _rc = ReactionCoordinate(_atom_pat,_mass_c,_type)
        _rc.GetRCLabel(self.system)

        self.reactionCoordinates.append(_rc)
        self.rcs +=1 
        

#==============================================================================#
#=====================END OF CLASS FILE========================================#
#==============================================================================#
