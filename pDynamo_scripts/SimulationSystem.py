#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = CoreInterface.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

import os, glob, sys
#--------------------------------------------------------------
#Loading own libraries
from commonFunctions import *
from Simulation import Simulation
from Analysis import *
from QuantumMethods import *
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
from ReactionCoordinate import ReactionCoordinate
#==========================================================================

#**************************************************************************
class SimulationSystem:
    '''
    Class to Wrapper information around Sysyem class from pDynamo3 Libraries
    '''  
    #.-------------------------------------------------------------------------
    def __init__(self,_label="No specified"):
        '''
        Class constructor
        Parameters:
        '''        
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
        '''
        Initialize project object from a pdynamo PKL file.
        '''
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
        '''
        Initialize project from force field topology and coordinate files.
        '''
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
        '''
        '''
        self = selfClass()
        parameters   = GromacsParameterFileReader.PathToParameters ( _topologyFile )
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
        '''
        Initialize project from coordinates 
        '''
        self = selfClass()        
        self.system   = ImportSystem(_coordinateFile)
        _name         = os.path.basename(_coordinateFile)
        self.baseName = _name[:-4]
        return(self) 
    #===================================================================================
    @classmethod
    def Protein_From_Coordinates(selfClass,_coordinateFile,_modelNumber=1):
        '''
        Initialize project from coordinate file with OPLS general force field
        '''
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
        '''
        Calculates single point energy.
        '''
        self.refEnergy = self.system.Energy(doGradients = False)
        self.system.Summary()
        return self.refEnergy    
    #====================================================================================
    def Spherical_Pruning(self,_centerAtom,_radius,_DEBUG=False):
        '''
        Perform a spherical pruning from a certain atom center coordinates.
        Parameters:
            _centerAtom:
            _radius    :
        '''
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
        '''
        Set the list of atoms to keep with the positions fixed through the next simulations
        Parameters:
            _centerAtom:
            _radius    :
        '''
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
        '''
        '''
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
        '''
            lig = AtomSelection.FromAtomPattern(proj.system,"*:LIG.248:*")

        '''
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
        '''
        '''
        _atom_pat = []
        for atom in atoms_rc:            
            _atom_pat.append( AtomSelection.FromAtomPattern(self.system, atom)[0] )
        
        
        _rc = ReactionCoordinate(_atom_pat,_mass_c,_type)
        _rc.GetRCLabel(self.system)

        self.reactionCoordinates.append(_rc)
        self.rcs +=1 
        

#================================================================================

