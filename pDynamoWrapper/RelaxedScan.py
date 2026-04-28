#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Relaxed Surface Scan Module for Reaction Coordinate Exploration.

This module implements relaxed surface scanning (RSS) methodology for exploring
potential energy surfaces along molecular reaction coordinates. It performs
iterative geometry optimization with fixed distance/angle constraints, enabling
Construction of 1D and 2D energy surfaces and reaction pathways.

Classes:
    SCAN: Main class for managing and executing relaxed scans.

Methods include:
    - 1D/2D relaxed scan execution
    - Adaptive convergence parameter adjustment
    - Support for distance and dihedral reaction coordinates.
"""

#FILE = RelaxedScan.py

#=============================================================
import pymp
import numpy as np 
#-----------------------------
#my library
from .GeometrySearcher import * 
#-----------------------------
#pDynamo library
from pMolecule import *
from pMolecule.QCModel import *
#**************************************************************
class SCAN:
    """Execute relaxed surface scan for reaction coordinate exploration.
    
    This class manages the setup and execution of relaxed scans along 1 or 2
    reaction coordinates. It performs iterative geometry optimization with
    distance/angle constraints to map potential energy surfaces.
    
    Attributes:
        molecule: pDynamo System object.
        baseName (str): Output directory for scan results.
        nDim (int): Number of active reaction coordinates (1 or 2).
        energiesMatrix (ndarray): 2D array storing energy values.
        reactionCoordinate1, reactionCoordinate2 (list): RC values.
        forceC (list): Harmonic restraint force constants.
        nprocs (int): Number of parallel processors.
        adaptative (bool): Enable adaptive convergence adjustment.
    """
    #---------------------------------------------------------------
    def __init__(self,_system,_baseFolder,_optimizer,ADAPTATIVE=False,RESTART=False):
        """Initialize SCAN object.
        
        Args:
            _system (pDynamo System): Molecular system for scanning.
            _baseFolder (str): Output directory for results.
            _optimizer (str): Geometry optimizer algorithm ('ConjugatedGradient', 'SteepestDescent').
            ADAPTATIVE (bool, optional): Enable adaptive convergence. Defaults to False.
            RESTART (bool, optional): Restart from previous scan. Defaults to False.
        """
        self.parameters         = None 
        
        self.baseName           = _baseFolder
        self.molecule           = _system 
        self.nDim               = 0                                  # Number of active reaction coordinates to Scan
        self.reactionCoordinate1= []                                 # array with the first reaction coordinate in angstroms
        self.reactionCoordinate2= []                                 # array with the second reaction coordinate in angstroms
        self.atoms              = []                                 # array of the atomic indices for the reaction coordinates
        self.nprocs             = 1                                  # Maximum virtual threads to use in parallel runs using pymp
        self.energiesMatrix     = None                               # Multidimensional array to store calculated energy values
        self.forceC             = [ 2500.0, 2500.0 ]                 # Force constant for restraint model
        self.forceCRef          = [ self.forceC[0], self.forceC[1] ] # Inital value for the force constant
        self.EnergyRef          = 0.0                               # Float to hold energy reference value for adaptative scheme
        self.nsteps             = [ 1, 1 ]                          # List of integer indicating the number of steps to be taken
        self.maxIt              = 800                               # Maximum number of iterations in goemtry search 
        self.rmsGT              = 0.1                               # Float with root mean square tolerance for geometry optimization 
        self.optmizer           = _optimizer                        # string with optimizer algorithm for geomtry optimization
        self.adaptative         = ADAPTATIVE                        # Boolean indicating if the scan can use the adptative scheme to change the convergence paramters
        self.text               = ""                                # Text container for energy log
        self.saveFormat         = None
        self.trajFolder         = "ScanTraj"
        self.restart            = RESTART
        self.RCs                = []                                 # List to hold the reaction coordinate objects       
        #------------------------------------------------------------------------  
        #set the parameters dict for the geometry search classes
        self.GeoOptPars =   { "maxIterations":self.maxIt  ,
                              "rmsGradient"  :self.rmsGT   }
    #===========================================================================================
    def ChangeDefaultParameters(self,_parameters):
        """Modify default scan parameters.
        
        Updates scanning parameters including force constants, convergence criteria,
        optimization method, and parallelization options.
        
        Args:
            _parameters (dict): Dictionary with parameter updates. Supported keys:
                - 'traj_folder_name': Trajectory folder name.
                - 'rmsGradient': RMS gradient convergence threshold.
                - 'maxIterations': Maximum optimization iterations.
                - 'log_frequency': Logging frequency.
                - 'optimizer': Geometry optimizer type.
                - 'NmaxThreads': Number of parallel threads.
                - 'save_format': Coordinate file format (.dcd, .ptGeo, etc.).
                - 'force_constants': List containing [FC1, FC2] for restraints.
        """
        self.parameters = _parameters
        self.parameters["system_name"]         = self.molecule.label
        self.parameters["initial_coordinates"] = self.molecule.coordinates3 
        #-----------------------------------------------------------
        if "traj_folder_name" in _parameters: self.trajFolder = _parameters["traj_folder_name"]
        if "rmsGradient"      in _parameters: self.GeoOptPars["rmsGradient"]   = _parameters["rmsGradient"]
        if "maxIterations"    in _parameters: self.GeoOptPars["maxIterations"] = _parameters["maxIterations"]
        if "log_frequency"    in _parameters: self.GeoOptPars["log_frequency"] = _parameters["log_frequency"]
        if "optimizer"        in _parameters: self.GeoOptPars["optmizer"]      = _parameters["optmizer"]
        if "NmaxThreads"      in _parameters: self.nprocs                      = _parameters["NmaxThreads"]        
        if "save_format"      in _parameters: self.saveFormat = _parameters["save_format"] 
        if "force_constants"  in _parameters:
            cnt=0
            for fc in self.parameters["force_constants"]:
                self.forceC[cnt] = fc
                cnt +=1
                

    #===========================================================================================
    def ChangeConvergenceParameters(self,_xframe,_yframe):
        """Adaptively adjust convergence parameters based on energy barriers.
        
        Implements adaptive convergence scheme that loosens SCF/geometry convergence
        in regions of high energy to improve computational efficiency while maintaining
        accuracy at low energies.
        
        Args:
            _xframe (int): Index of current X coordinate point.
            _yframe (int): Index of current Y coordinate point.
            
        Note:
            Convergence parameters are adjusted based on energy delta from reference.
            Force constants may also be reduced in high-energy regions.
        """
        if not self.energiesMatrix[_xframe,_yframe] == 0.0:
            delta = self.energiesMatrix[_xframe,_yframe]  
            if delta < 150.0:
               
                self.molecule.qcModel.converger.energyTolerance  = 0.0001
                self.molecule.qcModel.converger.densityTolerance = 3e-08
                self.molecule.qcModel.converger.diisDeviation    = 1e-06
            elif delta >= 150.0:
                self.forceC[0] = self.forceCRef[0] - self.forceCRef[0]*0.40
                self.forceC[1] = self.forceCRef[1] - self.forceCRef[1]*0.40
                self.molecule.qcModel.converger.energyTolerance  = 0.0003
                self.molecule.qcModel.converger.densityTolerance = 3e-08
                self.molecule.qcModel.converger.diisDeviation    = 1e-06
                if delta > 160.0 and delta < 170.0:
                    self.forceC[0] = self.forceCRef[0] - self.forceCRef[0]*0.50
                    self.forceC[1] = self.forceCRef[1] - self.forceCRef[1]*0.50
                    self.molecule.qcModel.converger.energyTolerance  = 0.0006
                    self.molecule.qcModel.converger.densityTolerance = 1e-07
                    self.molecule.qcModel.converger.diisDeviation    = 2e-06
                elif delta > 170.0 and delta <180.0 :
                    self.forceC[0] = self.forceCRef[0] - self.forceCRef[0]*0.50
                    self.forceC[1] = self.forceCRef[1] - self.forceCRef[1]*0.50
                    self.molecule.qcModel.converger.energyTolerance  = 0.001
                    self.molecule.qcModel.converger.densityTolerance = 3e-07
                    self.molecule.qcModel.converger.diisDeviation    = 5e-06
                elif delta > 180.0 and delta < 185.0:
                    self.forceC[0] = self.forceCRef[0] - self.forceCRef[0]*0.70
                    self.forceC[1] = self.forceCRef[1] - self.forceCRef[0]*0.70
                    self.molecule.qcModel.converger.energyTolerance  = 0.0015
                    self.molecule.qcModel.converger.densityTolerance = 1e-06
                    self.molecule.qcModel.converger.diisDeviation    = 1e-05
                elif delta > 185.0 and delta <200.0:
                    self.forceC[0] = self.forceCRef[0] - self.forceCRef[0]*0.70
                    self.forceC[1] = self.forceCRef[1] - self.forceCRef[1]*0.70
                    self.molecule.qcModel.converger.energyTolerance  = 0.003
                    self.molecule.qcModel.converger.densityTolerance = 1e-05
                    self.molecule.qcModel.converger.diisDeviation    = 5e-05
                elif delta > 200.0:
                    self.forceC[0] = self.forceCRef[0] - self.forceCRef[0]*0.80
                    self.forceC[1] = self.forceCRef[1] - self.forceCRef[1]*0.80
                    self.molecule.qcModel.converger.energyTolerance  = 0.003
                    self.molecule.qcModel.converger.densityTolerance = 1e-04
                    self.molecule.qcModel.converger.diisDeviation    = 5e-04
    
    #=============================================================================================
    def SetReactionCoord(self,_RC):
        """Define a reaction coordinate for the relaxed scan.
        
        Registers a reaction coordinate (distance, angle, or dihedral) and extracts
        relevant structural parameters for constraining during the scan.
        
        Args:
            _RC (ReactionCoordinate): Reaction coordinate object specifying atoms,
                constraints, and restraint properties.
                
        Sets:
            nDim: Increments dimension counter.
            atoms: Adds atom indices for constraint.
            DINCREMENT: Distance/angle increment value.
            DMINIMUM: Minimum distance/angle value.
            multipleDistance: Flag for 3-atom distance constraints.
            dihedral: Flag for dihedral constraints.
        """
        #------------------------------------------------------------
        self.RCs.append(_RC)
        self.nDim += 1
        self.atoms.append(_RC.atoms)  

    #===============================================================================================
    def Run1DScan(self,_nsteps):
        """Execute one-dimensional relaxed surface scan.
        
        Performs iterative geometry optimization along a single reaction coordinate,
        constraining the RC to fixed values while allowing all other coordinates to relax.
        
        Args:
            _nsteps (int): Number of steps along the reaction coordinate.
            
        Generates:
            Trajectory folder: Contains optimized structures at each RC point.
            Log file: Energy values and RC coordinates for all steps.
            
        Note:
            Handles both simple and multiple-distance constraints, as well as dihedral angles.
        """
        if not os.path.exists( os.path.join( self.baseName, self.trajFolder +".ptGeo" ) ):
            os.makedirs(  os.path.join( self.baseName,self.trajFolder +".ptGeo"  ) )

        text_line = "{0:>3s} {1:>15s} {2:>15s} {3:>15s}".format('x','RC1','Energy', "Energy(kcal/mol)" )
        self.text += text_line+"\n"
        
        if _nsteps == -1:
            print("Number of steps not provided. Defining steps based on RC parameters.")           
            self.RCs[0].DefineSteps()
            _nsteps = self.RCs[0].nsteps
        else: self.RCs[0].nsteps = _nsteps

        self.energiesMatrix      = pymp.shared.array( (_nsteps), dtype=float ) 
        self.reactionCoordinate1 = pymp.shared.array( (_nsteps), dtype=float )
        
        
        if self.RCs[0].Type == "Dihedral":  
            self.Run1DScanDihedral()
        else:
            if    self.RCs[0].Type == "multipleDistance":
                self.Run1DScanMultipleDistance(_nsteps)
            else: 
                self.Run1DScanSimpleDistance(_nsteps)
        for i in range(_nsteps):
            kcal = self.energiesMatrix[i]/4.184              
            text_line =  "{0:3d} {1:15.8f} {2:15.8f} {3:15.8f}".format( i,self.reactionCoordinate1[i], self.energiesMatrix[i], kcal)
            self.text += text_line+ '\n'
            
    #=================================================================================================
    def Run1DScanSimpleDistance(self,_nsteps):
        """Execute relaxed scan with simple 2-atom distance constraint.
        
        Performs geometry optimization with a fixed distance constraint between
        two atoms. This is the most common type of relaxed surface scan.
        
        Note:
            This method is automatically called by Run1DScan() when appropriate.
            Structures are saved as pickled coordinate objects at each RC point.
        """
        #-------------------------------------------------------------------------
        #Setting some local vars to ease the notation in the pDynamo methods
        #----------------------------------
        atom1 = self.RCs[0].atoms[0]
        atom2 = self.RCs[0].atoms[1]                
        _dminimum   = self.RCs[0].minimumD 
        _dincrement = self.RCs[0].increment 
        #---------------------------------
        restraints = RestraintModel()
        self.molecule.DefineRestraintModel( restraints )                     
        #----------------------------------------------------------------------------------------
        for i in range(_nsteps):       
            distance = _dminimum + ( _dincrement * float(i) )
            #--------------------------------------------------------------------
            rmodel            = RestraintEnergyModel.Harmonic(distance, self.forceC[0])
            restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom1, point2= atom2)
            restraints["RC1"] = restraint            
            #--------------------------------------------------------------------
            
            if i > 0:
                initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}.pkl".format(i-1) )
                if self.restart:
                    if os.path.exists( os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}.pkl".format(i) ) ):
                        initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}.pkl".format(i) )                
                self.molecule.coordinates3 = ImportCoordinates3(initCoordinateFile,log=None)
            relaxRun = GeometrySearcher(self.molecule,self.baseName)
            relaxRun.ChangeDefaultParameters(self.GeoOptPars)
            relaxRun.Minimization(self.optmizer)
            #--------------------------------------------------------------------
            if i == 0:
                self.EnergyRef = en0 = self.molecule.Energy(log=None)
                self.energiesMatrix[i] = 0.0
            else: self.energiesMatrix[i] = self.molecule.Energy(log=None) - en0 
            #--------------------------------------------------------------------
            self.reactionCoordinate1[i] = self.molecule.coordinates3.Distance( atom1 , atom2  )   
            Pickle( os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}.pkl".format(i) ), self.molecule.coordinates3 ) 
        #---------------------------------------
        self.molecule.DefineRestraintModel(None)
    #===================================================================================================
    def Run1DScanMultipleDistance(self,_nsteps):
        '''
        '''
        atom1 = self.RCs[0].atoms[0]
        atom2 = self.RCs[0].atoms[1]
        atom3 = self.RCs[0].atoms[2]
        weight1 = self.RCs[0].weight13
        weight2 = self.RCs[0].weight31 
        _dminimum   = self.RCs[0].minimumD 
        _dincrement = self.RCs[0].increment
        #---------------------------------
        restraints = RestraintModel()
        self.molecule.DefineRestraintModel( restraints )
        #---------------------------------
        for i in range(0,_nsteps):

            distance = _dminimum + ( _dincrement * float(i) )             
            #--------------------------------------------------------------------
            rmodel    = RestraintEnergyModel.Harmonic( distance, self.forceC[0] )
            restraint = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom2, atom1, weight1 ], [ atom2, atom3, weight2 ] ] )
            restraints["RC1"] =  restraint            
            #--------------------------------------------------------------------
            if i > 0:
                initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}.pkl".format(i-1) )
                if self.restart:
                    if os.path.exists( os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}.pkl".format(i) ) ):
                        initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}.pkl".format(i) )                
                self.molecule.coordinates3 = ImportCoordinates3(initCoordinateFile,log=None)
            relaxRun = GeometrySearcher(self.molecule, self.baseName)
            relaxRun.ChangeDefaultParameters(self.GeoOptPars)
            relaxRun.Minimization(self.optmizer)
            #--------------------------------------------------------------------
            if i == 0:
                self.EnergyRef = en0 = self.molecule.Energy(log=None)
                self.energiesMatrix[0] = 0.0
            else: self.energiesMatrix[i] = self.molecule.Energy(log=None) - en0 
            #--------------------------------------------------------------------
            self.reactionCoordinate1[i] = self.molecule.coordinates3.Distance( atom1 , atom2  ) - self.molecule.coordinates3.Distance( atom2, atom3  ) 
            Pickle( os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}.pkl".format(i) ), self.molecule.coordinates3 )
        self.molecule.DefineRestraintModel(None)
    #===================================================================================================
    def Run1DScanDihedral(self,_nsteps):
        '''
        Run scan in dihedral angles.
        '''
        atom1 = self.RCs[0].atoms[0]
        atom2 = self.RCs[0].atoms[1]
        atom3 = self.RCs[0].atoms[2]
        atom4 = self.RCs[0].atoms[3]     

        _dminimum   = self.RCs[0].minimumD 
        _dincrement = self.RCs[0].increment

        #---------------------------------
        restraints = RestraintModel()
        self.molecule.DefineRestraintModel( restraints )
        #---------------------------------
        if _dincrement == 0.0: _dincrement = 360.0/float(self.nsteps[0])
        #----------------------------------------------------------------------------------------
        for i in range(0,_nsteps):
            angle = _dminimum + _dincrement * float(i)
            #--------------------------------------------------------------------
            rmodel    = RestraintEnergyModel.Harmonic( angle, self.forceC[0], period = 360.0 )
            restraint = RestraintDihedral.WithOptions( energyModel = rmodel,
                                                       point1      = atom1 ,
                                                       point2      = atom2 ,
                                                       point3      = atom3 ,
                                                       point4      = atom4 )
            restraints["PHI"] =  restraint            
            #--------------------------------------------------------------------
            if i > 0:
                initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}.pkl".format(i-1) )
                if self.restart:
                    if os.path.exists( os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}.pkl".format(i) ) ):
                        initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}.pkl".format(i) )                
                self.molecule.coordinates3 = ImportCoordinates3(initCoordinateFile,log=None)
            relaxRun = GeometrySearcher(self.molecule, self.baseName  )
            relaxRun.ChangeDefaultParameters(self.GeoOptPars)
            relaxRun.Minimization(self.optmizer)
            #--------------------------------------------------------------------
            if i == 0:
                self.EnergyRef = en0 = self.molecule.Energy(log=None)
                self.energiesMatrix[0] = 0.0
            else: self.energiesMatrix[i]  = self.molecule.Energy(log=None) - en0 
            #--------------------------------------------------------------------
            self.reactionCoordinate1[i] = self.molecule.coordinates3.Dihedral( atom1 , atom2, atom3, atom4 ) 
            Pickle( os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}.pkl".format(i) ), self.molecule.coordinates3 )
        #----------------------------------------
        self.molecule.DefineRestraintModel(None)
    #===================================================================================================
    def Run2DScan(self,_nsteps_x,_nsteps_y):
        '''
        Run two-dimensional relaxed surface scan.
        '''
        if not os.path.exists( os.path.join( self.baseName, self.trajFolder +".ptGeo" ) ):  os.makedirs(  os.path.join( self.baseName,self.trajFolder +".ptGeo"  ) )
        #------------------------------------------------------
        #self.text += "x y RC1 RC2 Energy\n" 
        
        text_line = "{0:>3s} {1:>3s} {2:>15s} {3:>15s} {4:>15s} {5:>15s}".format('x', 'y', 'RC1', 'RC2', 'Energy', "Energy(kcal/mol)" )

        self.text += text_line+"\n"
        #------------------------------------------------------               
        if _nsteps_x == -1:             
            _nsteps_x = self.RCs[0].DefineSteps()    
        else: self.RCs[0].nsteps = _nsteps_x    
        if _nsteps_y == -1: 
            _nsteps_y = self.RCs[1].DefineSteps()
        else: self.RCs[1].nsteps = _nsteps_y
        
        self.RCs[0].Print()
        self.RCs[1].Print()

        print("Running 2D Scan with {} steps in RC1 and {} steps in RC2".format(_nsteps_x,_nsteps_y))
                
        self.energiesMatrix = pymp.shared.array( (_nsteps_x,_nsteps_y), dtype=float ) 
        self.reactionCoordinate1 = pymp.shared.array( (_nsteps_x,_nsteps_y), dtype=float )   
        self.reactionCoordinate2 = pymp.shared.array( (_nsteps_x,_nsteps_y), dtype=float )   
        if _nsteps_x > 0 and _nsteps_y > 0:
            if self.RCs[0].Type == "dihedral": self.Run2DScanDihedral(_nsteps_x,_nsteps_y)
            else:
                if self.RCs[0].Type == "multipleDistance" and self.RCs[1].Type == "multipleDistance": self.Run2DScanMultipleDistance(_nsteps_x,_nsteps_y)
                elif self.RCs[0].Type == "multipleDistance" and self.RCs[1].Type == "Distance": self.Run2DMixedDistance(_nsteps_x,_nsteps_y)
                else: self.Run2DSimpleDistance(_nsteps_x,_nsteps_y)
            #------------------------------------
            for i in range(_nsteps_x):
                for j in range(_nsteps_y):
                    kcal = self.energiesMatrix[i,j]/4.184
                    text_line =  "{0:3d} {1:3d} {2:15.8f} {3:15.8f} {4:15.8f} {5:15.8f}".format( i,j,self.reactionCoordinate1[i,j], self.reactionCoordinate2[i,j], self.energiesMatrix[i,j], kcal)
                    self.text += text_line+ '\n'
    #=============================================================================
    def Run2DSimpleDistance(self, X, Y ):
        '''
        Run two-dimensional simple distance relaxed surface scan
        '''
        atom1 = self.RCs[0].atoms[0]
        atom2 = self.RCs[0].atoms[1]
        atom3 = self.RCs[1].atoms[0]
        atom4 = self.RCs[1].atoms[1]   

        _dminimum_A   = self.RCs[0].minimumD 
        _dincrement_A = self.RCs[0].increment
        _dminimum_B   = self.RCs[1].minimumD 
        _dincrement_B = self.RCs[1].increment

        restraints = RestraintModel( )
        self.molecule.DefineRestraintModel( restraints )

        self.reactionCoordinate1[ 0,0 ] = self.molecule.coordinates3.Distance( atom1, atom2 ) 
        self.reactionCoordinate2[ 0,0 ] = self.molecule.coordinates3.Distance( atom3, atom4 ) 

        rmodel     =  RestraintEnergyModel.Harmonic( _dminimum_A, self.forceC[0] )
        restraint  =  RestraintDistance.WithOptions( energyModel = rmodel,  point1=atom1, point2=atom2  )
        restraints["RC1"] = restraint                
        #----------------------------------------------------------------------------------------------                
        rmodel      = RestraintEnergyModel.Harmonic( _dminimum_B, self.forceC[1] )
        restraint   = RestraintDistance.WithOptions( energyModel = rmodel, point1=atom3, point2=atom4 )                    
        restraints["RC2"] = restraint 
        #----------------------------------------------------------------------------------------------
        coordinateFile = os.path.join( self.baseName ,self.trajFolder+".ptGeo","frame{}_{}.pkl".format( 0, 0 ) )
        relaxRun = GeometrySearcher( self.molecule, self.baseName )
        relaxRun.ChangeDefaultParameters( self.GeoOptPars )
        relaxRun.Minimization(self.optmizer)
        #-----------------------------------------------------------------------------------------------
        self.EnergyRef = self.en0 = self.molecule.Energy(log=None)         
        Pickle( coordinateFile, self.molecule.coordinates3 )

        for i in range ( 1, X ):  
            #---------------------------------------------------------------------------------------------
            if self.adaptative: self.ChangeConvergenceParameters(i-1,0)
            #---------------------------------------------------------------------------------------------             
            distance_1 = _dminimum_A + ( _dincrement_A * float(i) ) 
            rmodel     =  RestraintEnergyModel.Harmonic( distance_1, self.forceC[0] )
            restraint  =  RestraintDistance.WithOptions( energyModel = rmodel,  point1=atom1, point2=atom2  )
            restraints["RC1"] = restraint                
            #----------------------------------------------------------------------------------------------                
            distance_2  = _dminimum_B                
            rmodel      = RestraintEnergyModel.Harmonic( distance_2, self.forceC[1] )
            restraint   = RestraintDistance.WithOptions( energyModel = rmodel, point1=atom3, point2=atom4 )                    
            restraints["RC2"] = restraint  
            #---------------------------------------------------------------------------------------------- 
            initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format( 0 , 0) ) 
            if i > 0:
                initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format(i-1,0) )
                if self.restart:
                    if os.path.exists( os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format(i,0) ) ):
                        initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format(i,0) )          
            self.molecule.coordinates3 = ImportCoordinates3( initCoordinateFile, log=None )               
            #----------------------------------------------------------------------------------------------
            relaxRun = GeometrySearcher( self.molecule, self.baseName )
            relaxRun.ChangeDefaultParameters(self.GeoOptPars)
            relaxRun.Minimization(self.optmizer)
            #----------------------------------------------------------------------------------------------
            self.energiesMatrix[ i,0 ]      = self.molecule.Energy(log=None) - self.en0
            self.reactionCoordinate1[ i,0 ] = self.molecule.coordinates3.Distance( atom1, atom2 ) 
            self.reactionCoordinate2[ i,0 ] = self.molecule.coordinates3.Distance( atom3, atom4 )   
            #-----------------------------------------------------------------------------------
            coordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format( i, 0 ) )
            Pickle( coordinateFile, self.molecule.coordinates3 )                
        #-------------------------------------------------------------------------------------------
        with pymp.Parallel(self.nprocs) as p:
            #Pergomr the calculations for the rest of the grid
            for i in p.range ( 0, X ):
                #----------------------------------------------------------------------------------------------
                distance_1 = _dminimum_A + ( _dincrement_A  * float(i) )
                rmodel     = RestraintEnergyModel.Harmonic( distance_1, self.forceC[0] )
                restraint  = RestraintDistance.WithOptions(energyModel =rmodel, point1=atom1, point2=atom2  )
                restraints["RC1"] = restraint
                #----------------------------------------------------------------------------------------------
                for j in range( 1, Y ):
                    distance_2  = _dminimum_B + ( _dincrement_B * float(j) )
                    rmodel      = RestraintEnergyModel.Harmonic( distance_2, self.forceC[1] )
                    restraint   = RestraintDistance.WithOptions(energyModel = rmodel, point1=atom3, point2=atom4  )
                    restraints["RC2"] = restraint                    
                    #------------------------------------------------------------------------------------------
                    initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo" , "frame{}_{}.pkl".format( i, j-1 ) )   
                    if self.restart:
                        if os.path.exists( os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format(i,j) ) ):
                            initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format(i,j) )                 
                    self.molecule.coordinates3 = ImportCoordinates3( initCoordinateFile, log=None )             
                    #----------------------------------------------------------------------------------------------
                    if self.adaptative: self.ChangeConvergenceParameters(i,j-1)
                    coordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format( i, j ) )
                    relaxRun = GeometrySearcher( self.molecule, self.baseName  )
                    relaxRun.ChangeDefaultParameters( self.GeoOptPars )
                    relaxRun.Minimization(self.optmizer) 
                    #----------------------------------------------------------------------------------------------
                    self.energiesMatrix[ i,j ] = self.molecule.Energy(log=None) - self.en0
                    self.reactionCoordinate1[ i,j ] = self.molecule.coordinates3.Distance( atom1, atom2 ) 
                    self.reactionCoordinate2[ i,j ] = self.molecule.coordinates3.Distance( atom3, atom4 )
                    #-----------------------------------------------------------------------------------
                    Pickle( coordinateFile, self.molecule.coordinates3 )
                    #-----------------------------------------------------------------------------------
        self.molecule.DefineRestraintModel(None)
    #=======================================================   
    def Run2DMixedDistance(self,X, Y ):
        '''
        Run two-dimensional simple distance relaxed surface scan
        '''
        atom1 = self.RCs[0].atoms[0]
        atom2 = self.RCs[0].atoms[1]
        atom3 = self.RCs[0].atoms[2]   
        atom4 = self.RCs[1].atoms[0]
        atom5 = self.RCs[1].atoms[1]

        weight1 = self.RCs[0].weight13
        weight2 = self.RCs[0].weight31

        _dminimum_A   = self.RCs[0].minimumD 
        _dincrement_A = self.RCs[0].increment
        _dminimum_B   = self.RCs[1].minimumD 
        _dincrement_B = self.RCs[1].increment
      
        
        restraints = RestraintModel( )
        self.molecule.DefineRestraintModel( restraints )
        
        self.reactionCoordinate1[ 0,0 ] = self.molecule.coordinates3.Distance( atom1, atom2 ) - self.molecule.coordinates3.Distance( atom3, atom2 )
        self.reactionCoordinate2[ 0,0 ] = self.molecule.coordinates3.Distance( atom4, atom5 )
       
        distance_1 = _dminimum_A 
        rmodel     = RestraintEnergyModel.Harmonic( distance_1, self.forceC[0] )
        restraint  = RestraintMultipleDistance.WithOptions(energyModel = rmodel, distances = [ [ atom2, atom1, weight1 ],[ atom2, atom3, weight2 ] ] )
        restraints["RC1"] = restraint                
        #--------------------------------------------------------------------------------
        distance_2  = _dminimum_B
        rmodel      = RestraintEnergyModel.Harmonic( distance_2, self.forceC[1] )
        restraint   = RestraintDistance.WithOptions( energyModel = rmodel, point1=atom4, point2=atom5 )                
        restraints["RC2"] = restraint 
        #--------------------------------------------------------------------------------
        coordinateFile = os.path.join( self.baseName ,self.trajFolder+".ptGeo","frame{}_{}.pkl".format( 0, 0 ) )
        relaxRun = GeometrySearcher( self.molecule, self.baseName )
        relaxRun.ChangeDefaultParameters( self.GeoOptPars )
        relaxRun.Minimization(self.optmizer)
        #-----------------------------------------------------------------------------------------------
        self.EnergyRef = self.en0 = self.molecule.Energy(log=None)         
        Pickle( coordinateFile, self.molecule.coordinates3 )

        for i in range ( 1, X ):
            if self.adaptative:
                try: self.ChangeConvergenceParameters(i-1,0) 
                except: pass
            distance_1 = _dminimum_A + _dincrement_A * float(i)
            rmodel     = RestraintEnergyModel.Harmonic( distance_1, self.forceC[0] )
            restraint  = RestraintMultipleDistance.WithOptions(energyModel = rmodel, distances = [ [ atom2, atom1, weight1 ],[ atom2, atom3, weight2 ] ] )
            restraints["RC1"] = restraint                
            #--------------------------------------------------------------------------------
            distance_2  = _dminimum_B
            rmodel      = RestraintEnergyModel.Harmonic( distance_2, self.forceC[1] )
            restraint   = RestraintDistance.WithOptions( energyModel = rmodel, point1=atom4, point2=atom5 )                
            restraints["RC2"] = restraint  
            #---------------------------------------------------------------------------------                    
            initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format(i-1,0) ) 
            if self.restart:
                if os.path.exists( os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format(i,0) ) ):
                    initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format(i,0) )                 

            coordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format( i, 0 ) )
            self.molecule.coordinates3 = ImportCoordinates3( initCoordinateFile )                
            #--------------------------------------------------------------------------------
            relaxRun = GeometrySearcher( self.molecule, self.baseName )
            relaxRun.ChangeDefaultParameters( self.GeoOptPars )
            relaxRun.Minimization(self.optmizer)
            #-------------------------------------------------------------------------------- 
            self.energiesMatrix[ i,0 ]      = self.molecule.Energy(log=None) - self.en0
            self.reactionCoordinate1[ i,0 ] = self.molecule.coordinates3.Distance( atom1, atom2 ) - self.molecule.coordinates3.Distance( atom3, atom2 )
            self.reactionCoordinate2[ i,0 ] = self.molecule.coordinates3.Distance( atom4, atom5 )              
            #-----------------------------------------------------------------------------------
            Pickle( coordinateFile, self.molecule.coordinates3 ) 
        #---------------------------------------------------------------------------------------------
        with pymp.Parallel(self.nprocs) as p:
            #Pergomr the calculations for the rest of the grid
            for i in p.range ( 0, X ):
                distance_1 = _dminimum_A + _dincrement_A * float(i) 
                rmodel     =  RestraintEnergyModel.Harmonic( distance_1, self.forceC[0] )
                restraint  =  RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom2, atom1, weight1 ],[ atom2, atom3, weight2 ] ] )
                restraints["RC1"] = restraint
                #-----------------------------------------------------------------------------------
                for j in range( 1, Y ):
                    distance_2  = _dminimum_B + ( _dincrement_B * float(j) )
                    rmodel      = RestraintEnergyModel.Harmonic( distance_2, self.forceC[1] )
                    restraint   = RestraintDistance.WithOptions( energyModel = rmodel, point1=atom4, point2=atom5 )
                    restraints["RC2"] = restraint  
                    #-----------------------------------------------------------------------------------
                    initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo" , "frame{}_{}.pkl".format( i, j-1 ) )
                    if self.restart:
                        if os.path.exists( os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format(i,j) ) ):
                            initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format(i,j) )             
                    #-----------------------------------------------------------------------------------
                    self.molecule.coordinates3 = ImportCoordinates3( initCoordinateFile, log=None )             
                    #-----------------------------------------------------------------------------------
                    if self.adaptative:
                        print(self.adaptative)
                        try: self.ChangeConvergenceParameters(i,j-1)
                        except:pass
                    #-----------------------------------------------------------------------------------
                    coordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format( i, j ) )
                    relaxRun = GeometrySearcher( self.molecule, self.baseName  )
                    relaxRun.ChangeDefaultParameters( self.GeoOptPars )
                    relaxRun.Minimization(self.optmizer)
                    #-----------------------------------------------------------------------------------
                    self.energiesMatrix[ i,j ]      = self.molecule.Energy(log=None) - self.en0
                    self.reactionCoordinate1[ i,j ] = self.molecule.coordinates3.Distance( atom1, atom2 ) - self.molecule.coordinates3.Distance( atom3, atom2 )
                    self.reactionCoordinate2[ i,j ] = self.molecule.coordinates3.Distance( atom4, atom5 )
                    #-----------------------------------------------------------------------------------
                    Pickle( coordinateFile, self.molecule.coordinates3 )
                    #...................................................
        self.molecule.DefineRestraintModel(None)
    
    #===========================================================
    def Run2DScanMultipleDistance(self, X, Y ):
        '''
        Run two-dimensional simple distance relaxed surface scan
        '''
        atom1 = self.RCs[0].atoms[0]
        atom2 = self.RCs[0].atoms[1]
        atom3 = self.RCs[0].atoms[2]
        atom4 = self.RCs[1].atoms[0]
        atom5 = self.RCs[1].atoms[1]
        atom6 = self.RCs[1].atoms[2]
        weight1 = self.RCs[0].weight13
        weight2 = self.RCs[0].weight31
        weight3 = self.RCs[1].weight13
        weight4 = self.RCs[1].weight31

        _dminimum_A   = self.RCs[0].minimumD 
        _dincrement_A = self.RCs[0].increment
        _dminimum_B   = self.RCs[1].minimumD 
        _dincrement_B = self.RCs[1].increment

        restraints = RestraintModel( )
        self.molecule.DefineRestraintModel( restraints )
        self.reactionCoordinate1[ 0,0 ] = self.molecule.coordinates3.Distance( atom1, atom2 ) - self.molecule.coordinates3.Distance( atom3, atom2 )
        self.reactionCoordinate2[ 0,0 ] = self.molecule.coordinates3.Distance( atom4, atom5 ) - self.molecule.coordinates3.Distance( atom6, atom5 )
        #-------------------------------------------------------------------------------------
        distance_1 = _dminimum_A
        rmodel     = RestraintEnergyModel.Harmonic( distance_1, self.forceC[0] )
        restraint  = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = [ [ atom2, atom1, weight1 ] , [ atom2, atom3, weight2 ] ] )
        restraints["RC1"] = restraint
        #---- ----------------------------------------------------------------------------        
        distance_2  = _dminimum_A
        rmodel      = RestraintEnergyModel.Harmonic( distance_2, self.forceC[1] )
        restraint   = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = [ [ atom5, atom4, weight3 ],[ atom5, atom6, weight4 ] ] )
        restraints["RC2"] = restraint  
        #---- ----------------------------------------------------------------------------   
        coordinateFile = os.path.join( self.baseName ,self.trajFolder+".ptGeo","frame{}_{}.pkl".format( 0, 0 ) )
        relaxRun = GeometrySearcher( self.molecule, self.baseName )
        relaxRun.ChangeDefaultParameters( self.GeoOptPars )
        relaxRun.Minimization(self.optmizer)
        #-----------------------------------------------------------------------------------------------
        self.EnergyRef = self.en0 = self.molecule.Energy(log=None)         
        Pickle( coordinateFile, self.molecule.coordinates3 )     

        for i in range ( 1, X ):  
            #.---- ----------------------------------------------------------------------------            
            distance_1 = _dminimum_A + _dincrement_A * float(i) 
            rmodel     = RestraintEnergyModel.Harmonic( distance_1, self.forceC[0] )
            restraint  = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = [ [ atom2, atom1, weight1 ] , [ atom2, atom3, weight2 ] ] )
            restraints["RC1"] = restraint
            #---------------------------------------------------------------------------------        
            distance_2  = _dminimum_B
            rmodel      = RestraintEnergyModel.Harmonic( distance_2, self.forceC[1] )
            restraint   = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = [ [ atom5, atom4, weight3 ],[ atom5, atom6, weight4 ] ] )
            restraints["RC2"] = restraint  
            #---------------------------------------------------------------------------------
            initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format(i-1,0) ) 
            if self.restart:
                if os.path.exists( os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format(i,0) ) ):
                    initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format(i,0) )   
            self.molecule.coordinates3 = ImportCoordinates3( initCoordinateFile, log = None )             
            #-----------------------------------------------------------------------------------
            if self.adaptative: self.ChangeConvergenceParameters(i-1,0)
            #-----------------------------------------------------------------------------------
            relaxRun = GeometrySearcher( self.molecule, self.baseName )
            relaxRun.ChangeDefaultParameters(self.GeoOptPars)
            relaxRun.Minimization(self.optmizer)
            #-----------------------------------------------------------------------------------
            self.energiesMatrix[ i,0 ]      = self.molecule.Energy(log=None) - self.en0
            self.reactionCoordinate1[ i,0 ] = self.molecule.coordinates3.Distance( atom1, atom2 ) - self.molecule.coordinates3.Distance( atom3, atom2 )
            self.reactionCoordinate2[ i,0 ] = self.molecule.coordinates3.Distance( atom4, atom5 ) - self.molecule.coordinates3.Distance( atom6, atom5 )
            #-----------------------------------------------------------------------------------
            coordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format( i, 0 ) )                       
            Pickle( coordinateFile, self.molecule.coordinates3 )        
        #........................................................................................
        with pymp.Parallel(self.nprocs) as p:
            for i in p.range ( 0, X ):
                distance_1  =  _dminimum_A + _dincrement_A * float(i) 
                rmodel      = RestraintEnergyModel.Harmonic( distance_1, self.forceC[0] )
                restraint   = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom2, atom1, weight1 ],[ atom2, atom3, weight2 ] ] )
                restraints["RC1"] = restraint                       
                #---------------------------------------------------------------------------------
                for j in range( 1, Y ):
                    distance_2  =  _dminimum_B + _dincrement_B * float(j) 
                    rmodel      = RestraintEnergyModel.Harmonic( distance_2, self.forceC[1] )
                    restraint   = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = [ [ atom5, atom4, weight3 ],[ atom5, atom6, weight4 ] ] )
                    restraints["RC2"] = restraint  
                    #---------------------------------------------------------------------------------
                    initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo" , "frame{}_{}.pkl".format( i, j-1 ) )  
                    if self.restart:
                        if os.path.exists( os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format(i,j) ) ):
                            initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format(i,j) )   
                    #----------------------------------------------------------------------------              
                    self.molecule.coordinates3 = ImportCoordinates3( initCoordinateFile, log = None )             
                    if self.adaptative: self.ChangeConvergenceParameters(i,j-1)
                    #----------------------------------------------------------------------------
                    relaxRun = GeometrySearcher( self.molecule, self.baseName  )
                    relaxRun.ChangeDefaultParameters(self.GeoOptPars)
                    relaxRun.Minimization( self.optmizer )
                    #----------------------------------------------------------------------------
                    self.energiesMatrix[ i,j ]      = self.molecule.Energy(log=None) - self.en0
                    self.reactionCoordinate1[ i,j ] = self.molecule.coordinates3.Distance( atom1, atom2 ) - self.molecule.coordinates3.Distance( atom3, atom2 )
                    self.reactionCoordinate2[ i,j ] = self.molecule.coordinates3.Distance( atom4, atom5 ) - self.molecule.coordinates3.Distance( atom6, atom5 )
                    #-----------------------------------------------------------------------------------
                    coordinateFile = os.path.join( self.baseName, self.trajFolder+".ptGeo", "frame"+str(i)+"_"+str(j)+".pkl" )
                    Pickle( coordinateFile, self.molecule.coordinates3 )                    
        #--------------------------------------                
        self.molecule.DefineRestraintModel(None)
    #=======================================================================================
    def Run2DScanDihedral(self, X, Y):
        '''
        Run two-dimensional dihedral relaxed surface scan
        '''
        atom1 = self.RCs[0].atoms[0]
        atom2 = self.RCs[0].atoms[1]
        atom3 = self.RCs[0].atoms[2]
        atom4 = self.RCs[0].atoms[3]
        atom5 = self.RCs[1].atoms[0]
        atom6 = self.RCs[1].atoms[1]
        atom7 = self.RCs[1].atoms[2]
        atom8 = self.RCs[1].atoms[3]

        _dincrement_A = self.RCs[0].increment
        _dincrement_B = self.RCs[1].increment
        _dminimum_A   = self.RCs[0].minimumD
        _dminimum_B   = self.RCs[1].minimumD

        restraints = RestraintModel( )
        self.molecule.DefineRestraintModel( restraints )

        self.reactionCoordinate1[ 0,0 ] = self.molecule.coordinates3.Dihedral( atom1, atom2, atom3, atom4 ) 
        self.reactionCoordinate2[ 0,0 ] = self.molecule.coordinates3.Dihedral( atom5, atom6, atom7, atom8 )
        if _dincrement_A == 0.0: _dincrement_A = 360.0/float(X)
        if _dincrement_B == 0.0: _dincrement_B = 360.0/float(Y)
        
        angle_1    = _dminimum_A
        rmodel     = RestraintEnergyModel.Harmonic( angle_1, self.forceC[0], period = 360.0 )
        restraint  = RestraintDihedral.WithOptions( energyModel = rmodel, 
                                                    point1 = atom1      ,
                                                    point2 = atom2      ,
                                                    point3 = atom3      ,
                                                    point4 = atom4      )
        restraints["RC1"] = restraint
        #---- ----------------------------------------------------------------------------        
        angle_2     = _dminimum_B
        rmodel      = RestraintEnergyModel.Harmonic( angle_2, self.forceC[1], period = 360.0 )
        restraint   = RestraintDihedral.WithOptions( energyModel = rmodel, 
                                                     point1 = atom1      ,
                                                     point2 = atom2      ,
                                                     point3 = atom3      ,
                                                     point4 = atom4      )
        restraints["RC2"] = restraint  
        #---- ----------------------------------------------------------------------------   
        coordinateFile = os.path.join( self.baseName ,self.trajFolder+".ptGeo","frame{}_{}.pkl".format( 0, 0 ) )
        relaxRun = GeometrySearcher( self.molecule, self.baseName )
        relaxRun.ChangeDefaultParameters( self.GeoOptPars )
        relaxRun.Minimization(self.optmizer)
        self.en0 = self.molecule.Energy(log=None)
        Pickle( coordinateFile, self.molecule.coordinates3 ) 
        #-------------------------------------------------------------------------------------
        for i in range ( 1, X ):  
        #.--------------------------------------------------------------------------------            
            angle_1    = _dminimum_A + float(i)*_dincrement_A 
            rmodel     = RestraintEnergyModel.Harmonic( angle_1, self.forceC[0], period = 360.0 )
            restraint  = RestraintDihedral.WithOptions( energyModel = rmodel,
                                                        point1 = atom1      ,
                                                        point2 = atom2      ,
                                                        point3 = atom3      ,
                                                        point4 = atom4      )
            restraints["PHI"] = restraint
            #--------------------------------------------------------------------------------        
            angle_2     = _dminimum_B
            rmodel      = RestraintEnergyModel.Harmonic( angle_2, self.forceC[1], period = 360.0 )
            restraint   = RestraintDihedral.WithOptions( energyModel = rmodel, 
                                                         point1 = atom5      ,
                                                         point2 = atom6      ,
                                                         point3 = atom7      ,
                                                         point4 = atom8      )
            restraints["PSI"] = restraint  
            #---------------------------------------------------------------------------------
            initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format(0,0) ) 
            self.molecule.coordinates3 = ImportCoordinates3( initCoordinateFile, log = None )                          
            #-----------------------------------------------------------------------------------
            relaxRun = GeometrySearcher( self.molecule, self.baseName )
            relaxRun.ChangeDefaultParameters(self.GeoOptPars)
            relaxRun.Minimization(self.optmizer)
            #-----------------------------------------------------------------------------------
            self.energiesMatrix[ i,0 ]      = self.molecule.Energy(log=None) - self.en0
            self.reactionCoordinate1[ i,0 ] = self.molecule.coordinates3.Dihedral( atom1, atom2, atom3, atom4 )
            self.reactionCoordinate2[ i,0 ] = self.molecule.coordinates3.Dihedral( atom5, atom6, atom7, atom8 )
            #-----------------------------------------------------------------------------------
            coordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo", "frame{}_{}.pkl".format( i, 0 ) )                       
            Pickle( coordinateFile, self.molecule.coordinates3 )        
        #........................................................................................
        with pymp.Parallel(self.nprocs) as p:
            for i in p.range ( 0, X ):
                angle_1     = _dminimum_A + float(i)*_dincrement_A 
                rmodel      = RestraintEnergyModel.Harmonic( angle_1, self.forceC[0], period = 360.0  )
                restraint   = RestraintDihedral.WithOptions( energyModel = rmodel, 
                                                             point1 = atom1      ,
                                                             point2 = atom2      ,
                                                             point3 = atom3      ,
                                                             point4 = atom4      )
                restraints["PHI"] = restraint                       
                #---------------------------------------------------------------------------------
                for j in range( 1, Y ):
                    angle_2  =  _dminimum_B + float(j)*_dincrement_B
                    rmodel      = RestraintEnergyModel.Harmonic( angle_2, self.forceC[1], period = 360.0  )
                    restraint   = RestraintDihedral.WithOptions( energyModel = rmodel, 
                                                                 point1 = atom5      ,
                                                                 point2 = atom6      ,
                                                                 point3 = atom7      ,
                                                                 point4 = atom8      )
                    restraints["PSI"] = restraint  
                    #---------------------------------------------------------------------------------
                    initCoordinateFile = os.path.join( self.baseName,self.trajFolder+".ptGeo" , "frame{}_{}.pkl".format( i, j-1 ) )  
                    #----------------------------------------------------------------------------              
                    self.molecule.coordinates3 = ImportCoordinates3( initCoordinateFile, log = None )      
                    #----------------------------------------------------------------------------
                    relaxRun = GeometrySearcher( self.molecule, self.baseName  )
                    relaxRun.ChangeDefaultParameters(self.GeoOptPars)
                    relaxRun.Minimization( self.optmizer )
                    #----------------------------------------------------------------------------
                    self.energiesMatrix[ i,j ]      = self.molecule.Energy(log=None) - self.en0
                    self.reactionCoordinate1[ i,j ] = self.molecule.coordinates3.Dihedral( atom1, atom2, atom3, atom4 ) 
                    self.reactionCoordinate2[ i,j ] = self.molecule.coordinates3.Dihedral( atom5, atom6, atom7, atom8 )
                    #-----------------------------------------------------------------------------------
                    coordinateFile = os.path.join( self.baseName, self.trajFolder+".ptGeo", "frame"+str(i)+"_"+str(j)+".pkl" )
                    Pickle( coordinateFile, self.molecule.coordinates3 )                    
        #--------------------------------------                
        self.molecule.DefineRestraintModel(None)

    #=======================================================================================
    def Finalize(self):
        '''
        Writing logs, making plots and saving trajectories
        '''       
        if self.nDim == 1:
            #..................................................
            if not self.saveFormat == None: 
                trajName = os.path.join( self.baseName, self.trajFolder+self.saveFormat )
                trajpath = os.path.join( self.baseName, self.trajFolder+".ptGeo" )
                Duplicate( trajpath, trajName, self.molecule )           
                

        textLog = open( os.path.join(self.baseName,self.trajFolder+".log"), "w" ) 
        textLog.write(self.text)
        textLog.close() 

        pymol_text = "preset.publication(selection='all')\n"
        pymol_text+= "set sticks\n"
        pymol_text+= "set label_size, 20\n"
        pymol_text+= "set sphere_scale, 0.2\n"
        pymol_text+= "set bg_rgb, white\n" 
        pymol_text+= "set stick_radius, 0.18\n"
        pymol_text+= "load {}".format( os.path.join( self.baseName, self.trajFolder+".ptGeo", "frame0.pdb" ) )
        pymol_text+= "\nload_traj {}, ".format( os.path.join( self.baseName, self.trajFolder+self.saveFormat ) )
        pymol_text+= "frame0, 1, start=1, stop=-1, interval=1"
        
        pymols_file = open( os.path.join(self.baseName,"traj1d.pym"), "w") 
        pymols_file.write(pymol_text)
        pymols_file.close()

        return( os.path.join(self.baseName,self.trajFolder+".log") )
        

    #========================================================================================
    def Print(self):
        '''
        Printing relaxed scan parameters
        '''
        pass


#==============================================================================#
#=====================END OF CLASS FILE========================================#
#==============================================================================#
