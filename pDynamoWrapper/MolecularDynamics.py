#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MolecularDynamics module - MD simulation control and orchestration.

This module provides the MD class for setting up and executing molecular
dynamics simulations using pDynamo. Supports multiple integrators (Velocity Verlet,
Leap Frog, Langevin), temperature and pressure control, trajectory management,
and soft constraints for restricted simulations.

Authors: Igor Barden Grillo, contributors
"""

#FILE = MolecularDynamics.py


#---------------------------------------
#importing libraries
import os
#import sys
#----------------------------------------
# pDynamo
from pBabel                    import *                                     
from pCore                     import *                                     
from pMolecule                 import *                              
from pMolecule.MMModel         import *
from pMolecule.NBModel         import *                                     
from pMolecule.QCModel         import *
from pScientific               import *                                     
from pScientific.Arrays        import *                                     
from pScientific.Geometry3     import *                                     
from pScientific.RandomNumbers import *                                     
from pScientific.Statistics    import *
from pScientific.Symmetry      import *                                     
from pSimulation               import *

#---------------------------------------

#**************************************************************************
class MD:
    """Manager for molecular dynamics simulations.
    
    Handles setup and execution of MD simulations with various integrators
    (Velocity Verlet, Leap Frog, Langevin). Manages heating, equilibration,
    and production stages with configurable temperature, pressure control,
    and soft constraints for restricted simulations.
    
    Attributes:
        molecule (System): Molecular system to simulate.
        baseName (str): Base directory for output.
        trajectory (Trajectory): Cartesian coordinate trajectory.
        trajectorySoft (Trajectory): Soft-constraint trajectory.
        algorithm (str): Integration algorithm.
        Nsteps (int): Total number of simulation steps.
        timeStep (float): Integration time step (ps).
        temperature (float): Target temperature (K).
        pressureControl (bool): Apply pressure coupling.
        samplingFactor (int): Save frequency.
    """
    #.---------------------------------------
    def __init__(self,_system,_baseFolder,_parameters):
        """Initialize MD simulation manager.
        
        Args:
            _system (System): Molecular system for MD.
            _baseFolder (str): Directory for trajectory output.
            _parameters (dict): MD configuration with keys like:
                - trajectory_name, MD_method, temperature, timeStep
                - pressure, pressure_coupling, seed, log_frequency
                - coll_freq, temperature_scale_option, etc.
        """        
        #Important parameters that are recurrently wanted to be change by the user
        self.molecule               = _system
        self.baseName               = _baseFolder     
        self.trajName               = _parameters["trajectory_name"]
        self.trajectoryNameSoft     = os.path.join(_baseFolder,self.trajName+".ptRes")
        self.trajectoryNameCurr     = os.path.join(_baseFolder,self.trajName+".ptGeo")
        self.trajectory             = None
        self.trajectorySoft         = None
        self.algorithm              = _parameters["MD_method"]
        self.saveFormat             = None # binary file format to save the trajectory
        self.Nsteps                 = _parameters["production_nsteps"] 
        self.timeStep               = _parameters["timeStep"]
        self.temperature            = _parameters["temperature"]
        self.pressureControl        = False
        self.samplingFactor         = 100
        self.logFreq                = _parameters["log_frequency"]
        self.seed                   = _parameters["seed"]
        self.softConstraint         = False # boolean flags signlizing whether the system has soft constraints or not              
        #Default constants less acessible by the users
        self.collFreq               = _parameters["coll_freq"]
        self.pressureCoupling       = _parameters["pressure_coupling"]    
        self.pressure               = _parameters["pressure"]  
        self.temperatureScaleOption = _parameters["temperature_scale_option"]
        self.temperatureScaleFreq   = 100
        self.startTemperature       = _parameters["start_temperature"]
        self.DEBUG                  = False
        #self.NDIM                   = _parameters["NmaxThreads"]
        #Setting parameters based on information that we collected on the instance construction
        self.RNG                    = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( self.seed ) )
        if not os.path.exists(_baseFolder): os.makedirs(_baseFolder)

        if _parameters["pressure_coupling"] == "True": self.pressureControl = True 

    #=============================================================================================    
    def HeatingSystem(self,_nsteps,_samplingFactor):
        """Run heating simulation to bring system to target temperature.
        
        Performs Velocity Verlet MD with linear temperature scaling
        from initial temperature to target temperature.
        
        Args:
            _nsteps (int): Number of MD steps for heating.
            _samplingFactor (int): Frequency for saving trajectory frames.
        """
        self.nsteps             = _nsteps
        self.trajectoryNameCurr = os.path.join(self.baseName,self.trajName+"heating.ptGeo")  
        self.trajectory         = ExportTrajectory( self.trajectoryNameCurr, self.molecule,log=None )
        self.samplingFactor     = _samplingFactor
        #---------------------------------------------------------------------------------------------
        if not os.path.exists( self.trajectoryNameCurr ): os.makedirs( self.trajectoryNameCurr )
        #---------------------------------------------------------------------------------------------
        VelocityVerletDynamics_SystemGeometry(self.molecule                                                      ,
                                              logFrequency              = self.logFreq                           ,
                                              normalDeviateGenerator    = self.RNG                               ,
                                              steps                     = self.nsteps                            ,
                                              timeStep                  = self.timeStep                          ,
                                              trajectories              = [(self.trajectory,self.samplingFactor)],
                                              temperatureScaleFrequency = self.temperatureScaleFreq              ,
                                              temperatureScaleOption    = "linear"            ,
                                              temperatureStart          = self.startTemperature                  ,
                                              temperatureStop           = self.temperature                       )            
    #===============================================================================================    
    def RunProduction(self,_prodSteps,_samplingFactor,_Restricted=False,_equi=False):
        """Execute production or equilibration MD simulation.
        
        Args:
            _prodSteps (int): Number of MD steps.
            _samplingFactor (int): Trajectory save frequency.
            _Restricted (bool): Apply soft constraints. Default: False
            _equi (bool): Whether this is equilibration (vs production). Default: False
        """
        self.softConstraint     = _Restricted
        self.nsteps             = _prodSteps
        self.samplingFactor     = _samplingFactor
        
        if _equi: self.trajectoryNameCurr = os.path.join(self.baseName,self.trajName+"equilibration.ptGeo") 
        else    : self.trajectoryNameCurr = os.path.join(self.baseName,self.trajName+"production.ptGeo") 

        self.trajectory = ExportTrajectory( self.trajectoryNameCurr, self.molecule, log=None )
        if not os.path.exists( self.trajectoryNameCurr ): os.makedirs( self.trajectoryNameCurr )
        if _Restricted: 
            self.trajectorySoft = ExportTrajectory( self.trajectoryNameSoft, self.molecule, log=None )
            if not os.path.exists( self.trajectoryNameSoft ): os.makedirs( self.trajectoryNameSoft )
        if   self.DEBUG: self.Print()        
        if   self.algorithm == "Verlet":      self.runVerlet()
        elif self.algorithm == "LeapFrog":    self.runLeapFrog()
        elif self.algorithm == "Langevin":    self.runLangevin()       
    #===================================================================================================
    def runVerlet(self):
        """Execute Velocity Verlet molecular dynamics.
        
        Runs Velocity Verlet integration with configurable temperature scaling,
        trajectories, and soft constraints for restricted simulations.
        
        Returns:
            None (trajectory saved to file)
        """
        
        trajectory_list = [( self.trajectory, self.samplingFactor)]
        if   self.softConstraint and self.samplingFactor >  0:  trajectory_list = [ ( self.trajectory, self.samplingFactor ), (self.trajectorySoft, 1) ]
        elif self.softConstraint and self.samplingFactor == 0:  trajectory_list = [ (self.trajectorySoft, 1) ]  
        VelocityVerletDynamics_SystemGeometry(  self.molecule                                           ,
                                                logFrequency                = self.logFreq              ,
                                                normalDeviateGenerator      = self.RNG                  ,
                                                steps                       = self.nsteps               ,
                                                timeStep                    = self.timeStep             ,
                                                temperatureScaleFrequency   = self.temperatureScaleFreq ,
                                                temperatureScaleOption      = self.temperatureScaleOption,
                                                trajectories                = trajectory_list           ,
                                                temperatureStart            = self.temperature          )

    #====================================================================================================
    def runLeapFrog(self):
        """Execute Leap Frog molecular dynamics integration.
        
        Performs Leap Frog integration with temperature and pressure control,
        trajectories with optional soft constraints.
        
        Returns:
            None (trajectory saved to file)
        """
        #--------------------------------------------------------------------------------
        trajectory_list = [( self.trajectory, self.samplingFactor)]
        if   self.softConstraint and self.samplingFactor >  0:  trajectory_list = [ ( self.trajectory, self.samplingFactor ), (self.trajectorySoft, 1) ]
        elif self.softConstraint and self.samplingFactor == 0:  trajectory_list = [ (self.trajectorySoft, 1) ]
        LeapFrogDynamics_SystemGeometry(self.molecule                                   ,
                                        trajectories            = trajectory_list       ,
                                        logFrequency            = self.logFreq          ,
                                        normalDeviateGenerator  = self.RNG              ,
                                        pressure                = self.pressure         ,
                                        pressureCoupling        = self.pressureCoupling ,
                                        steps                   = self.nsteps           ,
                                        timeStep                = self.timeStep         ,
                                        temperatureControl      = True                  ,
                                        temperature             = self.temperature      , 
                                        temperatureCoupling     = 0.1                   ) 
    #===================================================================================================
    def runLeapFrog_MT(self):
        """Execute parallel Leap Frog MD using multiple threads.
        
        Distributed version of Leap Frog integration for multi-threaded
        trajectory generation (currently incomplete implementation).
        
        Returns:
            None
        """
        with pymp.Parallel(self.NDIM) as p:
            for i in p.range(self.NDIM):
                trajectory_name = self.trajectoryNameCurr[:-6] + "_" +str(p)+ ".ptGeo" 
                trajectory_soft_name = self.trajectoryNameSoft[:-6] + "_" +str(p)+ ".ptGeo" 

                trajectory = ExportTrajectory( trajectory_name, self.molecule, log=None )
                if not os.path.exists( self.trajectory_name ): os.makedirs( self.trajectory_name )

            if _Restricted: 
                trajectorySoft = ExportTrajectory( trajectory_soft_name, self.molecule, log=None )
                if not os.path.exists( self.trajectory_soft_name ): os.makedirs( self.trajectory_soft_name )


            LeapFrogDynamics_SystemGeometry(self.molecule                                   ,
                                        trajectories            = trajectory_list       ,
                                        logFrequency            = self.logFreq          ,
                                        normalDeviateGenerator  = self.RNG              ,
                                        pressure                = self.pressure         ,
                                        pressureCoupling        = self.pressureCoupling ,
                                        steps                   = self.nsteps           ,
                                        timeStep                = self.timeStep         ,
                                        temperatureControl      = True                  ,
                                        temperature             = self.temperature      , 
                                        temperatureCoupling     = self.temperatureCoupling   ) 
    
    #======================================================================================
    def runLangevin(self):
        """Execute Langevin molecular dynamics integration.
        
        Performs Langevin/stochastic MD with collision frequency control,
        maintaining specified temperature through random forces.
        
        Returns:
            None (trajectory saved to file)
        """
        #-----------------------------------------------------------------------------
        trajectory_list = [( self.trajectory, self.samplingFactor)]
        if   self.softConstraint and self.samplingFactor >  0:  trajectory_list = [ ( self.trajectory, self.samplingFactor ), (self.trajectorySoft, 1) ]
        elif self.softConstraint and self.samplingFactor == 0:  trajectory_list = [ (self.trajectorySoft, 1) ]
        LangevinDynamics_SystemGeometry ( self.molecule                             ,
                                          collisionFrequency     = self.collFreq    ,
                                          logFrequency           = self.logFreq     ,
                                          normalDeviateGenerator = self.RNG         ,
                                          steps                  = self.nsteps      ,
                                          temperature            = self.temperature ,
                                          timeStep               = self.timeStep    ,
                                          trajectories           = trajectory_list  )
    
    #====================================================================================
    def Finalize(self):
        """Finalize trajectory and convert to output format.
        
        Converts trajectory to specified output format (DCD, MDCRD)
        if different from native pDynamo format.
        
        Returns:
            None
        """        
        if self.saveFormat == ".dcd" or self.saveFormat == ".mdcrd":
            if self.saveFormat != self.trajName:
                traj_save = os.path.join(self.baseName,self.trajName + self.saveFormat)
                Duplicate(self.trajectoryNameCurr,traj_save,self.molecule)

        
    #=====================================================================================
    def Print(self):
        """Display MD simulation parameters and configuration.
        
        Prints to console: working directories, simulation time, timestep,
        saving frequency, integrator choice, temperature, pressure settings,
        and constraint information.
        
        Returns:
            None
        """
        ps_time = self.nsteps * self.timeStep
        print( "Molecular Dynamics working folder:{}".format(self.baseName) )
        print( "Molecular Dynamics production trajectory folder:{}".format(self.trajectoryNameCurr) )
        print( "Simulation time (ps): {}".format(ps_time) )
        print( "Simulation steps (n): {}".format(self.nsteps) )
        print( "Frame save frequency: {}".format(self.samplingFactor) )
        print( "Selected Integrator:  {}".format(self.algorithm) )
        print( "temperature (K): {}".format(self.temperature) )
        if self.pressureControl:
            print( "Pressure control applied!")
            print( "Pressure (bar): {}".format(self.pressure) )
        if self.softConstraint:
            print( "Restrictions applied!")
            print( "Molecular Dynamics restraints trajectory folder:{}".format(self.trajectorySoft) )

#================================================================================================#
#======================================END OF THE FILE===========================================#
#================================================================================================#


    
