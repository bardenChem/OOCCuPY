#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = GeometrySearcher.py

#==============================================================================
import os, sys
#importing our library functions
from . import commonFunctions
from .logging_config import get_logger, log_step_start, log_step_end, log_checkpoint
# pDynamo
from pBabel                    import *                                     
from pCore                     import *                                     
from pMolecule                 import *                  
from pScientific               import *                                     
from pScientific.Arrays        import *                                     
from pScientific.Geometry3     import *                 
from pSimulation               import *

from .TrajectoryAnalysis  import TrajectoryAnalysis
#*********************************************************************************
class GeometrySearcher:
    '''
    Class to handle with pDynamo methods that search geometries for the system, such as global/local minimuns
    as saddle points and reaction path trajectories. 
    '''
    #.-------------------------------------------------------------------------   
    def __init__(self,_system,_baseFolder,_trajName=None, _enable_debug_file=False, _verbosity="INFO"):
        '''
        Class constructor.
        Parameters:
        -----------
        _system : pDynamo system object
        _baseFolder : str, base directory path
        _trajName : str, optional trajectory name
        _enable_debug_file : bool, if True creates separate DEBUG log file
        _verbosity : str, logging verbosity level ("DEBUG", "INFO", "WARNING", "ERROR")
        '''
        self.molecule       = _system
        self.baseName       = _baseFolder
        self.optAlg         = "ConjugatedGradient"
        self.InitCrd3D      = Clone(_system.coordinates3)
        self.finalCrd3D     = None
        self.massWeighting  = False
        self.logFreq        = 50 # deafult value for otimizations, must to be changed through the specific class method
        self.trajectoryName = None
        self.savePdb        = False
        self.saveFormat     = None    
        self.rmsGrad        = 0.1
        self.maxIt          = 500
        self.saveFrequency  = 0
        self.DEBUG          = False
        if not _trajName == None: self.trajectoryName = os.path.join(_baseFolder,_trajName)
        
        # Initialize logger
        self.logger = get_logger(
            "GeometrySearcher",
            output_dir=_baseFolder,
            enable_debug_file=_enable_debug_file,
            verbosity=_verbosity
        )
        self.logger.debug(f"GeometrySearcher initialized for {_baseFolder}")
    #=========================================================================
    def ChangeDefaultParameters(self,_parameters):
        '''
        Class method to modify default parameters for the minimization runs
        '''       
        if "save_pdb"       in _parameters: self.savePdb        = _parameters["save_pdb"]
        if "maxIterations"  in _parameters: self.maxIt          = _parameters['maxIterations']            
        if "log_frequency"  in _parameters: self.logFreq        = _parameters["log_frequency"] 
        if "save_format"    in _parameters: self.saveFormat     = _parameters["save_format"]
        if "save_frequency" in _parameters: self.saveFrequency  = _parameters["save_frequency"]        
        if "rmsGradient"    in _parameters: self.rmsGrad        = _parameters["rmsGradient"]
        if "Debug"          in _parameters: self.DEBUG          = _parameters["Debug"]
        
        # Log parameter changes
        self.logger.debug(f"Parameters updated: savePdb={self.savePdb}, maxIt={self.maxIt}, "
                         f"logFreq={self.logFreq}, rmsGrad={self.rmsGrad}, DEBUG={self.DEBUG}")
    #======================================================================================
    # Main minimization class method
    def Minimization(self,_optmizer):
        '''
        Execute the minimization routine for search of geometry corresponding to local minima
        '''
        #------------------------------------------------------------------
        log_step_start(self.logger, "Minimization", algorithm=_optmizer, max_iterations=self.maxIt)
        self.optAlg = _optmizer                 
        # run the minimization for the chosen algorithm
        try:
            if   self.optAlg == "ConjugatedGradient": self.RunConjugatedGrad()
            elif self.optAlg == "SteepestDescent"   : self.RunSteepestDescent()
            elif self.optAlg == "LFBGS"             : self.RunLFBGS()
            elif self.optAlg == "QuasiNewton"       : self.RunQuasiNewton()
            elif self.optAlg == "FIRE"              : self.RunFIREmin()
            else:
                self.logger.warning(f"Unknown algorithm: {self.optAlg}")
            
            log_checkpoint(self.logger, "Minimization completed", status="OK")
        except Exception as e:
            self.logger.error(f"Minimization failed with error: {str(e)}")
            raise
        
        self.finalCrd3D = Clone(self.molecule.coordinates3)
        if self.DEBUG:
            self.Print()
            pdbFileA = os.path.join(self.baseName, "initialCoord_{}.pdb".format(self.optAlg) )
            pdbFileB = os.path.join(self.baseName, "finalCoord_{}.pdb".format(self.optAlg) )
            self.molecule.coordinates3 = Clone(self.InitCrd3D)
            ExportSystem(pdbFileA,self.molecule)
            self.molecule.coordinates3 = Clone(self.finalCrd3D)
            ExportSystem(pdbFileB,self.molecule)
        
        log_step_end(self.logger, "Minimization")
    #=============================================================================
    #Minimizers methods
    def RunConjugatedGrad(self):
        '''
        Class method to apply the conjugated gradient minimizer
        '''
        self.logger.debug("Starting Conjugated Gradient minimization...")
        try:
            if self.trajectoryName == None:
                self.logger.debug("Running without trajectory output")
                ConjugateGradientMinimize_SystemGeometry(self.molecule                      ,                
                                                         logFrequency           = self.logFreq  ,
                                                         maximumIterations      = self.maxIt    ,
                                                         rmsGradientTolerance   = self.rmsGrad  )
            else:
                self.logger.debug(f"Running with trajectory output to {self.trajectoryName}")
                trajectory = ExportTrajectory( self.trajectoryName, self.molecule, log=None )
                ConjugateGradientMinimize_SystemGeometry(self.molecule                        ,                
                                                     logFrequency           = self.logFreq    ,
                                                     trajectories           = [(trajectory, self.saveFrequency)],
                                                     maximumIterations      = self.maxIt      ,
                                                     rmsGradientTolerance   = self.rmsGrad    )
            self.logger.debug("Conjugated Gradient minimization completed")
        except Exception as e:
            self.logger.error(f"Conjugated Gradient minimization failed: {str(e)}")
            raise

    #=====================================================================================
    def RunSteepestDescent(self):
        '''
        Class method to apply the steepest descent minimizer
        '''
        self.logger.debug("Starting Steepest Descent minimization...")
        try:
            if self.trajectoryName == None:
                self.logger.debug("Running without trajectory output")
                SteepestDescentMinimize_SystemGeometry(self.molecule                       ,               
                                                    logFrequency            = self.logFreq ,
                                                    maximumIterations       = self.maxIt   ,
                                                    rmsGradientTolerance    = self.rmsGrad )
            else:
                self.logger.debug(f"Running with trajectory output to {self.trajectoryName}")
                trajectory = ExportTrajectory( self.trajectoryName, self.molecule, log=None  )
                SteepestDescentMinimize_SystemGeometry(self.molecule                       ,               
                                                    logFrequency            = self.logFreq ,
                                                    trajectories            = [(trajectory, self.saveFrequency)],
                                                    maximumIterations       = self.maxIt   ,
                                                    rmsGradientTolerance    = self.rmsGrad )
            self.logger.debug("Steepest Descent minimization completed")
        except Exception as e:
            self.logger.error(f"Steepest Descent minimization failed: {str(e)}")
            raise
    #============================================================================
    def RunLFBGS(self):
        '''
        Class method to apply the LFBGS minimizer
        '''
        self.logger.debug("Starting LFBGS minimization...")
        try:
            if self.trajectoryName == None:
                self.logger.debug("Running without trajectory output")
                LBFGSMinimize_SystemGeometry(self.molecule                          ,                
                                        logFrequency         = self.logFreq         ,
                                        maximumIterations    = self.maxIt           ,
                                        rmsGradientTolerance = self.rmsGrad         )
            else:
                self.logger.debug(f"Running with trajectory output to {self.trajectoryName}")
                trajectory = ExportTrajectory( self.trajectoryName, self.molecule, log=None )
                LBFGSMinimize_SystemGeometry(self.molecule                          ,                
                                        logFrequency         = self.logFreq         ,
                                        trajectories         = [(trajectory, self.saveFrequency)],
                                        maximumIterations    = self.maxIt           ,
                                        rmsGradientTolerance = self.rmsGrad         )
            self.logger.debug("LFBGS minimization completed")
        except Exception as e:
            self.logger.error(f"LFBGS minimization failed: {str(e)}")
            raise    
    #=============================================================================
    def RunQuasiNewton(self):
        '''
        Class method to apply the Quaisi-Newton minimizer
        '''
        self.logger.debug("Starting Quasi-Newton minimization...")
        try:
            if self.trajectoryName == None: 
                self.logger.debug("Running without trajectory output")
                QuasiNewtonMinimize_SystemGeometry( self.molecule                       ,                
                                                    logFrequency         = self.logFreq ,
                                                    maximumIterations    = self.maxIt   ,
                                                    rmsGradientTolerance = self.rmsGrad )
            else:
                self.logger.debug(f"Running with trajectory output to {self.trajectoryName}")
                trajectory = ExportTrajectory( self.trajectoryName, self.molecule, log=None )
                QuasiNewtonMinimize_SystemGeometry( self.molecule                       ,                
                                                    logFrequency         = self.logFreq ,
                                                    trajectories         = [(trajectory, self.saveFrequency)],
                                                    maximumIterations    = self.maxIt   ,
                                                    rmsGradientTolerance = self.rmsGrad )
            self.logger.debug("Quasi-Newton minimization completed")
        except Exception as e:
            self.logger.error(f"Quasi-Newton minimization failed: {str(e)}")
            raise
    #==============================================================================
    def RunFIREmin(self):
        '''
        Class method to apply the FIRE minimizer
        '''
        self.logger.debug("Starting FIRE minimization...")
        try:
            if self.trajectoryName == None:
                self.logger.debug("Running without trajectory output")
                FIREMinimize_SystemGeometry( self.molecule                  ,                
                                             logFrequency         = self.logFreq ,
                                             maximumIterations    = self.maxIt   ,
                                             rmsGradientTolerance = self.rmsGrad )
            else:
                self.logger.debug(f"Running with trajectory output to {self.trajectoryName}")
                trajectory = ExportTrajectory( self.trajectoryName, self.molecule, log=None )
                FIREMinimize_SystemGeometry( self.molecule                                            ,                
                                             logFrequency         = self.logFreq                      ,
                                             trajectories         = [(trajectory, self.saveFrequency)],
                                             maximumIterations    = self.maxIt                        ,
                                             rmsGradientTolerance = self.rmsGrad                      )
            self.logger.debug("FIRE minimization completed")
        except Exception as e:
            self.logger.error(f"FIRE minimization failed: {str(e)}")
            raise        
    #=============================================================================
    # Reaction path searchers
    def NudgedElasticBand(self,_parameters):
        '''
        Nudget Elastic Band procedure to estimate a reaction path
        '''
        #-------------------------------------------------------------------------
        has_traj_source = "traj_source" in _parameters
        log_step_start(self.logger, "NudgedElasticBand (NEB)", 
                      trajectory_source="PROVIDED" if has_traj_source else "CREATE_FROM_COORDS")
        
        try:
            rmdGIS          = 1
            springCF        = 500.0
            fixedTerminal   = False
            useSpline       = False
            spline_tol      = 1.5
            traj_bins       = 0
            overwrite       = False 
            
            # Parse parameters
            self.logger.debug("Parsing NEB parameters...")
            if "spring_constant_force" in _parameters: 
                springCF = _parameters["spring_constant_force"]
            if "fixed_terminal_images" in _parameters: 
                fixedTerminal = _parameters["fixed_terminal_images"]
            if "RMS_growing_intial_string" in _parameters: 
                rmsGIS = _parameters["RMS_growing_intial_string"]
            if "spline_redistribution" in _parameters: 
                useSpline = _parameters["spline_redistribution"]
            if "overwrite_traj" in _parameters: 
                if _parameters["overwrite_traj"] == "yes":
                    overwrite = True
            if "traj_bins" in _parameters: 
                traj_bins = _parameters["traj_bins"]
            
            msg = f"NEB parameters: spring_constant={springCF}, fixed_terminal={fixedTerminal}, "
            msg += f"spline_tolerance={spline_tol}, max_iterations={self.maxIt}, rms_gradient={self.rmsGrad}"
            self.logger.debug(msg)

            self.trajectoryName = os.path.join(self.baseName,"NEB.ptGeo")        
            trajectory = None
            log_checkpoint(self.logger, "NEB trajectory path setup", status="OK", 
                          message=f"Path: {self.trajectoryName}")
            
            #-----------------------------------------------------------------------------------------
            # Two modes: create from coordinates OR load from trajectory source
            if not has_traj_source:
                log_checkpoint(self.logger, "NEB Mode: Create from initial and final coordinates", status="OK")
                
                if "init_coord" in _parameters: 
                    self.logger.info(f"Loading initial coordinates from {_parameters['init_coord']}")
                    self.InitCrd3D = ImportCoordinates3(_parameters["init_coord"], log=None)
                    log_checkpoint(self.logger, "Initial coordinates loaded", status="OK")
                    
                if "final_coord" in _parameters: 
                    self.logger.info(f"Loading final coordinates from {_parameters['final_coord']}")
                    self.finalCrd3D = ImportCoordinates3(_parameters["final_coord"], log=None)
                    log_checkpoint(self.logger, "Final coordinates loaded", status="OK")
                
                self.logger.info(f"Creating growing string path with {traj_bins} intermediate frames...")
                GrowingStringInitialPath(self.molecule,
                                        _parameters["traj_bins"],
                                        self.InitCrd3D,
                                        self.finalCrd3D,
                                        self.trajectoryName,
                                        rmsGradientTolerance=rmsGIS)
                log_checkpoint(self.logger, "Growing string path created", status="OK")
                
                self.logger.debug("Exporting trajectory for optimization...")
                trajectory = ExportTrajectory(self.trajectoryName, self.molecule, append=True)
                log_checkpoint(self.logger, "Trajectory exported", status="OK")
                
            else:
                # TRAJECTORY SOURCE MODE - This is where the segmentation fault likely occurs
                log_checkpoint(self.logger, "NEB Mode: Load from trajectory source", status="OK")
                self.logger.info(f"Trajectory source: {_parameters['traj_source']}")
                
                import shutil
                from pathlib import Path

                src = Path(_parameters["traj_source"])
                self.logger.debug(f"Source path exists: {src.exists()}")
                log_checkpoint(self.logger, "Source path validation", 
                              status="OK" if src.exists() else "ERROR", 
                              message=f"Source exists: {src.exists()}")
                
                self.trajectoryName = os.path.join(self.baseName,"NEB.ptGeo")
                self.logger.debug(f"Target trajectory path: {self.trajectoryName}")
                
                if not os.path.exists(self.trajectoryName): 
                    self.logger.info(f"Creating trajectory directory: {self.trajectoryName}")
                    os.makedirs(self.trajectoryName)
                    log_checkpoint(self.logger, "Trajectory directory created", status="OK")
                else:
                    self.logger.debug("Trajectory directory already exists")
                    if overwrite:
                        self.logger.info("Overwrite flag set - copying source to target...")
                        shutil.copytree(src, self.trajectoryName, 
                                      ignore=shutil.ignore_patterns('*.png'), dirs_exist_ok=True)
                        log_checkpoint(self.logger, "Trajectory directory overwritten", status="OK")
                    else:
                        self.logger.debug("Overwrite flag not set - reusing existing trajectory directory")
                
                self.logger.debug("Exporting trajectory from source...")
                trajectory = ExportTrajectory(self.trajectoryName, self.molecule, append=True)
                traj_bins = trajectory.numberOfFrames
                self.logger.info(f"Trajectory loaded: {traj_bins} frames")
                log_checkpoint(self.logger, "Trajectory from source loaded", status="OK", 
                              message=f"Frames: {traj_bins}")
            
            #------------------------------------------------------------------------------------------
            # Chain of States Optimization
            self.logger.info("Starting Chain of States (NEB) optimization...")
            msg = f"Optimization parameters: fixedTerminal={fixedTerminal}, "
            msg += f"springForceConstant={springCF}, useSpline={useSpline}"
            self.logger.debug(msg)
            
            ChainOfStatesOptimizePath_SystemGeometry(self.molecule,
                                                    trajectory,
                                                    logFrequency=1,
                                                    maximumIterations=self.maxIt,
                                                    fixedTerminalImages=fixedTerminal,
                                                    springForceConstant=springCF,
                                                    splineRedistributionTolerance=spline_tol,
                                                    forceSplineRedistributionCheckPerIteration=useSpline,
                                                    rmsGradientTolerance=self.rmsGrad)
            
            log_checkpoint(self.logger, "Chain of States optimization completed", status="OK")
            self.logger.info(f"NEB optimization complete with {traj_bins} frames")
            log_step_end(self.logger, "NudgedElasticBand (NEB)", frames=traj_bins)
            
            return traj_bins
            
        except Exception as e:
            self.logger.error(f"NEB failed with exception: {type(e).__name__}: {str(e)}")
            import traceback
            self.logger.debug(f"Traceback: {traceback.format_exc()}")
            raise
        
    #========================================================================================
    def SelfAvoidWalking(self,_parameters):
        '''
        Self-Avoid-Walking procedure to estimate a reaction path
        '''       
        self.trajectoryName = self.baseName + "SAW.ptGeo"
        self.traj = ExportTrajectory( self.trajectoryName, self.molecule, append=True ) 
        ExpandByLinearInterpolation( _parameters["traj_source"], self.trajectoryName, self.molecule, _parameters["traj_bins"])
        Gamma = 100.0
        Rho   = 2.0
        Kappa = 5000.0
        if "gamma" in _parameters: Gamma = _parameters["gamma"]
        if "rho"   in _parameters: Rho   = _parameters["rho"]
        if "kappa" in _parameters: Kappa = _parameters["kappa"]
        SAWOptimize_SystemGeometry ( self.molecule, self.traj, gamma=Gamma, kappa=Kappa )
    #========================================================================================
    def SteepestDescentPathSearch(self,_parameters):
        '''
        '''
        massW    = True
        funcStep = 2.0
        pathStep = 0.025 

        if "mass_weighting" in _parameters: massw    = _parameters["mass_weighting"]
        if "function_step"  in _parameters: funcStep = _parameters["function_step"]
        if "path_step"      in _parameters: pathStep = _parameters["path_step"]

        self.molecule.coordinates3 = _parameters["saddle_conformation"]
        self.trajectoryName = self.baseName + ".steepPath.ptGeo"
        self.traj = ExportTrajectory( self.trajectoryName, self.molecule )
        SteepestDescentPath_SystemGeometry( self.molecule                           ,
                                            functionStep      = funcStep            ,
                                            logFrequency      = self.logFrequency   ,
                                            maximumIterations = self.maxIt          ,
                                            pathStep          = pathStep            ,
                                            saveFrequency     = self.save_frequency ,
                                            trajectory        = self.traj           ,
                                            useMassWeighting  = massW               )

    #========================================================================================
    def BakerSaddleOptimizer(self,_parameters):
        '''
        Class method to search saddle-points transition structure
        '''

        self.InitCrd3D = ImportCoordinates3(_parameters["saddle_coord"] )
        self.molecule.coordinates3 = Clone(self.InitCrd3D)
        BakerSaddleOptimize_SystemGeometry( self.molecule                       ,
                                            logFrequency         =      1       ,
                                            maximumIterations    = self.maxIt   ,
                                            rmsGradientTolerance = self.rmsGrad )

        self.finalCrd3D = Clone(self.molecule.coordinates3)
        Pickle(self.baseName+"_BakerOpt.pkl",self.finalCrd3D)
        if savePdb: 
            ExportSystem(self.baseName+"_BakerOpt.pdb",self.finalCrd3D)
            savePdb = False

    #=========================================================================================
    def CalculateRMS(self):
        '''
        Calculate the root mean square of deviation of the final coordinate found with the first set given.
        '''
        try:
            self.logger.debug("Calculating RMS with mass weighting...")
            masses = Array.FromIterable ( [ atom.mass for atom in self.molecule.atoms ] )
            self.InitCrd3D.Superimpose ( self.finalCrd3D, weights = masses )
            rms = self.InitCrd3D.RootMeanSquareDeviation ( self.finalCrd3D, weights = masses )
            msg = f"Root Mean Square of Deviation of the optimized structure from the initial: {rms:.6f}"
            print(msg)
            self.logger.info(msg)
        except Exception as e:
            self.logger.error(f"RMS calculation failed: {str(e)}")
            raise
    #===========================================================================================
    def Finalize(self):
        '''
        Finaluze the Geometry searcher procedures, save structures and/or trajectories.
        This is a critical section where memory issues and segmentation faults can occur.
        '''
        log_step_start(self.logger, "Finalize", savePdb=self.savePdb, saveFrequency=self.saveFrequency)
        
        try:
            # Calculate and log RMS deviation
            self.logger.debug("Calculating RMS deviation...")
            self.CalculateRMS()
            log_checkpoint(self.logger, "RMS calculation completed", status="OK")
            
            #----------------------------------------------------------------------
            # Save structures and/or trajectories
            if self.savePdb:
                self.logger.info("Saving PDB structure...")
                try:
                    pdbFile = self.baseName + "opt_{}.pdb".format(self.optAlg)
                    i = 0
                    while os.path.exists(pdbFile):
                        pdbFile = self.baseName + "_#{}_opt_{}.pdb".format(i,self.optAlg)
                        i += 1
                    self.logger.debug(f"Writing PDB to: {pdbFile}")
                    ExportSystem(pdbFile,self.molecule)
                    log_checkpoint(self.logger, "PDB structure saved", status="OK", message=pdbFile)
                except Exception as e:
                    self.logger.error(f"Failed to save PDB: {str(e)}")
                    log_checkpoint(self.logger, "PDB save failed", status="ERROR", message=str(e))
                    raise
            
            #----------------------------------------------------------------------
            # Save trajectory if frequency > 0
            if self.saveFrequency > 0:
                self.logger.debug(f"Save frequency: {self.saveFrequency}")
                if self.trajectoryName:
                    self.logger.info(f"Processing trajectory: {self.trajectoryName}")
                    log_checkpoint(self.logger, "Trajectory processing started", status="OK")
                    
                    if self.saveFormat == ".dcd" or self.saveFormat == ".mdcrd":
                        if self.saveFormat != self.trajectoryName:
                            try:
                                self.logger.debug(f"Converting trajectory to {self.saveFormat} format...")
                                traj_save = os.path.splitext(self.trajectoryName)[0] + self.saveFormat
                                self.logger.debug(f"Target trajectory file: {traj_save}")
                                log_checkpoint(self.logger, "Trajectory conversion starting", status="OK")
                                
                                Duplicate(self.trajectoryName,traj_save,self.molecule)
                                log_checkpoint(self.logger, "Trajectory converted", status="OK", 
                                              message=f"Format: {self.saveFormat}")
                                
                            except Exception as e:
                                self.logger.warning(f"Trajectory conversion failed: {str(e)}")
                                log_checkpoint(self.logger, "Trajectory conversion failed", status="WARNING", 
                                              message=str(e))
                    
                    # Analyze trajectory
                    try:
                        self.logger.debug("Starting trajectory analysis...")
                        import glob
                        pkl_files = glob.glob( os.path.join( self.trajectoryName,"*.pkl") )
                        xsi = len(pkl_files)
                        self.logger.debug(f"Found {xsi} frames in trajectory")
                        log_checkpoint(self.logger, "Trajectory frame count", status="OK", message=f"Frames: {xsi}")
                        
                        if xsi > 0:
                            self.logger.info(f"Creating trajectory analysis object for {xsi} frames...")
                            trajAn = TrajectoryAnalysis(self.trajectoryName,self.molecule,xsi)
                            log_checkpoint(self.logger, "Trajectory analysis object created", status="OK")
                            
                            self.logger.info("Calculating Rg and RMSD...")
                            trajAn.CalculateRG_RMSD(qc_mm=False)
                            log_checkpoint(self.logger, "Rg and RMSD calculated", status="OK")
                            
                            self.logger.info("Plotting Rg and RMS...")
                            trajAn.PlotRG_RMS()
                            log_checkpoint(self.logger, "Trajectory plots created", status="OK")
                        else:
                            self.logger.warning("No frames found in trajectory!")
                            log_checkpoint(self.logger, "No trajectory frames found", status="WARNING")
                            
                    except Exception as e:
                        self.logger.error(f"Trajectory analysis failed: {str(e)}")
                        import traceback
                        self.logger.debug(f"Traceback: {traceback.format_exc()}")
                        log_checkpoint(self.logger, "Trajectory analysis failed", status="ERROR", message=str(e))
                else:
                    self.logger.debug("No trajectory name defined")
            
            log_step_end(self.logger, "Finalize")
            self.logger.info("="*70)
            self.logger.info("FINALIZE COMPLETED SUCCESSFULLY")
            self.logger.info("="*70)
            
        except Exception as e:
            self.logger.error(f"CRITICAL: Finalize failed with exception: {type(e).__name__}: {str(e)}")
            import traceback
            self.logger.error(f"Full traceback:\n{traceback.format_exc()}")
            log_checkpoint(self.logger, "Finalize failed", status="ERROR", message=f"{type(e).__name__}: {str(e)}")
            raise

    #===========================================================================================
    def Print(self):
        '''
        Print to screen basic info for the simulation. 
        '''        
        info_msg = f"\nGeometry Searcher Configuration:\n"
        info_msg += f"  Trajectory folder: {self.trajectoryName}\n"
        info_msg += f"  RMS gradient tolerance: {self.rmsGrad}\n"
        info_msg += f"  Optimization Algorithm: {self.optAlg}\n"
        info_msg += f"  Maximum iterations: {self.maxIt}\n"
        info_msg += f"  Log frequency: {self.logFreq}\n"
        print(info_msg)
        self.logger.info(info_msg)
        
#================================================================================================#
#======================================END OF THE FILE===========================================#
#================================================================================================#
