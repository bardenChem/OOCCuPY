#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Analysis module for molecular dynamics and quantum chemistry calculations.

This module provides comprehensive analysis tools for molecular simulation data,
including trajectory analysis, energy calculations, and potential of mean force (PMF)
computation. It integrates multiple analysis methods to process simulation results
from pDynamo-based calculations.

Classes:
    Analysis: Main analysis class for coordinating various analysis workflows.

Attributes:
    EnergyAnalysis: Energy visualization and 1D/2D plot generation.
    TrajectoryAnalysis: Trajectory-based structural analysis.
    PMF: Potential of mean force calculations from umbrella sampling.
"""

#FILE = Analysis.py
###
#--------------------------------------------------------------
import os, glob, sys
import numpy as np
#--------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#Loading own libraries
#-------------------------------------------------------------
from .EnergyAnalysis     	import EnergyAnalysis
from .TrajectoryAnalysis 	import TrajectoryAnalysis
from .PotentialOfMeanForce   import PMF
from pBabel                    	import *                                     
from pCore                     	import *                                     
from pMolecule                 	import *            
from pScientific               	import *                 
         
from pSimulation               	import *
#-------------------------------------------------------------
#=======================================================================
class Analysis:
	"""Coordinate and execute various molecular analysis workflows.
	
	This class orchestrates multiple analysis types including trajectory analysis,
	energy plotting, PMF calculations, and trajectory splitting. It acts as a central
	hub for accessing specialized analysis modules.
	
	Attributes:
		parameters (dict): Configuration dictionary for analysis parameters.
		molecule: Molecular system object from pDynamo.
		baseFolder (str): Base directory for output files.
	"""
	def __init__(self,_parameters):
		"""Initialize Analysis object with required parameters.
		
		Args:
			_parameters (dict): Configuration dictionary containing:
				- 'active_system': Molecular system (required).
				- 'analysis_type': Type of analysis to perform (required).
				- 'folder': Output folder path (optional, defaults to current directory).
		
		Raises:
			KeyError: If 'analysis_type' parameter is missing.
		"""
		self.parameters = _parameters
		self.molecule   = _parameters["active_system"]

		if "folder" in self.parameters: self.baseFolder = _parameters["folder"]
		else: self.baseFolder = os.getcwd()

		#check parameters
		if not "analysis_type" in self.parameters:
			raise KeyError("Missing required parameter: analysis_type")
		

	#=========================================================================
	def Execute(self):
		"""Execute the specified analysis type.
		
		Dispatches to appropriate analysis method based on the 'analysis_type'
		parameter set during initialization.
		
		Supported analysis types:
			- 'Trajectory_Analysis': Analyze MD trajectory data.
			- 'Energy_Plots': Generate energy landscapes.
			- 'PMF': Calculate potential of mean force.
			- 'Split_Traj': Split trajectory into subsets.
		"""
		Type = self.parameters["analysis_type"]
		if   Type == "Trajectory_Analysis": self.TrajectoryPlots() 		
		elif Type == "Energy_Plots":		self.EnergyPlots()
		elif Type == "PMF":                 self.PMFAnalysis()
		elif Type == "Split_Traj": 			self.SplitTraj()

	#=========================================================================	
	def TrajectoryPlots(self) :
		"""Analyze and plot trajectory statistics (RG, RMSD, distances).
		
		Generates radius of gyration and root mean square deviation plots from
		MD trajectory data. Can optionally calculate and plot reaction coordinate
		distances.
		
		Required parameters:
			nsteps (int): Total simulation time in picoseconds.
		
		Optional parameters:
			show (bool): Display plots interactively.
			calculate_distances (bool): Calculate RC distances.
			ATOMS_RC1 (list): Atom indices for first reaction coordinate.
			ATOMS_RC2 (list): Atom indices for second reaction coordinate.
		"""
		RCs  = None
		if "show" in self.parameters: show = self.parameters["show"]
		t_time = self.parameters["nsteps"]*0.001
		DA = TrajectoryAnalysis(MDrun.trajectoryNameCurr,self.molecule,t_time)
		DA.CalculateRG_RMSD()
		DA.PlotRG_RMS(show)						
		if "calculate_distances" in self.parameters:
			if self.parameters["calculate_distances"] == True:
				rc1 = ReactionCoordinate(self.parameters["ATOMS_RC1"],False,0)
				rc1.GetRCLabel(self.molecule)
				RCs = [rc1]
				rc2 = None
				if "ATOMS_RC2" in self.parameters:
					rc2 = ReactionCoordinate(self.parameters["ATOMS_RC2"],False,0)						
					rc2.GetRCLabel(self.molecule)
					RCs.append(rc2)
				DA.DistancePlots(RCs,show)
	#=========================================================================
	def EnergyPlots(self):
		"""Generate 1D and 2D energy landscape plots from simulation logs.
		
		Creates potential energy surface (PES) visualizations with optional contour
		lines and multiple plot support for comparing different calculation methods.
		
		Required parameters:
			xsize (int): Number of points in first coordinate.
			type (str): Plot type ('1D', '2D', etc.).
			log_name (str): Path to energy log file.
		
		Optional parameters:
			ysize (int): Number of points in second coordinate for 2D plots.
			contour_lines (int): Number of contour levels.
			xlim, ylim (list): Axis limits.
			show (bool): Display plots.
			multiple_plot (bool): Plot multiple methods.
			retrieve_path (str): Extract lowest-energy path.
		"""
		#-------------------------------------------------------------------
		#detect xsize and ysize if not provided
		xindx= []
		yindx= []
		if self.parameters["xsize"] < 0:
			print("Auto-detecting xsize from log file...")
			with open(self.parameters["log_name"], 'r') as log_file:				
				lines = log_file.readlines()
				if lines[0].split()[0] == "x":
					for line in lines:
						line2 = line.split()
						try: xindx.append(int(line2[0]))
						except: pass
				if lines[0].split()[1] == "y":
					for line in lines:
						line2 = line.split()
						try: yindx.append(int(line2[1]))
						except: pass
		self.parameters["xsize"] = max(xindx) + 1
		if len(yindx) > 0: self.parameters["ysize"] = max(yindx) + 1
		#-------------------------------------------------------------------
		#detect xlim and ylim if not provided
		if not "xlim" in self.parameters: self.parameters["xlim"] = [0,0]
		if not "ylim" in self.parameters: self.parameters["ylim"] = [0,0]
		if self.parameters["xlim"] == [0,0]:
			rc1 = [] 
			rc2 = []
			with open(self.parameters["log_name"], 'r') as log_file:				
				lines = log_file.readlines()
				if lines[0].split()[2] == "RC1":
					for line in lines:
						line2 = line.split()
						try: rc1.append(float(line2[2]))
						except: pass
				if lines[0].split()[3] == "RC2":
					for line in lines:
						line2 = line.split()
						try: rc2.append(float(line2[3]))
						except: pass

			self.parameters["xlim"] = [ rc1[0] , rc1[-1] ]
			if self.parameters["ysize"] > 0:
				self.parameters["ylim"] = [rc2[0] , rc2[-1] ]
			
		multiPlot = False
		ndim      = 1
		try: crd1_label= self.molecule.reactionCoordinates[0].label
		except: 
			crd1_label = "Reaction Path Frames (n)"
			pass
		try: crd2_label= self.molecule.reactionCoordinates[1].label
		except:
			crd2_label = "No label"
			pass

		cnt_lines = 0 
		ysize     = 0
		if "ysize" in self.parameters: ysize = self.parameters["ysize"]
		xlim      = [ 0, self.parameters["xsize"] ]
		ylim 	  = [ 0, ysize ]
		show 	  = False

		in_point  = [0,0]
		fin_point = [0,0]
		FindPath  = False
		#--------------------------------------------------------
		if "contour_lines" 	in self.parameters: cnt_lines  = self.parameters["contour_lines"]		
		if "xlim" 		in self.parameters: xlim  	   = self.parameters["xlim"		    ]
		if "ylim" 		in self.parameters: ylim       = self.parameters["ylim"   		]
		if "show" 			in self.parameters: show       = self.parameters["show"         ]
		if "in_point"       in self.parameters: in_point   = self.parameters["in_point"     ]
		if "fin_point"      in self.parameters: fin_point  = self.parameters["fin_point"    ]
		if "multiple_plot" 	in self.parameters: multiPlot  = True 		
		if ysize > 0: ndim = 2
		#--------------------------------------------------------
		EA = EnergyAnalysis(self.parameters["xsize"],ysize,_type=self.parameters["type"] )
		
		EA.ReadLog(self.parameters["log_name"] )
		if multiPlot:
			print("Multiplot_required")
			if   ndim == 1: EA.MultPlot1D(label=crd1_label)
			elif ndim == 2: EA.MultPlot2D(cnt_lines,crd1label=crd1_label,crd2label=crd2_label,_xlim=xlim,_ylim=ylim,SHOW=show) 
		#--------------------------------------------------------
		elif ndim == 1: EA.Plot1D(crd1_label,XLIM=xlim,SHOW=show)
		elif ndim == 2:	EA.Plot2D(cnt_lines,crd1_label,crd2_label,xlim,ylim,show)

		if  "retrieve_path" in self.parameters: 
			if fin_point[0] < 0: fin_point[0] = self.parameters["xsize"] - 1
			if fin_point[1] < 0: fin_point[1] = self.parameters["ysize"] - 1
			EA.Path_From_PES(in_point,fin_point,self.parameters["retrieve_path"],self.baseFolder,self.molecule.system)

	#=========================================================================
	def PMFAnalysis(self):
		"""Calculate and plot potential of mean force (PMF) from umbrella sampling data.
		
		Computes PMF and free energy surfaces using the WHAM algorithm from
		restricted MD or umbrella sampling simulations.
		
		Required parameters:
			source_folder (str): Folder containing umbrella sampling windows.
			xnbins (int): Number of bins in first dimension.
			temperature (float): Simulation temperature in Kelvin.
		
Optional parameters:
			ynbins (int): Number of bins in second dimension (default=0 for 1D).
			contour_lines (int): Number of contour levels for 2D plots.
			xlim_list, ylim_list (list): Axis limits.
			crd1_label, crd2_label (str): Axis labels.
			show (bool): Display plots.
			oneDimPlot (bool): Show 1D projection of 2D PMF.
		"""
		ynbins = 0 
		if "ynbins" in self.parameters: ynbins = self.parameters["ynbins"]
		potmean = PMF( self.molecule, self.parameters["source_folder"], self.baseFolder )
		potmean.CalculateWHAM(self.parameters["xnbins"],ynbins,self.parameters["temperature"])
		#================================================================
		#Set default plot parameters
		cnt_lines  = 12
		crd1_label = ""
		crd2_label = ""
		nRC2       = ynbins
		show       = False
		xwin       = 0
		ywin       = 0 
		#-----------------------------------------------------------------
		nDims = 1
		if ynbins > 0: nDims = 2
		xlims = [ 0,  self.parameters['xnbins'] ]
		ylims = [ 0,  ynbins ]
		OneDimPlot = False
		#-------------------------------------------------------------
		#check parameters for plot
		if "contour_lines" 	in self.parameters: cnt_lines  = self.parameters["contour_lines"]		
		if "xlim_list" 		in self.parameters: xlims 	   = self.parameters["xlim_list"]
		if "ylim_list" 		in self.parameters:	ylims 	   = self.parameters["ylim_list"]
		if "show" 			in self.parameters:	show 	   = self.parameters["show"]
		if "crd1_label" 	in self.parameters:	crd1_label = self.parameters["crd1_label"]
		if "crd2_label" 	in self.parameters:	crd2_label = self.parameters["crd2_label"]
		if "xwindows" 		in self.parameters: xwin 	   = self.parameters["xwindows"]
		if "ywindows" 		in self.parameters:	ywin 	   = self.parameters["ywindows"]
		if "oneDimPlot"     in self.parameters: OneDimPlot = self.parameters["oneDimPlot"]
		#------------------------------------------------------------
		if   nDims == 2: TYPE = "WHAM2D"
		elif nDims == 1: TYPE = "WHAM1D"	
		#------------------------------------------------------------
		# Plot PMF graphs
		EA = EnergyAnalysis(self.parameters['xnbins'],nRC2,_type=TYPE)
		EA.ReadLog( os.path.join(potmean.baseName,"PotentialOfMeanForce.dat") ) 
		#-------------------------------------------------------------
		xlims = [ np.min(EA.RC1), np.max(EA.RC1) ]
		if   nDims == 2: 
			ylims = [ np.min(EA.RC2), np.max(EA.RC2) ]
			EA.Plot2D(cnt_lines,crd1_label,crd2_label,xlims,ylims,show)
		elif nDims == 1: EA.Plot1D(crd1_label,SHOW=show)
		#-------------------------------------------
		#Plot Free energy of the calculated windows
		if OneDimPlot: TYPE = "FE1D"
		elif nDims 	  == 2: 	  TYPE = "FE2D"
		elif nDims 	  == 1: 	  TYPE = "FE1D"
		

		xlims = [ np.min(EA.RC1), np.max(EA.RC1) ]

		if nDims  == 2:  ylims = [ np.min(EA.RC2), np.max(EA.RC2) ]	
		#------------------------------------------
		EAfe = EnergyAnalysis(self.parameters["xnbins"],ynbins,_type=TYPE)
		EAfe.ReadLog( os.path.join(potmean.baseName,"FreeEnergy.log") ) 
		#-------------------------------------------------------------
		if nDims == 2: 
			if OneDimPlot: EAfe.Plot1D_FreeEnergy(crd1_label,crd2_label,show)
			else 		 : EAfe.Plot2D(cnt_lines,crd1_label,crd2_label,xlims,ylims,show)
		elif nDims == 1: EAfe.Plot1D(crd1_label,XLIM=xlims,SHOW=show)
	#====================================================================	
	def SplitTraj(self):
		'''
		'''
		trj = TrajectoryAnalysis(self.parameters["trajectory_name"],self.molecule,0)
		trj.Split_Traj(self.parameters["break_point"])

#==================================================================================