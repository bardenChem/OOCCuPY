#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
	'''
	'''
	def __init__(self,_parameters):
		'''
		Default Constructor
		'''
		self.parameters = _parameters
		self.molecule   = _parameters["active_system"]

		if "folder" in self.parameters: self.baseFolder = _parameters["folder"]
		else: self.baseFolder = os.getcwd()

		#check parameters
		if not "analysis_type" in self.parameters:
			raise KeyError("Missing required parameter: analysis_type")
		

	#=========================================================================
	def Execute(self):
		'''
		'''
		Type = self.parameters["analysis_type"]
		if   Type == "Trajectory_Analysis": self.TrajectoryPlots() 		
		elif Type == "Energy_Plots":		self.EnergyPlots()
		elif Type == "PMF":                 self.PMFAnalysis()
		elif Type == "Split_Traj": 			self.SplitTraj()

	#=========================================================================	
	def TrajectoryPlots(self) :
		'''
		Mandatory keys in self.parameters:
		Optional keys in self.parameters:
		'''
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
		'''
		Produce Energy plots from previus simulations log files
		Mandatory keys in self.parameters:
		Optional keys in self.parameters:
		'''		
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
		if "xlim_list" 		in self.parameters: xlim  	   = self.parameters["xlim_list"    ]
		if "ylim_list" 		in self.parameters: ylim       = self.parameters["ylim_list"    ]
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
			EA.Path_From_PES(in_point,fin_point,self.parameters["retrieve_path"],self.baseFolder,self.molecule.system)

	#=========================================================================
	def PMFAnalysis(self):
		'''
		Calculate potential of mean force and Free energy from restricted molecular dynamics
		Mandatory keys: 
			"source_folder"	:
			"xbins"			:
			"ybins"			:
			"temperature"	:
		Optinal keys        :
		plot keys           :				
		'''
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
		EAfe = EnergyAnalysis(xwin,ywin,_type=TYPE)
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