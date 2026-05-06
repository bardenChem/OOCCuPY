#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Energy Analysis and Visualization Module.

This module contains tools for analyzing and plotting energy data from molecular
simulations. It supports 1D and 2D energy surfaces, multiple calculation methods,
and various energy plot types including potentials, free energies, and PMF surfaces.

Classes:
    EnergyAnalysis: Central class for energy data management and visualization.

Supported plot types:
    - 1D/2D: Standard potential energy surfaces.
    - WHAM1D/WHAM2D: PMF from WHAM analysis.
    - FE1D/FE2D: Free energy surfaces.
"""

#FILE = Analysis.py

#==============================================================================

import os, sys, glob, shutil
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.colors as colors
from matplotlib.colors import BoundaryNorm
from scipy.interpolate import RBFInterpolator

from .commonFunctions import *

from pBabel                    import *                                     
from pCore                     import *                                     
from pMolecule                 import *                   
from pScientific               import *                                                             
from pScientific.Statistics    import *
from pScientific.Arrays        import *
from pSimulation               import *

#*********************************************************************
class EnergyAnalysis:
	"""Manage and visualize energy data from molecular simulations.
	
	This class handles reading, storing, and plotting 1D/2D energy surfaces from
	simulation data. It supports multiple output formats and can generate contour
	plots, free energy curves, and PMF surfaces using WHAM data.
	
	Attributes:
		energies1D (list): 1D energy values.
		energiesMatrix (ndarray): 2D energy grid.
		RC1, RC2 (list): Reaction coordinate values.
		xlen, ylen (int): Dimensions of energy grid.
		Type (str): Plot type identifier.
		dimensions (int): Number of coordinates (1 or 2).
	"""
	#------------------------------------------------
	def __init__(self, x, y, _type="1D"):
		"""Initialize EnergyAnalysis object.
		
		Args:
			x (int): Number of steps in first coordinate.
			y (int): Number of steps in second coordinate (0 for 1D).
			_type (str, optional): Plot type. Options: '1D', '2D', 'WHAM1D', 'WHAM2D',
			    'FE1D', 'FE2D', 'PMF', etc. Defaults to '1D'.
		"""
		self.energies1D 	= []						  # array for energy values from one-dimension coordinate simulation 
		self.energiesMatrix = np.zeros( (y, x), dtype=float ) # array for energy values from two-dimension coordinate simulation
		self.multiple1Dplot = []							  # List of one-dimension energy arrays, each one for a different energy method
		self.multiple2Dplot = []							  # List of two-dimension energy arrays, 
		self.RC1            = []							  # Array with the first reaction coordinate values
		self.RC2            = []							  # Array with the second reaction coordinate values
		self.dimensions     = 0								  # Number of coordinates used in the analysis
		self.nplots1D       = 0								  # Number of one-dimension energy plots to do 
		self.nplots2D 		= 0								  # Number of two-dimension energy plots to do
		self.xlen 			= x 							  # Number of steps in the first coordinate
		self.ylen  			= y         					  # Number of steps in the seconf coordinate
		self.Type 			= _type                           # Type of the plot, could be of: Free Energy, PMF, Potential Energy 
		self.labely         = ""                              # String holding the name of Y axis  
		self.baseName 		= ""							  # string holding the path of the folder
		self.identifiers    = []                              # List of string identifiers 
		self.fig_size_x     = 0 
		self.fig_size_y     = 0 
		self.resumedTrajectory = None
		self.plot1d_name   = "1d.png"
		self.plot2d_name   = "2d.png"
		if self.ylen > 0:
			self.dimensions = 2
		else:
			self.dimensions = 1

		
	#================================================
	def ReadLog(self, _fileName):
		"""Parse energy data from log files.
		
		Reads and processes energy log files in various formats (1D, 2D, WHAM, FE).
		Automatically detects and handles different data organization schemes.
		
		Args:
			_fileName (str): Path to energy log file.
			
		Sets:
			RC1, RC2 (list): Reaction coordinate values from log.
			energies1D, energiesMatrix: Parsed energy data.
		"""
		self.baseName = _fileName[:-4]
		reading       = open(_fileName,'r')
		i = 0 
		energyTmp = []
		#----------------------------------
		if self.Type == "1D":
			self.energies1D = []
			for line in reading:
				if i > 0:
					lns = line.split()
					self.RC1.append( float(lns[1] ) )
					energyTmp.append( float(lns[2] ) )
					self.energies1D.append( float(lns[2]) )
				i += 1
			self.multiple1Dplot = None
			self.labely = "Potential Energy (kJ/mol)"
		#-----------------------------------
		elif self.Type == "1DRef":
			i = 0
			oldMethod = "none"
			method    = "none"
			for line in reading:
				if i > 0:
					lns = line.split()
					if oldMethod == "none":
						oldMethod = lns[3]
					method = lns[3]
					if not method == oldMethod:
						self.multiple1Dplot.append(energyTmp)
						self.identifiers.append(oldMethod)
						oldMethod = method
						self.nplots1D += 1	
						energyTmp = []		
					energyTmp.append( float(lns[1]) )
					self.energies1D.append( float(lns[1]) )
				i+=1
			self.multiple1Dplot.append(energyTmp)
			self.identifiers.append(method)
			self.labely = "Potential Energy (kJ/mol)"
			self.plot1d_name = "1D_Ref.png"
		#----------------------------------
		elif self.Type == "2D":
			records = []
			max_m = -1
			max_n = -1
			for line in reading:
				if i > 0:
					lns = line.split()
					m = int(lns[0])
					n = int(lns[1])
					max_m = max(max_m, m)
					max_n = max(max_n, n)
					self.RC1.append( float(lns[2] ) )
					self.RC2.append( float(lns[3] ) )
					records.append((m, n, float(lns[4])))
				i += 1
			if self.xlen <= 0 or self.ylen <= 0 or self.xlen <= max_m or self.ylen <= max_n:
				self.xlen = max_m + 1
				self.ylen = max_n + 1
				self.energiesMatrix = np.zeros((self.ylen, self.xlen), dtype=float)
			elif self.energiesMatrix is None:
				self.energiesMatrix = np.zeros((self.ylen, self.xlen), dtype=float)
			for m, n, value in records:
				self.energiesMatrix[n][m] = value

		#----------------------------------
		elif self.Type == "2DRef":
			oldMethod = "none"
			method    = "none"
			i = 0 
			for line in reading:
				if i > 0:
					lns = line.split()
					if oldMethod == "none":
						oldMethod = lns[4]
					method = lns[4]
					if not method == oldMethod:
						self.multiple2Dplot.append(self.energiesMatrix)
						self.identifiers.append(oldMethod)
						oldMethod = method
						self.nplots2D += 1
						self.energiesMatrix = np.zeros( (self.ylen, self.xlen), dtype=float )
					m = int(lns[0])		
					n = int(lns[1])				
					self.energiesMatrix[n][m] = float(lns[2])
				i += 1
			
			self.multiple2Dplot.append(self.energiesMatrix)
			self.identifiers.append(method)
			self.nplots2D += 1
		#----------------------------------
		elif self.Type == "WHAM1D":
			MaX = 0.0
			for line in reading:
				lns = line.split()
				pmf = float(lns[1])
				if pmf > MaX: MaX = pmf
				self.RC1.append( float(lns[0]) )
				if lns[1] == "inf": self.energies1D.append( 43434.0000 )
				else:       		self.energies1D.append( float(lns[1]) )

			for i in range(len(self.energies1D)):
				if self.energies1D[i] == 43434.0000:
					self.energies1D[i] = MaX
			self.plot1d_name = "WHAM1D.png"
		#----------------------------------
		elif self.Type == "WHAM2D":
			m = 0
			n = 0
			MaX = 0.0
			for line in reading:
				lns = line.split()				
				self.RC1.append( float(lns[0]) )
				self.RC2.append( float(lns[1]) )
				if lns[2] == "inf":
					self.energiesMatrix[m][n] = 43434.0000
				else:
					self.energiesMatrix[m][n] = float(lns[2])	
					pmf = float(lns[2])
					if pmf > MaX:
						MaX = pmf			
				i +=1
				n +=1 				
				if i % self.xlen == 0:
					m += 1
					n = 0	
			for j in range(self.xlen):
				for i in range(self.ylen):
					if self.energiesMatrix[i][j] == 43434.0000:
						self.energiesMatrix[i][j] = MaX
			self.plot2d_name = "WHAM2D.png"
		#----------------------------------
		elif self.Type == "FE1D":
			energyTmp = np.zeros( (self.xlen), dtype=float )
			for line in reading:
				lns = line.split()
				m = int( lns[0] )
				energyTmp[m]  = float(lns[1]) 
			self.energies1D = energyTmp
			self.plot1d_name = "FE1D.png"		
		#----------------------------------
		elif self.Type == "FE2D":
			for line in reading:
				lns = line.split()
				m = int( lns[0])				
				n = int( lns[1])
				self.energiesMatrix[n][m] = float(lns[2])
			self.plot2d_name = "FE2D.png"
		#----------------------------------
		self.nplots1D += 1	
	#================================================
	def ReadLogs(self,_folder):
		"""Read multiple energy log files from a directory.
		
		Collects and processes multiple log files from a folder for batch processing.
		
		Args:
			_folder (str): Path to folder containing log files.
			
		Note:
			Log files should have .log extension.
		"""
		_path = os.path.join(_folder,"")
		logs = glob.glob( _path + "*.log" )
		for log in logs:
			self.ReadLog(log)
			self.identifiers.append( os.path.basename(log[:-4]) )
	#==============================================
	def NormalizeEnergies(self):
		'''
		Normalize energy arrays
		'''
		#------------------------------------------
		if self.Type == "1D":
			Min = 0
			if self.nplots1D == 1:
				Min = self.energies1D[0]
				for i in range( len(self.energies1D) ):
					self.energies1D[i] = self.energies1D[i] - Min
			elif self.nplots1D > 2:
				for k in range( self.nplots1D ):
					Min = self.multiple1Dplot[k][0]
					for i in range(len(self.multiple1Dplot)):
						self.multiple1Dplot[k][i] = self.multiple1Dplot[k][i] - Min		
		#------------------------------------------
		if self.Type == "2D" or self.Type == "WHAM2D" or self.Type == "FE2D" or self.Type == "2DRef":
			if not self.energiesMatrix[0][0] == 0.0:
				self.energiesMatrix = self.energiesMatrix - np.min(self.energiesMatrix)
	#===============================================
	def FES_HL_SMO(self, logPES, logSMO, logFE):
		'''
		Free energy surface from a combination of High level QC method PES and semiempirical free energy
		Parameters:
			logPES:
			logSMO:
			logFE:
		'''
		pass
	#===============================================
	def Plot1D(self,label,XLIM=None,SHOW=False):
		'''
		Plot one dimensional energy plot.
		'''

		self.NormalizeEnergies()
		if self.Type == "1DRef":
			if XLIM == None: self.RC1 = np.linspace( 0,len(self.energies1D),len(self.energies1D) )
			else:	         self.RC1 = np.linspace( XLIM[0],XLIM[1],len(self.energies1D) )
			self.labely = "Potential Energy (kJ/mol)"	
		elif self.Type == "FE1D":
			if XLIM == None: self.RC1 = np.linspace( 0,len(self.energies1D),len(self.energies1D) )
			else:	         self.RC1 = np.linspace( XLIM[0],XLIM[1],len(self.energies1D) )
			self.labely = "Free Energy (kJ/mol)"
		elif self.Type == "WHAM1D":
			self.RC1 = np.linspace( np.min(self.RC1), np.max(self.RC1), len(self.RC1) )
			self.labely = "Potential of Mean Field (kJ/mol)"
		
		#--------------------------------------------
		plt.plot(self.RC1,self.energies1D,'-ok')
		if self.fig_size_x > 0:
			plt.set_size_inches(self.fig_size_x,self.fig_size_y)
		plt.xlabel(label)
		plt.ylabel(self.labely)	
		#--------------------------------------------
		plt.savefig(self.baseName+self.plot1d_name,dpi=1000)
		#---------------------------------------------
		if SHOW: plt.show()
		plt.clf()
		plt.close()


	#===============================================
	def MultPlot1D(self,label,SHOW=False):
		'''
		Plot one-dimensinal energy plot for several methods
		'''
		#---------------------------------------------
		self.NormalizeEnergies()
		x = np.linspace(0, self.xlen, self.xlen )
		for i in range(self.nplots1D):
			plt.plot(x,self.multiple1Dplot[i],label=self.identifiers[i])
		#---------------------------------------------
		plt.xlabel(label)
		plt.ylabel(self.labely)
		plt.legend()
		plt.savefig(self.baseName+".png",dpi=1000)
		#---------------------------------------------
		if SHOW: plt.show()	
		plt.clf()	
		plt.close()
	#==================================================================
	def Plot2D(self,contourlines,
			crd1label,
			crd2label,
			_xlim=None,_ylim=None,
			SHOW=False,
			_figS=[7,5],
			_reverserc1=False,
			_reverserc2=False,
			path_points=None):
		'''
		Plot contour plot for potential, free energy and potential of mean field.
		
		Args:
			contourlines: Number of contour lines to plot
			crd1label: Label for x-axis (first coordinate)
			crd2label: Label for y-axis (second coordinate)
			_xlim: X-axis limits [min, max]
			_ylim: Y-axis limits [min, max]
			SHOW: Show plot window
			_figS: Figure size [width, height]
			_reverserc1: Reverse x-axis
			_reverserc2: Reverse y-axis
			path_points: Tuple of (pathx, pathy) lists for plotting minimum energy path
		'''			
		#-----------------------------------------------------
		X = []
		Y = []
		self.NormalizeEnergies()
		# Force xlen and ylen to be integers if they're not
		if not isinstance(self.xlen, int):
			print(f"DEBUG: xlen is not int, converting. Original value: {self.xlen}")
			self.xlen = int(self.xlen)
		if not isinstance(self.ylen, int):
			print(f"DEBUG: ylen is not int, converting. Original value: {self.ylen}")
			self.ylen = int(self.ylen)
		#------------------------------------------------------
		if _xlim == None:
			_xlim = [ 0, self.xlen ]
			if len(self.RC1) > 0:
				X = np.linspace( np.min(self.RC1) , np.max(self.RC1), self.xlen )
			else:
				if _reverserc1: 
					X = np.linspace(self.xlen,0,self.xlen)
				else:
					X = np.linspace(0,self.xlen,self.xlen)				
		else:			
			X = np.linspace(_xlim[0],_xlim[1],self.xlen)
		#------------------------------------------------------

		if _ylim == None:
			_ylim = [ 0, self.ylen ]
			if len(self.RC2) > 0:
				Y = np.linspace( np.min(self.RC2) , np.max(self.RC2), self.ylen )
			else:
				if _reverserc2: 
					Y = np.linspace(self.ylen,0,self.ylen)
				else: 
					Y = np.linspace(0,self.ylen,self.ylen)
		#------------------------------------------------------
		else:
			Y = np.linspace(_ylim[0],_ylim[1],self.ylen)
		#------------------------------------------------------		
		z = self.energiesMatrix
		#------------------------------------------------------
		fig, (ax0) = plt.subplots( nrows=1, figsize=(_figS[0],_figS[1]) )
		vmin=z.min()
		vmax=z.max()
		#------------------------------------------------------
		levels = MaxNLocator(nbins=20).tick_values( z.min(), z.max() )
		cmap = plt.get_cmap("jet")
		#------------------------------------------------------
		norm = BoundaryNorm(levels, ncolors=cmap.N,	clip=True)
		norm= colors.PowerNorm(gamma=1./2.)
		norm= colors.Normalize(vmin=vmin, vmax=vmax)
		#------------------------------------------------------
		X_grid, Y_grid = np.meshgrid(X, Y)
		im = ax0.pcolormesh(X,Y,z, cmap=cmap, norm=norm, shading = "gouraud")		
		am = ax0.contour(X_grid,Y_grid,z,contourlines, colors='k')		
		ax0.clabel(am, inline=1, fontsize=8, fmt='%1.1f',colors="k")
		
		# Plot minimum energy path if provided
		if path_points is not None:
			pathx, pathy = path_points
			# Convert grid indices to coordinate values
			if len(self.RC1) > 0 and len(self.RC2) > 0:
				# Map grid indices to actual coordinate values
				path_x_coords = [X[min(int(x), len(X)-1)] for x in pathx]
				path_y_coords = [Y[min(int(y), len(Y)-1)] for y in pathy]
			else:
				path_x_coords = pathx
				path_y_coords = pathy
			
			# Plot path as line with markers
			ax0.plot(path_x_coords, path_y_coords, 'r-', linewidth=2, alpha=0.7, label='Minimum Energy Path')
			ax0.plot(path_x_coords, path_y_coords, 'ro', markersize=6, alpha=0.8)
			# Mark start and end points
			ax0.plot(path_x_coords[0], path_y_coords[0], 'g^', markersize=10, label='Start')
			ax0.plot(path_x_coords[-1], path_y_coords[-1], 'bs', markersize=10, label='End')
			ax0.legend(loc='best', fontsize=10)
		
		cbar = fig.colorbar(im, ax=ax0)
		cbar.ax.tick_params()
		#---------------------------------------------
		# Set the tick labels font
		#axis_font = {'fontname':'Michroma', 'size':14}
		for tick in (ax0.xaxis.get_major_ticks()):
			#tick.label1.set_fontname('Arial')
			tick.label1.set_fontsize(14)
		for tick in (ax0.yaxis.get_major_ticks()):
			#tick.label1.set_fontname('Dejavu')
			tick.label1.set_fontsize(14) 
		#---------------------------------------------				
		ax0.set_xlabel(crd1label)
		ax0.set_ylabel(crd2label)
		fig.tight_layout()
		_method = ""
		if len(self.identifiers) > 0: 
			_method = "_" + self.identifiers[-1]

		plotName = self.baseName + _method[:5]
		if path_points is not None:
			plotName += "_MEP"
		plt.savefig(plotName+".png",dpi=1000)
		if SHOW: plt.show()
		plt.clf()
		plt.close()
		fig.clf()
	#----------------------------------------------------------------------------------------
	def MultPlot2D(self,contourlines,
				crd1label,
				crd2label,
				_xlim=None,
				_ylim=None,
				SHOW=False,
				_reverserc1=False,
				_reverserc2=False):
		'''
		'''
		for i in range(self.nplots2D):
			self.identifiers.append( self.identifiers[i] )
			self.energiesMatrix = self.multiple2Dplot[i]
			self.Plot2D(contourlines,crd1label,crd2label,_xlim=_xlim,_ylim=_ylim,SHOW=False,_reverserc1=_reverserc1,_reverserc2=_reverserc2)
	#----------------------------------------------------------------------------------------
	def Plot1D_FreeEnergy(self,crd1label,crd2label,SHOW=False):
		'''
		'''
		self.NormalizeEnergies()
		if 	 self.Type == "FE2D"   or self.Type == "FE1D"  : self.labely = "Free Energy (kJ/mol)"
		elif self.Type == "WHAM1D" or self.Type == "WHAM2D": self.labely = "Potential of Mean Field (kJ/mol)"
		
		rc0 = np.linspace( 0, len(self.energies1D)-1, len(self.energies1D) )
		#--------------------------------------------
		plt.plot(rc0,self.energies1D,'-ok')
		plt.xlabel("Frame Window (n)")
		plt.ylabel(self.labely)	
		plt.savefig(self.baseName+self.plot1d_name,dpi=1000)
		#---------------------------------------------
		if SHOW: plt.show()	
		plt.clf()
		plt.close()
	
	#---------------------------------------------------------------------------------------
	def Calulate_Kcat(self, energy_barrier, T=300, kT=0.593):
		'''
		Calculate kcat from the free energy barrier using Eyring equation
		Parameters:
			energy_barrier: Free energy barrier in kJ/mol
			T: Temperature in Kelvin
			kT: Boltzmann constant in kcal/mol
		Returns:
			kcat in s^-1
		'''
		R = 8.31446261815324 # J/(mol*K)
		kB = 1.380649e-23 # J/K
		h = 6.62607015e-34 # J*s
		energy_barrier_J = energy_barrier * 1000  # Convert kJ/mol to J/mol
		prefactor = (kB * T) / h # in s^-1
		exponent = - (energy_barrier_J ) / (R * T)
		kcat = prefactor * math.exp(exponent)
		return kcat

	#----------------------------------------------------------------------------------------
	def Path_From_PES(self, in_point,fin_point,_path_pkl,_folder_dst,_system,min_points=15,max_points=21,method="hessian"):
		'''
		Calculate minimum energy path from starting point to final point on PES.
		
		Returns:
			tuple: (pathx, pathy) - lists of x and y coordinates of the minimum energy path
		'''
		#setting current point as initial

		z = self.energiesMatrix	
		dirs = [ [1,0] ,[0,1], [1,1] ]
		pathx = [in_point[0]] 
		pathy = [in_point[1]]
		self.energies1D.append( z[ pathx[0],pathy[0] ] )
		cp = [ in_point[0],in_point[1] ]
		final_fp = [ fin_point[0], fin_point[1] ]

		if not os.path.exists( os.path.join(_folder_dst,"traj1d.ptGeo") ): 
			os.makedirs( os.path.join(_folder_dst,"traj1d.ptGeo") )

		if method == "Dijkstra":
			print(cp, fin_point)
			paths = self.FindMinEnergyPath_Dijkstra(in_point,fin_point)
			for i in range(1,len(paths[0])): 
				pathx.append(paths[0][i][0])
				pathy.append(paths[0][i][1])
		elif method =="MiniMax":
			print(cp, fin_point)
			paths = self.FindMinimaxPath(in_point,fin_point)
			for i in range(1,len(paths[0])): 
				pathx.append(paths[0][i][0])
				pathy.append(paths[0][i][1])
		elif method == "hessian" or method =="simpleMEP":
			if method == "hessian":
				Results = self.AnalyzePES2D( output_folder= os.path.join(_folder_dst,"path_analysis"), in_point=in_point, fin_point=fin_point )
				saddle_p = Results["saddle_points"]
				reactants = Results["reactants"]
				products  = Results["products"]
				cp = [ reactants[0], reactants[1] ]
				final_fp  = [ products[0], products[1] ]			
				if len(saddle_p) > 0:
					fin_point = [ saddle_p[0][0], saddle_p[0][1] ]
		
			print(cp, fin_point)
			#search path until the TS
			while not cp == fin_point:

				print( "Current Point is: {} {} ".format(cp[0], cp[1] ) )
				print( "Energy of the Current Point: {}".format( z[ cp[1],cp[0] ] ) )

				A = math.inf
				B = math.inf
				C = math.inf			

				if ( cp[0] + 1 ) < self.xlen: 
					if  (cp[0] + 1) <= fin_point[0]:
						A = z[ cp[1], cp[0] ] + z[ cp[1], (cp[0] + 1) ]
						print( "Increment in X:  {}".format(A) )
				else: print( "No Increment in X:  {}".format(A) )
				if ( cp[1] + 1 ) < self.ylen:
					if  (cp[1] + 1) <= fin_point[1] : 
						B = z[ cp[1], cp[0] ] + z[ (cp[1] + 1), cp[0] ] 
						print( "Increment in Y:  {}".format(B) )
				else: print( "No Increment in Y: {}".format(B) )

				if ( cp[0] + 1 ) < self.xlen and ( cp[1] + 1 ) < self.ylen: 
					if  (cp[0] + 1) <=fin_point[0] and (cp[1] + 1) <= fin_point[1]:
						C = z[ cp[1], cp[0] ] + z[ (cp[1] + 1), (cp[0] + 1) ] 
						print( "Increment in both directions:  {}".format(C) )
				else: print("No Increment in both directions: {}".format(C) )

				D = [ A, B, C ]
				ind = D.index(min(D))			
				cp[0] += dirs[ind][0]
				cp[1] += dirs[ind][1]			
				self.energies1D.append( z[ cp[1], cp[0] ] )
				pathx.append(cp[0])
				pathy.append(cp[1])

			if not fin_point == final_fp:
				print(final_fp)
				print("Going from the sadle to the final point")
				while not cp == final_fp:

					print( "Current Point is: {} {} ".format(cp[0], cp[1] ) )
					print( "Energy of the Current Point: {}".format( z[ cp[1],cp[0] ] ) )

					A = math.inf
					B = math.inf
					C = math.inf			

					if ( cp[0] + 1 ) < self.xlen: 
						if  (cp[0] + 1) <= final_fp[0]:
							A = z[ cp[1], cp[0] ] + z[ cp[1], (cp[0] + 1) ]
							print( "Increment in X:  {}".format(A) )
					else: print( "No Increment in X:  {}".format(A) )
					if ( cp[1] + 1 ) < self.ylen:
						if  (cp[1] + 1) <= final_fp[1] : 
							B = z[ cp[1], cp[0] ] + z[ (cp[1] + 1), cp[0] ] 
							print( "Increment in Y:  {}".format(B) )
					else: print( "No Incrementin Y: {}".format(B) )

					if ( cp[0] + 1 ) < self.xlen and ( cp[1] + 1 ) < self.ylen: 
						if  (cp[0] + 1) <=final_fp[0] and (cp[1] + 1) <= final_fp[1]:
							C = z[ cp[1], cp[0] ] + z[ (cp[1] + 1), (cp[0] + 1) ] 
							print( "Increment in both directions:  {}".format(C) )
					else: print("No Increment in both directions: {}".format(C) )

					D = [ A, B, C ]
					ind = D.index(min(D))			
					cp[0] += dirs[ind][0]
					cp[1] += dirs[ind][1]			
					self.energies1D.append( z[ cp[1], cp[0] ] )
					pathx.append(cp[0])
					pathy.append(cp[1])
		
		#---------------------------------------------------------------
		
		kcats = [0.0]
		new_idx = 0
		min_energy = self.energies1D[0]

		for indx in range( 1,len(pathx) ):			
			pkl = _path_pkl + "/frame{}_{}.pkl".format(pathx[indx],pathy[indx])			
			finalPath = os.path.join( _folder_dst , "traj1d.ptGeo/frame{}.pkl".format(new_idx) )			
			_system.coordinates3 = ImportCoordinates3(pkl,log=None)
			pdb_file = os.path.join( _folder_dst , "frame{}.pdb".format(new_idx) )
			ExportSystem( pdb_file,_system,log=None)
			shutil.copy(pkl,finalPath)
			new_idx +=1
			if self.energies1D[indx] < min_energy:
				min_energy = self.energies1D[indx]				
			if indx > 0 and indx < (len(pathx)-2):
				if (self.energies1D[indx] > (self.energies1D[indx-1] + 1.0) ) and \
					(self.energies1D[indx] >  (self.energies1D[indx+1]+ 1.0 ) ):
					print("Found potential barrier at frame {} with energy {}".format(indx, self.energies1D[indx] ) )
					kcats.append( self.Calulate_Kcat( self.energies1D[indx] - min_energy ) )
					print("Calculated kcat for this barrier is: {} s^-1".format(kcats[-1]) )
				else:
					kcats.append(0.0)
			else:
				kcats.append(0.0)	

		#---------------------------------------------------------------
		# Export the path as a trajectory and log file for visualization in PyMOL and other tools
		#---------------------------------------------------------------
		trajName = os.path.join( _folder_dst, "traj1d.dcd" )
		trajpath = os.path.join( _folder_dst, "traj1d.ptGeo" )
		try: Duplicate( trajpath, trajName, _system ) 
		except: pass
		log_text = "x Energy method Energy_kcal kcat \n"
		new_log  = open( os.path.join(_folder_dst,"traj1D.log"), 'w' )
		for i in range(len(self.energies1D)):						
			energy_kcal = self.energies1D[i] / 4.184
			log_text += "{} {} pickPath {} {}\n".format(i,self.energies1D[i],energy_kcal,kcats[i])
		new_log.write(log_text)
		self.Type = "1DRef"
		self.Plot1D("Reaction Path frames (n)")

		pymol_text = "preset.publication(selection='all')\n"
		pymol_text+= "set sticks\n"
		pymol_text+= "set label_size, 20\n"
		pymol_text+= "set sphere_scale, 0.2\n"
		pymol_text+= "set bg_rgb, white\n" 
		pymol_text+= "set stick_radius, 0.18\n"
		pymol_text+= "load {}".format( "frame0.pdb" )
		pymol_text+= "\nload_traj {}, ".format( "traj1d.dcd" )
		pymol_text+= "frame0, 1, start=1, stop=-1, interval=1"

		pymols_file = open( os.path.join(_folder_dst,"traj1d.pym"), "w") 
		pymols_file.write(pymol_text)
		pymols_file.close()

		#-----------------------------------
		#Resample the path to reduce number of frames while preserving key features
		#-----------------------------------
		pathIndices = list(zip(pathx, pathy))
		Kept_indices = self.ResamplePath(pathIndices, min_points=min_points, max_points=max_points, include_curvature=True)

		if not os.path.exists( os.path.join(_folder_dst,"traj1d_resampled.ptGeo") ): 
			os.makedirs( os.path.join(_folder_dst,"traj1d_resampled.ptGeo") )

		#write the resampled path log
		log_text = "x Energy method Energy_kcal kcat \n"
		new_log  = open( os.path.join(_folder_dst,"traj1D_resampled.log"), 'w' )
		for i in Kept_indices:						
			energy_kcal = self.energies1D[i] / 4.184
			log_text += "{} {} pickPath {} {}\n".format(i,self.energies1D[i],energy_kcal,kcats[i])
		new_log.write(log_text)
		new_log.close()

		#export resampled path frames
		self.energies1D = [self.energies1D[i] for i in Kept_indices]
		self.Type = "1DRef"
		self.plot1d_name = "1D_Ref_resampled.png"
		self.Plot1D("Reaction Path frames (n)")
		for idx in range( len(Kept_indices) ):
			pkl = _path_pkl + "/frame{}_{}.pkl".format(pathx[idx],pathy[idx])			
			finalPath = os.path.join( _folder_dst , "traj1d_resampled.ptGeo/frame{}.pkl".format(idx) )			
			_system.coordinates3 = ImportCoordinates3(pkl,log=None)
			pdb_file = os.path.join( _folder_dst , "frame{}_rsd.pdb".format(idx) )
			ExportSystem( pdb_file,_system,log=None)
			shutil.copy(pkl,finalPath)

		trajName = os.path.join( _folder_dst, "traj1d_resampled.dcd" )
		trajpath = os.path.join( _folder_dst, "traj1d_resampled.ptGeo" )
		try: Duplicate( trajpath, trajName, _system ) 
		except: pass
		
		pymol_text = "preset.publication(selection='all')\n"
		pymol_text+= "set sticks\n"
		pymol_text+= "set label_size, 20\n"
		pymol_text+= "set sphere_scale, 0.2\n"
		pymol_text+= "set bg_rgb, white\n" 
		pymol_text+= "set stick_radius, 0.18\n"
		pymol_text+= "load {}".format( "frame0.pdb" )
		pymol_text+= "\nload_traj {}, ".format( "traj1d_resampled.dcd" )
		pymol_text+= "frame0, 1, start=1, stop=-1, interval=1"


		pymols_file = open( os.path.join(_folder_dst,"traj1d_resampled.pym"), "w") 
		pymols_file.write(pymol_text)
		pymols_file.close()

		# Return the minimum energy path coordinates
		return pathx, pathy
	#----------------------------------------------------------------------------------------
	def ResamplePath(self, path_indices,min_points=15, max_points=20, include_curvature=True):
		"""
		Reduce number of points in a reaction path while preserving key features.
    
    	Args:
        	path_indices (list of tuples): (x, y) grid indices along path.
        	max_points (int): Maximum number of points to keep.
        	include_curvature (bool): Keep points where path changes direction.
    
    	Returns:
        	list of indices into the original path to keep.
    	"""
		n = len(self.energies1D)
		keep = set()
		keep.add(0)                 # first point
		keep.add(n-1)               # last point
    
    	# 1. Local extrema (minima and maxima)
		for i in range(1, n-1):
			e_prev, e_curr, e_next = self.energies1D[i-1], self.energies1D[i], self.energies1D[i+1]
			if (e_curr > e_prev and e_curr > e_next) or (e_curr < e_prev and e_curr < e_next):
				keep.add(i)
    
		# 2. Direction changes (curvature)
		if include_curvature and n > 3:
			for i in range(1, n-1):
				v1 = np.array(path_indices[i]) - np.array(path_indices[i-1])
				v2 = np.array(path_indices[i+1]) - np.array(path_indices[i])
				n1 = np.linalg.norm(v1)
				n2 = np.linalg.norm(v2)
				if n1 > 0 and n2 > 0:
					cos_angle = np.dot(v1, v2) / (n1 * n2)
					if cos_angle < 0.707:   # angle > 45°
						keep.add(i)

		kept = sorted(keep)

		# 3. If mandatory set is smaller than min_points, add evenly spaced optional points
		optional = []
		if len(kept) < min_points:
		# Always keep the original first and last
			needed = min_points - len(kept)
			# Create a set of additional indices to add, choosing from those not already kept
			candidates = set(range(n)) - keep
			if candidates:
				# Sort candidates (they are naturally in order)
				candidates = sorted(candidates)
				# Choose evenly spaced indices from the candidate list
				step = len(candidates) / (needed + 1)
				optional = [candidates[int(i*step)] for i in range(1, needed+1)]

		all_kept = sorted(keep.union(optional))
    
    	# 3. If still too many, subsample evenly
		if len(all_kept) > max_points:
			max_optional = max_points - len(keep)
			if max_optional <=0:
				return all_kept
			else:
				if len(optional) > max_optional:
					step = len(optional) / max_optional
					subsampled_optional = [optional[int(i*step)] for i in range(max_optional)]
				else:
					subsampled_optional = optional
				all_kept = sorted(keep.union(subsampled_optional))

		return all_kept 
	#----------------------------------------------------------------------------------------
	def AnalyzePES2D(self, energy_threshold=0.5, output_folder=None, in_point=None, fin_point=None):
		"""
		Analyze 2D potential energy surface to identify critical points.
		
		Detects and classifies critical points on the surface:
		- Local minima: Points with lower energy than all neighbors
		- Saddle points: Points where Hessian has both positive and negative eigenvalues
		- Intermediates: Points between reactant/product regions with moderate energy
		- Reactants/Products: Nearest minima to initial/final points (if provided)
		
		Args:
			energy_threshold (float): Energy difference threshold for intermediate detection (kcal/mol)
			output_folder (str): Folder to save analysis results. If None, uses self.baseName
			in_point (tuple): (x, y) coordinates of initial/reactant region. If provided, finds nearest minimum
			fin_point (tuple): (x, y) coordinates of final/product region. If provided, finds nearest minimum
			
		Returns:
			dict: Dictionary containing:
				- 'minima': list of (x, y, energy) tuples for local minima
				- 'saddle_points': list of (x, y, energy, eigenvalues) tuples sorted by energy
				- 'intermediates': list of (x, y, energy) tuples
				- 'reactants': (x, y, energy) tuple of nearest minimum to in_point (if provided)
				- 'products': (x, y, energy) tuple of nearest minimum to fin_point (if provided)
				- 'analysis_report': string with formatted report
		"""
		if self.dimensions != 2:
			raise ValueError("Two-dimensional energy surface required for PES analysis")
			
		z = self.energiesMatrix
		results = {
			'minima': [],
			'saddle_points': [],
			'intermediates': [],
			'gradient_points': [],
			'inflection_points': [],
			'reactants': None,
			'products': None
		}
		
		# Ensure output folder exists
		if output_folder is None:
			output_folder = self.baseName
		if not os.path.exists(output_folder):
			os.makedirs(output_folder)
		
		# 1. Find local minima (8-neighbor connectivity)
		print("\n" + "="*60)
		print("ANALYZING 2D POTENTIAL ENERGY SURFACE")
		print("="*60)
		print("\n[1] Searching for local minima...")
		
		minima = []
		for i in range(1, self.ylen - 1):
			for j in range(1, self.xlen - 1):
				center = z[i, j]
				# Check 8 neighbors
				neighbors = [
					z[i-1, j-1], z[i-1, j], z[i-1, j+1],
					z[i, j-1],                z[i, j+1],
					z[i+1, j-1], z[i+1, j], z[i+1, j+1]
				]
				if center < min(neighbors):
					minima.append((j, i, center))
					results['minima'].append((j, i, center))
		
		print(f"  Found {len(minima)} local minima")
		for x, y, e in minima[:10]:  # Print first 10
			print(f"    - Position ({x:3d}, {y:3d}): Energy = {e:10.4f} kJ/mol")
		if len(minima) > 10:
			print(f"    ... and {len(minima)-10} more")
		
		# Find reactants and products if initial/final points are provided
		if in_point is not None and len(minima) > 0:
			print(f"\n  Finding nearest minimum to initial point {in_point}...")
			nearest_reactant = min(minima, key=lambda m: (m[0]-in_point[0])**2 + (m[1]-in_point[1])**2)
			results['reactants'] = nearest_reactant
			dist_to_reactant = np.sqrt((nearest_reactant[0]-in_point[0])**2 + (nearest_reactant[1]-in_point[1])**2)
			if dist_to_reactant > 3.0:
				print(f"    → Products: ({in_point[0]}, {in_point[1]}), Energy = {z[in_point[1],in_point[0]]:.6f} kJ, Distance = {dist_to_reactant:.2f}")
				results['reactants'] = ( in_point[0], in_point[1], z[in_point[1],in_point[0]] )
			else:
				print(f"    → Reactants: ({nearest_reactant[0]}, {nearest_reactant[1]}), Energy = {nearest_reactant[2]:.6f} kJ, Distance = {dist_to_reactant:.2f}")
		
			
			print(f"    → Reactants: ({nearest_reactant[0]}, {nearest_reactant[1]}), Energy = {nearest_reactant[2]:.6f} kJ, Distance = {dist_to_reactant:.2f}")
		
		if fin_point is not None and len(minima) > 0:
			print(f"\n  Finding nearest minimum to final point {fin_point}...")
			nearest_product = min(minima, key=lambda m: (m[0]-fin_point[0])**2 + (m[1]-fin_point[1])**2)
			results['products'] = nearest_product
			dist_to_product = np.sqrt((nearest_product[0]-fin_point[0])**2 + (nearest_product[1]-fin_point[1])**2)
			if dist_to_product > 3.0:
				print(f"    → Products: ({fin_point[0]}, {fin_point[1]}), Energy = {z[fin_point[1],fin_point[0]]:.6f} kJ, Distance = {dist_to_product:.2f}")
				results['products'] = ( fin_point[0], fin_point[1], z[fin_point[1],fin_point[0]] )
			else:
				print(f"    → Products: ({nearest_product[0]}, {nearest_product[1]}), Energy = {nearest_product[2]:.6f} kJ, Distance = {dist_to_product:.2f}")
		
		
		
		# First pass: identify candidates with strong saddle character
		saddle_candidates = []
		
		from scipy.interpolate import RBFInterpolator
		# Collect all points and energies
		points = np.column_stack((self.RC1, self.RC2))   # (N, 2)
		values = self.energiesMatrix.flatten()           # energies
		rbf = RBFInterpolator(points, values, kernel='thin_plate_spline', smoothing=0.001)
		
		# Define a uniform grid covering the range of RC1 and RC2
		x_uniform = np.linspace(min(self.RC1), max(self.RC1), self.xlen)
		y_uniform = np.linspace(min(self.RC2), max(self.RC2), self.ylen)
		X, Y = np.meshgrid(x_uniform, y_uniform, indexing='xy')
		Xg, Yg = np.meshgrid(x_uniform, y_uniform, indexing='xy')
		grid_points = np.column_stack((Xg.ravel(), Yg.ravel()))
		# Interpolate (use 'cubic' for smoothness, 'linear' for speed)
		z_smooth = rbf(grid_points).reshape(self.ylen, self.xlen)
		# Now Z_uniform is on a perfect rectangular grid with constant dx, dy
		dx = x_uniform[1] - x_uniform[0]
		dy = y_uniform[1] - y_uniform[0]

		for i in range(1, self.ylen - 1):
			for j in range(1, self.xlen - 1):
				# Get local spacings (central differences use average of surrounding spacings)
										
				# Approximate Hessian using finite differences on SMOOTHED energy with PHYSICAL spacing
				# H_xx = (E(i,j+1) - 2*E(i,j) + E(i,j-1)) / dx²
				# H_yy = (E(i+1,j) - 2*E(i,j) + E(i-1,j)) / dy²
				# H_xy = (E(i+1,j+1) - E(i+1,j-1) - E(i-1,j+1) + E(i-1,j-1)) / (4*dx*dy)
				
				center = z_smooth[i, j]
				h_xx = (z_smooth[i, j+1] - 2*center + z_smooth[i, j-1]) / (dx ** 2)
				h_yy = (z_smooth[i+1, j] - 2*center + z_smooth[i-1, j]) / (dy ** 2)
				h_xy = (z_smooth[i+1, j+1] - z_smooth[i+1, j-1] - z_smooth[i-1, j+1] + z_smooth[i-1, j-1]) / (4.0 * dx * dy)
				
				# Create Hessian matrix (now in proper physical units: kJ/Ų)
				hessian = np.array([[h_xx, h_xy], [h_xy, h_yy]])
				
				# Compute eigenvalues (now in kJ/Ų)
				eigenvalues = np.linalg.eigvalsh(hessian)
				
				# Saddle point detection with PHYSICAL thresholds:
				# - Mixed signs (one positive, one negative curvature)
				# - Both eigenvalues SIGNIFICANT in physical units
				#   For a smooth PES, typical saddle points have |λ| > 0.01 kJ/Ų
				#   This eliminates gentle ripples but keeps real transition states
				has_mixed_signs = eigenvalues[0] * eigenvalues[1] < 0
				both_substantial = abs(eigenvalues[0]) > 5.0 and abs(eigenvalues[0]) < 250.0 and abs(eigenvalues[1]) > 10 and abs(eigenvalues[1]) < 250.0 # Physical threshold
				saddle_strength = abs(eigenvalues[0] * eigenvalues[1])  # Product gives saddle strength
				
				if has_mixed_signs and both_substantial:
					saddle_candidates.append({
						'x': j, 'y': i, 'energy': z[i, j], 'eigs': eigenvalues,  # Use original energy, not smoothed
						'strength': saddle_strength
					})
		
		print(f"  Found {len(saddle_candidates)} initial saddle candidates (after physical filtering + smoothing)")
		
		# Second pass: cluster nearby saddle points and keep strongest in each cluster
		if saddle_candidates:
			saddle_points = self._ClusterSaddlePoints(saddle_candidates, cluster_radius=3)
			print(f"  After clustering (radius=3): {len(saddle_points)} saddle points")
			
			# Third pass: energy-based filtering - keep only saddles that are significantly elevated
			min_energy = min(m[2] for m in minima) if minima else min(z.flatten())
			energy_threshold_for_saddle = min_energy + 37.0  # At least 30.0 kJ above minimum
			
			filtered_saddles = [s for s in saddle_points if s['energy'] > energy_threshold_for_saddle]
			print(f"  After energy filtering (E > {energy_threshold_for_saddle:.2f} kJ): {len(filtered_saddles)} saddle points")
			
			# Sort by energy
			filtered_saddles.sort(key=lambda x: x['energy'])
			saddle_points = filtered_saddles
			
			# Store in results
			for s in saddle_points:
				results['saddle_points'].append((s['x'], s['y'], s['energy'], s['eigs'].tolist()))
		else:
			saddle_points = []
			print(f"  No significant saddle points found")
		
		print(f"\n  Final saddle points: {len(saddle_points)}")
		print(f"  (Physical eigenvalue threshold: 10 kJ/Ų)")
		for x_info in saddle_points[:10]:
			print(f"    - Position ({x_info['x']:3d}, {x_info['y']:3d}): Energy = {x_info['energy']:10.4f} kJ, Eigenvalues = [{x_info['eigs'][0]:9.5f}, {x_info['eigs'][1]:9.5f}] kJ/Ų")
		if len(saddle_points) > 10:
			print(f"    ... and {len(saddle_points)-10} more")
		
		print("\n  Saddle points sorted by energy (ascending):")
		
		# 3. Find points with significant gradients (along minimum energy path)
		print("\n[3] Searching for gradient points (transition regions)...")
		gradient_points = []
		
		for i in range(1, self.ylen - 1):
			for j in range(1, self.xlen - 1):
				# Get local spacings for gradient calculation
				
				# Compute gradient magnitude with PHYSICAL spacing (kJ/Å) from smoothed surface
				gx = (z_smooth[i, j+1] - z_smooth[i, j-1]) / (2.0 * dx)
				gy = (z_smooth[i+1, j] - z_smooth[i-1, j]) / (2.0 * dy)
				grad_magnitude = np.sqrt(gx**2 + gy**2)
				
				# High gradient indicates transition region (threshold in physical units: kJ/Å)
				if grad_magnitude > 1.0:  # Physical threshold: significant steepness
					gradient_points.append((j, i, z[i, j], grad_magnitude))
					results['gradient_points'].append((j, i, z[i, j], grad_magnitude))
		
		print(f"  Found {len(gradient_points)} transition region points")
		print(f"  (Physical gradient threshold: 1.0 eV/{chr(197)})")
		# Print highest gradient points
		top_gradient = sorted(gradient_points, key=lambda x: x[3], reverse=True)[:5]
		for x, y, e, grad in top_gradient:
			print(f"    - Position ({x:3d}, {y:3d}): Energy = {e:10.4f} eV, |∇E| = {grad:7.4f} eV/{chr(197)}")
		
		# 4. Identify intermediary structures (between minima with moderate energy)
		print("\n[4] Classifying intermediary structures...")
		if len(minima) >= 2:
			# Get min and max energies
			min_energy = min(m[2] for m in minima)
			max_search_energy = max(m[2] for m in minima) + energy_threshold
			
			intermediates = []
			for i in range(self.ylen):
				for j in range(self.xlen):
					e = z[i, j]
					if min_energy < e < max_search_energy:
						# Check if not already in minima or saddle list
						is_minimum = (j, i, e) in minima
						is_saddle = any((j, i) == (s[0], s[1]) for s in saddle_points)
						if not is_minimum and not is_saddle:
							intermediates.append((j, i, e))
			
			# Deduplicate and sort by energy
			intermediates = list(set(intermediates))
			intermediates.sort(key=lambda x: x[2])
			results['intermediates'] = intermediates
			
			print(f"  Found {len(intermediates)} intermediate points")
			print(f"  Energy range: {min_energy:.4f} to {max_search_energy:.4f} kJ/mol")
			for x, y, e in intermediates[:10]:
				print(f"    - Position ({x:3d}, {y:3d}): Energy = {e:10.4f} kJ/mol")
			if len(intermediates) > 10:
				print(f"    ... and {len(intermediates)-10} more")
		
		# 5. Generate analysis report
		print("\n" + "="*60)
		print("ANALYSIS SUMMARY")
		print("="*60)
		
		report = self._GeneratePESAnalysisReport(results, energy_threshold)
		
		# Save report to file
		report_file = os.path.join(output_folder, "PES_Analysis_Report.txt")
		with open(report_file, 'w') as f:
			f.write(report)
		
		print(f"\nDetailed report saved to: {report_file}")
		
		# 6. Generate visualization with critical points marked
		self._VisualizeCriticalPoints(results, output_folder)
		
		results['analysis_report'] = report
		return results
	
	#----------------------------------------------------------------------------------------
	def _ClusterSaddlePoints(self, candidates, cluster_radius=3):
		"""
		Cluster nearby saddle point candidates and keep the strongest in each cluster.
		
		This prevents detection of multiple nearly-identical saddle points from noise.
		For each cluster, keeps only the point with the highest saddle strength (largest |λ₁*λ₂|).
		
		Args:
			candidates (list): List of saddle point candidate dicts with 'x', 'y', 'strength', etc.
			cluster_radius (int): Distance within which to merge saddle points
		
		Returns:
			list: Clustered saddle points (strongest in each cluster)
		"""
		if not candidates:
			return []
		
		# Sort by saddle strength (highest first)
		sorted_candidates = sorted(candidates, key=lambda c: c['strength'], reverse=True)
		
		clusters = []
		used = set()
		
		for candidate in sorted_candidates:
			if (candidate['x'], candidate['y']) in used:
				continue
			
			# Start new cluster with this candidate
			cluster = [candidate]
			used.add((candidate['x'], candidate['y']))
			
			# Find all nearby unassigned candidates
			for other in sorted_candidates:
				if (other['x'], other['y']) in used:
					continue
				
				# Check distance
				dist = np.sqrt((candidate['x'] - other['x'])**2 + (candidate['y'] - other['y'])**2)
				if dist <= cluster_radius:
					cluster.append(other)
					used.add((other['x'], other['y']))
			
			# Add strongest from this cluster
			strongest = max(cluster, key=lambda c: c['strength'])
			clusters.append(strongest)
		
		return clusters
	
	#----------------------------------------------------------------------------------------
	def _GeneratePESAnalysisReport(self,results, energy_threshold=0.5 ):
		
		"""Generate formatted analysis report."""
		report = ""
		report += "POTENTIAL ENERGY SURFACE (PES) ANALYSIS REPORT\n"
		report += "=" * 60 + "\n\n"
		
		report += "1. LOCAL MINIMA\n"
		report += "-" * 60 + "\n"
		if results['minima']:
			report += f"Total minima found: {len(results['minima'])}\n"
			for x, y, e in results['minima'][:20]:
				report += f"  ({x:4d}, {y:4d}): {e:12.6f} kJ ({e/4.184:10.4f} kcal/mol)\n"
			if len(results['minima']) > 20:
				report += f"  ... and {len(results['minima'])-20} more\n"
		else:
			report += "No local minima found\n"
		report += "\n"
		
		report += "2. SADDLE POINTS (TRANSITION STATES)\n"
		report += "-" * 60 + "\n"
		if results['saddle_points']:
			report += f"Total saddle points found: {len(results['saddle_points'])}\n"
			report += "Position        Energy (eV)      Eigenvalues\n"
			for x, y, e, eigs in results['saddle_points'][:20]:
				report += f"({x:4d}, {y:4d}): {e:12.6f}    [{eigs[0]:8.4f}, {eigs[1]:8.4f}]\n"
			if len(results['saddle_points']) > 20:
				report += f"... and {len(results['saddle_points'])-20} more\n"
		else:
			report += "No saddle points found\n"
		report += "\n"
		
		report += "3. INTERMEDIATE STRUCTURES\n"
		report += "-" * 60 + "\n"
		if results['intermediates']:
			report += f"Total intermediate points found: {len(results['intermediates'])}\n"
			report += f"Energy threshold: {energy_threshold} kcal/mol\n"
			for x, y, e in results['intermediates'][:20]:
				report += f"  ({x:4d}, {y:4d}): {e:12.6f} kJ\n"
			if len(results['intermediates']) > 20:
				report += f"  ... and {len(results['intermediates'])-20} more\n"
		else:
			report += "No intermediate structures found\n"
		report += "\n"
		
		report += "4. TRANSITION REGIONS (HIGH GRADIENT ZONES)\n"
		report += "-" * 60 + "\n"
		if results['gradient_points']:
			report += f"Total high-gradient points: {len(results['gradient_points'])}\n"
			top_gradients = sorted(results['gradient_points'], key=lambda x: x[3], reverse=True)[:10]
			for x, y, e, grad in top_gradients:
				report += f"  ({x:4d}, {y:4d}): |∇E| = {grad:8.4f}, E = {e:10.6f} eV\n"
		report += "\n"
		
		return report
	
	#----------------------------------------------------------------------------------------
	def _VisualizeCriticalPoints(self, results, output_folder):
		"""Create visualization with all critical points marked."""
		try:
			import matplotlib.pyplot as plt
			from matplotlib.patches import Circle
			
			fig, ax = plt.subplots(figsize=(12, 10))
			
			# Plot energy surface
			z = self.energiesMatrix
			X, Y = np.meshgrid(range(self.xlen), range(self.ylen))
			
			# Use log scale for better visualization
			z_plot = z - np.min(z) + 0.1
			levels = np.logspace(np.log10(z_plot.min()), np.log10(z_plot.max()), 15)
			
			contour = ax.contourf(X, Y, z, levels=20, cmap='RdYlBu_r', alpha=0.7)
			ax.contour(X, Y, z, levels=10, colors='black', linewidths=0.5, alpha=0.3)
			
			# Mark minima
			for x, y, e in results['minima']:
				ax.plot(x, y, 'g*', markersize=15, label='Minima' if x == results['minima'][0][0] and y == results['minima'][0][1] else '')
			
			# Mark saddle points
			for x, y, e, eigs in results['saddle_points'][:50]:  # Limit for clarity
				ax.plot(x, y, 'rs', markersize=8, label='Saddle' if x == results['saddle_points'][0][0] and y == results['saddle_points'][0][1] else '')
			
			# Mark intermediates (sample them)
			if results['intermediates']:
				sample_step = max(1, len(results['intermediates']) // 50)
				for x, y, e in results['intermediates'][::sample_step]:
					ax.plot(x, y, 'bo', markersize=3, alpha=0.5)
			
			ax.set_xlabel('Reaction Coordinate 1 (steps)')
			ax.set_ylabel('Reaction Coordinate 2 (steps)')
			ax.set_title('2D PES with Critical Points')
			
			# Add colorbar
			cbar = plt.colorbar(contour, ax=ax)
			cbar.set_label('Energy (kJ)')
			
			# Clean legend
			handles, labels = ax.get_legend_handles_labels()
			by_label = dict(zip(labels, handles))
			ax.legend(by_label.values(), by_label.keys(), loc='upper right')
			
			# Save figure
			fig_path = os.path.join(output_folder, 'PES_Critical_Points.png')
			plt.savefig(fig_path, dpi=150, bbox_inches='tight')
			plt.close()
			
			print(f"Visualization saved to: {fig_path}")
		except Exception as e:
			print(f"Warning: Could not create visualization: {e}")
	
	#----------------------------------------------------------------------------------------
	def FindMinEnergyPath_Dijkstra(self, start_xy, goal_xy, cost_mode='neighbor_energy'):
		"""
		Find minimum cumulative-energy path using Dijkstra's algorithm on the grid.
		
		Args:
			start_xy: (x, y) grid indices of starting point.
			goal_xy:  (x, y) grid indices of goal point.
			cost_mode: 'neighbor_energy' -> cost = energy of destination cell.
					   'edge_average'     -> cost = average of current and destination.
		
		Returns:
			path: list of (x, y) indices from start to goal.
			energies: list of energy values along the path.
			barrier: maximum energy along path minus start energy (kJ/mol).
		"""
		import heapq
		
		print("\n" + "="*70)
		print("FINDING MINIMUM ENERGY PATH (Dijkstra Algorithm)")
		print("="*70)
		print(f"Start: {start_xy}, Goal: {goal_xy}")
		print(f"Start energy: {self.energiesMatrix[start_xy[1]][start_xy[0]]:.4f} kJ/mol")
		print(f"Goal energy:  {self.energiesMatrix[goal_xy[1]][goal_xy[0]]:.4f} kJ/mol")
		print(f"Cost mode: {cost_mode}\n")
		
		# Get grid dimensions and the energy matrix (rows = y, cols = x)
		rows, cols = self.ylen, self.xlen
		z = self.energiesMatrix
		
		# CRITICAL FIX: Offset energies to make all costs positive (Dijkstra requirement)
		# Negative energies cause cumulative costs to become unbounded negative, breaking algorithm
		z_min = np.min(z)
		z_offset = z - z_min + 0.001  # Add small epsilon to ensure strict positivity
		
		print(f"Grid dimensions: {cols} x {rows}")
		print(f"Energy range: {z_min:.6f} to {np.max(z):.6f} kJ/mol")
		print(f"Offset applied: +{(-z_min + 0.001):.6f} kJ/mol (ensures positive step costs)\n")
		
		# dist[y][x] = best cumulative cost found so far to reach (x, y)
		dist = np.full((rows, cols), np.inf)
		prev = np.full((rows, cols), None, dtype=object)
		
		# Start point coordinates
		sx, sy = start_xy
		dist[sy][sx] = 0.0
		heap = [(0.0, sx, sy)]
		
		print(f"Starting Dijkstra search...")
		nodes_explored = 0
		nodes_added = 0
		
		# Main loop: pop smallest cost node, explore neighbours
		while heap:
			current_cost, x, y = heapq.heappop(heap)
			nodes_explored += 1
			
			# Stop if we reached the goal
			if (x, y) == goal_xy:
				print(f"\n✓ Goal reached at ({x}, {y})")
				print(f"  Cumulative cost: {current_cost:.4f}")
				print(f"  Nodes explored: {nodes_explored}")
				break
			
			# If this entry is outdated, skip it
			if current_cost > dist[y][x]:
				continue
			
			# Periodic status output
			if nodes_explored % 100 == 0:
				print(f"  Progress: {nodes_explored} nodes explored, heap size: {len(heap)}, current position: ({x}, {y}), cost: {current_cost:.4f}")
			
			# Explore the 8 neighbours (includes diagonals)
			for dx in (-1, 0, 1):
				for dy in (-1, 0, 1):
					if dx == 0 and dy == 0:
						continue
					nx, ny = x + dx, y + dy
					
					# Check grid boundaries
					if 0 <= nx < cols and 0 <= ny < rows:
						# Compute step cost using OFFSET energies (always positive!)
						if cost_mode == 'neighbor_energy':
							step_cost = z_offset[ny][nx]
						elif cost_mode == 'edge_average':
							step_cost = 0.5 * (z_offset[y][x] + z_offset[ny][nx])
						else:
							step_cost = z_offset[ny][nx]
						
						new_cost = current_cost + step_cost
						
						# If we found a cheaper way to reach (nx, ny), update and push to heap
						if new_cost < dist[ny][nx]:
							dist[ny][nx] = new_cost
							prev[ny][nx] = (x, y)
							heapq.heappush(heap, (new_cost, nx, ny))
							nodes_added += 1
		
		print(f"  Total nodes added to heap: {nodes_added}\n")
		
		# Reconstruct the path from start to goal by walking backwards
		if prev[goal_xy[1]][goal_xy[0]] is None:
			print("✗ Error: Goal not reachable from start.")
			return [], [], np.inf
		
		print(f"Reconstructing path...")
		path = []
		cur = goal_xy
		while cur is not None:
			path.append(cur)
			cur = prev[cur[1]][cur[0]]
		
		path.reverse()
		print(f"✓ Path reconstructed: {len(path)} waypoints\n")
		
		# Extract energies along the path
		energies = [z[y][x] for (x, y) in path]
		self.energies1D = energies
		
		# Calculate barrier height (kJ/mol, convert to kcal/mol for display)
		barrier = max(energies) - energies[0]
		barrier_kcal = barrier / 4.184  # Convert kJ/mol to kcal/mol
		
		print(f"Path Statistics:")
		print(f"  Path length: {len(path)} steps")
		print(f"  Total cumulative cost: {dist[goal_xy[1]][goal_xy[0]]:.4f} kJ/mol")
		print(f"  Start energy: {energies[0]:.6f} kJ/mol")
		print(f"  End energy: {energies[-1]:.6f} kJ/mol")
		print(f"  Max energy on path: {max(energies):.6f} kJ/mol")
		print(f"  Min energy on path: {min(energies):.6f} kJ/mol")
		print(f"  Barrier height: {barrier:.6f} kJ/mol = {barrier_kcal:.4f} kcal/mol")
		print("="*70 + "\n")
		
		return path, energies, barrier_kcal
	
	#----------------------------------------------------------------------------------------
	def FindMinimaxPath(self, start_xy, goal_xy):
		"""
		Path that minimises the maximum energy (lowest barrier).
		Uses a variant of Dijkstra where the cost = max(max_so_far, energy_of_cell).
		This finds the true minimum-energy barrier path (MEP).
		
		Args:
			start_xy: (x, y) grid indices of starting point.
			goal_xy:  (x, y) grid indices of goal point.
		
		Returns:
			path: list of (x, y) indices from start to goal.
			energies: list of energy values along the path.
			barrier: maximum energy along path minus start energy (kcal/mol).
		"""
		import heapq
		
		print("\n" + "="*70)
		print("FINDING MINIMUM BARRIER PATH (Minimax Algorithm)")
		print("="*70)
		print(f"Start: {start_xy}, Goal: {goal_xy}")
		print(f"Start energy: {self.energiesMatrix[start_xy[1]][start_xy[0]]:.4f} kJ/mol")
		print(f"Goal energy:  {self.energiesMatrix[goal_xy[1]][goal_xy[0]]:.4f} kJ/mol")
		print(f"Algorithm: Minimize maximum energy along path (true MEP)\n")
		
		rows, cols = self.ylen, self.xlen
		z = self.energiesMatrix
		
		print(f"Grid dimensions: {cols} x {rows}")
		
		# best_max[y][x] = lowest possible maximum energy to reach (x,y)
		best_max = np.full((rows, cols), np.inf)
		prev = np.full((rows, cols), None, dtype=object)
		
		sx, sy = start_xy
		best_max[sy][sx] = z[sy][sx]        # initial max = start energy
		heap = [(z[sy][sx], sx, sy)]        # (current_max, x, y)
		
		print(f"Starting Minimax search...")
		nodes_explored = 0
		nodes_added = 0
		
		while heap:
			current_max, x, y = heapq.heappop(heap)
			nodes_explored += 1
			
			# Stop if we reached the goal
			if (x, y) == goal_xy:
				print(f"\n✓ Goal reached at ({x}, {y})")
				print(f"  Minimum barrier: {current_max:.4f} kJ/mol")
				print(f"  Nodes explored: {nodes_explored}")
				break
			
			# If this entry is outdated, skip it
			if current_max > best_max[y][x]:
				continue
			
			# Periodic status output
			if nodes_explored % 100 == 0:
				print(f"  Progress: {nodes_explored} nodes explored, heap size: {len(heap)}, current position: ({x}, {y}), max barrier: {current_max:.4f}")
			
			# Explore the 8 neighbours (includes diagonals)
			for dx in (-1, 0, 1):
				for dy in (-1, 0, 1):
					if dx == 0 and dy == 0:
						continue
					nx, ny = x + dx, y + dy
					
					# Check grid boundaries
					if 0 <= nx < cols and 0 <= ny < rows:
						# Minimax: new barrier = max of current barrier and this cell's energy
						new_max = max(current_max, z[ny][nx])
						
						# If we found a path with lower barrier to (nx, ny), update and push to heap
						if new_max < best_max[ny][nx]:
							best_max[ny][nx] = new_max
							prev[ny][nx] = (x, y)
							heapq.heappush(heap, (new_max, nx, ny))
							nodes_added += 1
		
		print(f"  Total nodes added to heap: {nodes_added}\n")
		
		# Reconstruct path from start to goal by walking backwards
		if prev[goal_xy[1]][goal_xy[0]] is None:
			print("✗ Error: Goal not reachable from start.")
			return [], [], np.inf
		
		print(f"Reconstructing path...")
		path = []
		cur = goal_xy
		while cur is not None:
			path.append(cur)
			cur = prev[cur[1]][cur[0]]
		
		path.reverse()
		print(f"✓ Path reconstructed: {len(path)} waypoints\n")
		
		# Extract energies along the path
		energies = [z[y][x] for (x, y) in path]
		self.energies1D = energies
		
		# Calculate barrier height (kJ/mol, convert to kcal/mol for display)
		barrier = best_max[goal_xy[1]][goal_xy[0]] - z[sy][sx]
		barrier_kcal = barrier / 4.184  # Convert kJ/mol to kcal/mol
		
		print(f"Path Statistics:")
		print(f"  Path length: {len(path)} steps")
		print(f"  Start energy: {energies[0]:.6f} kJ/mol")
		print(f"  End energy: {energies[-1]:.6f} kJ/mol")
		print(f"  Max energy on path: {max(energies):.6f} kJ/mol")
		print(f"  Min energy on path: {min(energies):.6f} kJ/mol")
		print(f"  Barrier height: {barrier:.6f} kJ/mol = {barrier_kcal:.4f} kcal/mol")
		print("="*70 + "\n")
		
		return path, energies, barrier_kcal
	#----------------------------------------------------------------------------------------
	def Rewrite_Log(self,_fileName):
		'''
		Rewrite log file with kcat values for each frame
		'''
		min_energy = self.energies1D[0]
		if not self.multiple1Dplot:
			kcats = [0]
			for i in range(1, len(self.energies1D)-1):
				e_prev, e_curr, e_next = self.energies1D[i-1], self.energies1D[i], self.energies1D[i+1]
				if (e_curr < min_energy):
					min_energy = e_curr					
				if (e_curr > e_prev and e_curr > e_next):
					kcats.append( self.Calulate_Kcat( e_curr - min_energy) ) 
				kcats.append(0.0)
			kcats.append(0.0)
			log_text = "x Energy method Energy_kcal kcat \n"
			new_log  = open( _fileName, 'w' )
			for i in range(len(self.energies1D)):										
				energy_kcal = self.energies1D[i] / 4.184
				log_text += "{} {} pickPath {} {}\n".format(i,self.energies1D[i],energy_kcal,kcats[i])
			new_log.write(log_text)
			new_log.close()
		elif self.multiple1Dplot:
			for k in range(self.nplots1D):
				kcats = [0]
				min_energy = self.multiple1Dplot[k][0]
				for i in range(1, len(self.multiple1Dplot[k])-1):
					e_prev, e_curr, e_next = self.multiple1Dplot[k][i-1], self.multiple1Dplot[k][i], self.multiple1Dplot[k][i+1]
					if (e_curr < min_energy):
						min_energy = e_curr	
					if (e_curr > e_prev and e_curr > e_next):
						kcats.append( self.Calulate_Kcat( e_curr - min_energy))
					kcats.append(0.0)
				kcats.append(0.0)
				log_text = "x Energy method Energy_kcal kcat \n"
				new_log  = open( _fileName, 'w' )
				for i in range(len(self.multiple1Dplot[k])):						
					energy_kcal = self.multiple1Dplot[k][i] / 4.184
					log_text += "{} {} pickPath {} {}\n".format(i,self.multiple1Dplot[k][i],energy_kcal,kcats[i])
			new_log.write(log_text)

#================================================================================================#
#======================================END OF THE FILE===========================================#
#================================================================================================#