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
	def Path_From_PES(self, in_point,fin_point,_path,_folder_dst,_system,min_points=15,max_points=21):
		'''
		Calculate minimum energy path from starting point to final point on PES.
		
		Returns:
			tuple: (pathx, pathy) - lists of x and y coordinates of the minimum energy path
		'''
		#setting current point as initial
		cp = in_point

		z = self.energiesMatrix	

		pathx = [in_point[0]] 
		pathy = [in_point[1]]
		self.energies1D.append( z[ pathx[0],pathy[0] ] )
		
		dirs = [ [1,0] ,[0,1], [1,1] ]

		print(cp, fin_point)
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
			else: print( "No Incrementin Y: {}".format(B) )

			if ( cp[0] + 1 ) < self.xlen and ( cp[1] + 1 ) < self.ylen: 
				if  (cp[0] + 1) <=fin_point[0] and (cp[1] + 1) <= fin_point[1]:
					C = z[ cp[1], cp[0] ] + z[ (cp[1] + 1), (cp[0] + 1) ] 
					print( "Increment in both directions:  {}".format(C) )
			else: print("No Incrementin both directions: {}".format(C) )

			D = [ A, B, C ]
			ind = D.index(min(D))			
			cp[0] += dirs[ind][0]
			cp[1] += dirs[ind][1]			
			self.energies1D.append( z[ cp[1], cp[0] ] )
			pathx.append(cp[0])
			pathy.append(cp[1])
		
		
		#---------------------------------------------------------------
		if not os.path.exists( os.path.join(_folder_dst,"traj1d.ptGeo") ): 
			os.makedirs( os.path.join(_folder_dst,"traj1d.ptGeo") )

		self.AnalyzePES2D( output_folder= os.path.join(_folder_dst,"path_analysis") )
		kcats = [0.0]
		new_idx = 0
		min_energy = self.energies1D[0]
		for indx in range( 1,len(pathx) ):
			pkl = _path + "/frame{}_{}.pkl".format(pathx[indx],pathy[indx])			
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
			energy_kcal = self.energies1D[i] * 0.239006
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
			energy_kcal = self.energies1D[i] * 0.239006
			log_text += "{} {} pickPath {} {}\n".format(i,self.energies1D[i],energy_kcal,kcats[i])
		new_log.write(log_text)
		new_log.close()

		#export resampled path frames
		self.energies1D = [self.energies1D[i] for i in Kept_indices]
		self.Type = "1DRef"
		self.plot1d_name = "1D_Ref_resampled.png"
		self.Plot1D("Reaction Path frames (n)")
		for idx in range( len(Kept_indices) ):
			pkl = _path + "/frame{}_{}.pkl".format(pathx[idx],pathy[idx])			
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
	def AnalyzePES2D(self, energy_threshold=0.5, output_folder=None):
		"""
		Analyze 2D potential energy surface to identify critical points.
		
		Detects and classifies critical points on the surface:
		- Local minima: Points with lower energy than all neighbors
		- Saddle points: Points where Hessian has both positive and negative eigenvalues
		- Intermediates: Points between reactant/product regions with moderate energy
		
		Args:
			energy_threshold (float): Energy difference threshold for intermediate detection (kcal/mol)
			output_folder (str): Folder to save analysis results. If None, uses self.baseName
			
		Returns:
			dict: Dictionary containing:
				- 'minima': list of (x, y, energy) tuples for local minima
				- 'saddle_points': list of (x, y, energy, eigenvalues) tuples
				- 'intermediates': list of (x, y, energy) tuples
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
			'inflection_points': []
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
			print(f"    - Position ({x:3d}, {y:3d}): Energy = {e:10.4f} eV")
		if len(minima) > 10:
			print(f"    ... and {len(minima)-10} more")
		
		# 2. Find saddle points using Hessian approximation
		print("\n[2] Searching for saddle points...")
		saddle_points = []
		
		for i in range(1, self.ylen - 1):
			for j in range(1, self.xlen - 1):
				# Approximate Hessian using finite differences
				# H_xx = (E(i,j+1) - 2*E(i,j) + E(i,j-1)) / dx^2
				# H_yy = (E(i+1,j) - 2*E(i,j) + E(i-1,j)) / dy^2
				# H_xy = (E(i+1,j+1) - E(i+1,j-1) - E(i-1,j+1) + E(i-1,j-1)) / 4
				
				center = z[i, j]
				h_xx = z[i, j+1] - 2*center + z[i, j-1]
				h_yy = z[i+1, j] - 2*center + z[i-1, j]
				h_xy = (z[i+1, j+1] - z[i+1, j-1] - z[i-1, j+1] + z[i-1, j-1]) / 4.0
				
				# Create Hessian matrix
				hessian = np.array([[h_xx, h_xy], [h_xy, h_yy]])
				
				# Compute eigenvalues
				eigenvalues = np.linalg.eigvalsh(hessian)
				
				# Saddle point: one positive, one negative eigenvalue
				if eigenvalues[0] * eigenvalues[1] < 0 and abs(eigenvalues[0]) > 0.001 and abs(eigenvalues[1]) > 0.001:
					saddle_points.append((j, i, center, eigenvalues))
					results['saddle_points'].append((j, i, center, eigenvalues.tolist()))
		
		print(f"  Found {len(saddle_points)} saddle points")
		for x, y, e, eigs in saddle_points[:10]:
			print(f"    - Position ({x:3d}, {y:3d}): Energy = {e:10.4f} eV, Eigenvalues = [{eigs[0]:7.3f}, {eigs[1]:7.3f}]")
		if len(saddle_points) > 10:
			print(f"    ... and {len(saddle_points)-10} more")
		
		# 3. Find points with significant gradients (along minimum energy path)
		print("\n[3] Searching for gradient points (transition regions)...")
		gradient_points = []
		
		for i in range(1, self.ylen - 1):
			for j in range(1, self.xlen - 1):
				# Compute gradient magnitude
				gx = (z[i, j+1] - z[i, j-1]) / 2.0
				gy = (z[i+1, j] - z[i-1, j]) / 2.0
				grad_magnitude = np.sqrt(gx**2 + gy**2)
				
				# High gradient indicates transition region
				if grad_magnitude > 0.5:  # Threshold for significant gradient
					gradient_points.append((j, i, z[i, j], grad_magnitude))
					results['gradient_points'].append((j, i, z[i, j], grad_magnitude))
		
		print(f"  Found {len(gradient_points)} transition region points")
		# Print highest gradient points
		top_gradient = sorted(gradient_points, key=lambda x: x[3], reverse=True)[:5]
		for x, y, e, grad in top_gradient:
			print(f"    - Position ({x:3d}, {y:3d}): Energy = {e:10.4f} eV, |∇E| = {grad:7.4f}")
		
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
			print(f"  Energy range: {min_energy:.4f} to {max_search_energy:.4f} eV")
			for x, y, e in intermediates[:10]:
				print(f"    - Position ({x:3d}, {y:3d}): Energy = {e:10.4f} eV")
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
	def _GeneratePESAnalysisReport(self, results, energy_threshold):
		"""Generate formatted analysis report."""
		report = ""
		report += "POTENTIAL ENERGY SURFACE (PES) ANALYSIS REPORT\n"
		report += "=" * 60 + "\n\n"
		
		report += "1. LOCAL MINIMA\n"
		report += "-" * 60 + "\n"
		if results['minima']:
			report += f"Total minima found: {len(results['minima'])}\n"
			for x, y, e in results['minima'][:20]:
				report += f"  ({x:4d}, {y:4d}): {e:12.6f} eV ({e*0.239006:10.4f} kcal/mol)\n"
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
				report += f"  ({x:4d}, {y:4d}): {e:12.6f} eV\n"
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
			cbar.set_label('Energy (eV)')
			
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
	def GenerateNEBPathsFromAnalysis(self, analysis_results, in_point, fin_point, path_folder, 
	                                   output_folder, system, max_distance=10, interpolation_points=None):
		"""
		Generate multiple NEB starting paths using analyzed critical points.
		
		Creates trajectories connecting critical points identified from PES analysis:
		- Minima → Saddle → Minima (traditional TS pathways)
		- Minima → Intermediates → Minima (through intermediate structures)
		- Minima → Gradient change points → Minima (through high-gradient regions)
		
		Each path is exported as a ptGeo folder with renumbered .pkl files for NEB calculations.
		This provides diverse initial guesses by exploring different reaction mechanisms.
		
		Args:
			analysis_results (dict): Results from AnalyzePES2D with minima, saddle_points, intermediates, gradient_points
			in_point (tuple): (x, y) initial point coordinates
			fin_point (tuple): (x, y) final point coordinates
			path_folder (str): Folder containing the frame .pkl files
			output_folder (str): Destination folder for NEB trajectory folders
			system: Molecular system object for coordinate operations
			max_distance (int): Maximum distance from in/fin point to consider as "close" (pixels)
			interpolation_points (int): Number of interpolated points between critical points. If None, uses existing frames
		
		Returns:
			dict: Information about generated paths with keys:
				- 'paths': list of generated path info dictionaries
				- 'summary_report': formatted text report
				- 'path_types': breakdown of path types generated
		"""
		print("\n" + "="*70)
		print("GENERATING NEB PATHS FROM CRITICAL POINT ANALYSIS")
		print("="*70)
		
		if 'minima' not in analysis_results or 'saddle_points' not in analysis_results:
			print("Error: Analysis results must contain 'minima' and 'saddle_points'")
			return None
		
		minima = analysis_results['minima']
		saddles = analysis_results['saddle_points']
		intermediates = analysis_results.get('intermediates', [])
		gradient_points = analysis_results.get('gradient_points', [])
		
		if not os.path.exists(output_folder):
			os.makedirs(output_folder)
		
		# 1. Find minima close to initial point
		print("\n[1] Finding minima close to initial point...")
		initial_minima = []
		for x, y, e in minima:
			dist = np.sqrt((x - in_point[0])**2 + (y - in_point[1])**2)
			if dist <= max_distance:
				initial_minima.append((x, y, e, dist))
		
		initial_minima.sort(key=lambda p: p[3])  # Sort by distance
		print(f"    Found {len(initial_minima)} minima near initial point")
		for x, y, e, d in initial_minima[:5]:
			print(f"      ({x:3d}, {y:3d}): Energy = {e:10.6f} eV, Distance = {d:6.2f}")
		
		# 2. Find minima close to final point
		print("\n[2] Finding minima close to final point...")
		final_minima = []
		for x, y, e in minima:
			dist = np.sqrt((x - fin_point[0])**2 + (y - fin_point[1])**2)
			if dist <= max_distance:
				final_minima.append((x, y, e, dist))
		
		final_minima.sort(key=lambda p: p[3])  # Sort by distance
		print(f"    Found {len(final_minima)} minima near final point")
		for x, y, e, d in final_minima[:5]:
			print(f"      ({x:3d}, {y:3d}): Energy = {e:10.6f} eV, Distance = {d:6.2f}")
		
		# 3. Find saddle points between regions
		print("\n[3] Identifying transition states (saddle points)...")
		relevant_saddles = []
		for x, y, e, eigs in saddles:
			relevant_saddles.append((x, y, e, eigs))
		
		print(f"    Found {len(relevant_saddles)} saddle points")
		for x, y, e, eigs in relevant_saddles[:5]:
			print(f"      ({x:3d}, {y:3d}): Energy = {e:10.6f} eV, λ = [{eigs[0]:7.3f}, {eigs[1]:7.3f}]")
		
		# 4. Sort and filter intermediates and gradient points
		print("\n[4] Processing intermediary structures...")
		relevant_intermediates = sorted(intermediates, key=lambda p: p[2])  # Sort by energy
		print(f"    Found {len(relevant_intermediates)} intermediate points")
		for x, y, e in relevant_intermediates[:5]:
			print(f"      ({x:3d}, {y:3d}): Energy = {e:10.6f} eV")
		
		print("\n[5] Processing high-gradient transition regions...")
		# Sort gradient points by magnitude (highest first)
		relevant_gradients = sorted(gradient_points, key=lambda p: p[3], reverse=True)
		print(f"    Found {len(relevant_gradients)} gradient change points")
		for x, y, e, grad in relevant_gradients[:5]:
			print(f"      ({x:3d}, {y:3d}): |∇E| = {grad:7.4f}, Energy = {e:10.6f} eV")
		
		# 5. Generate paths
		print("\n[6] Generating paths connecting critical points...")
		generated_paths = []
		path_types = {'min_saddle_min': 0, 'min_intermediate_min': 0, 'min_gradient_min': 0, 
		              'min_intermediate_saddle_min': 0, 'min_gradient_saddle_min': 0}
		path_id = 0
		
		# Type 1: Traditional min → saddle → min paths
		print("    [Type 1] Min → Saddle → Min paths...")
		for init_x, init_y, init_e, init_d in initial_minima[:3]:  # Limit to 3 closest
			for final_x, final_y, final_e, final_d in final_minima[:3]:  # Limit to 3 closest
				for sad_x, sad_y, sad_e, sad_eigs in relevant_saddles[:5]:  # Check multiple saddles
					
					# Create path: initial minimum → saddle → final minimum
					path_coords_x = [init_x, sad_x, final_x]
					path_coords_y = [init_y, sad_y, final_y]
					path_energies = [init_e, sad_e, final_e]
					
					barrier_height = (sad_e - min(init_e, final_e)) * 0.239006  # Convert to kcal/mol
					
					path_info = {
						'id': path_id,
						'type': 'min_saddle_min',
						'initial': (init_x, init_y, init_e),
						'middle': (sad_x, sad_y, sad_e),
						'final': (final_x, final_y, final_e),
						'barrier': barrier_height,
						'coords_x': path_coords_x,
						'coords_y': path_coords_y,
						'energies': path_energies
					}
					
					generated_paths.append(path_info)
					path_types['min_saddle_min'] += 1
					path_id += 1
		
		# Type 2: min → intermediate → min paths
		if relevant_intermediates:
			print("    [Type 2] Min → Intermediate → Min paths...")
			for init_x, init_y, init_e, init_d in initial_minima[:3]:
				for final_x, final_y, final_e, final_d in final_minima[:3]:
					for int_x, int_y, int_e in relevant_intermediates[:5]:  # Use top 5 intermediates
						
						path_coords_x = [init_x, int_x, final_x]
						path_coords_y = [init_y, int_y, final_y]
						path_energies = [init_e, int_e, final_e]
						
						barrier_height = (int_e - min(init_e, final_e)) * 0.239006
						
						path_info = {
							'id': path_id,
							'type': 'min_intermediate_min',
							'initial': (init_x, init_y, init_e),
							'middle': (int_x, int_y, int_e),
							'final': (final_x, final_y, final_e),
							'barrier': barrier_height,
							'coords_x': path_coords_x,
							'coords_y': path_coords_y,
							'energies': path_energies
						}
						
						generated_paths.append(path_info)
						path_types['min_intermediate_min'] += 1
						path_id += 1
		
		# Type 3: min → gradient point → min paths
		if relevant_gradients:
			print("    [Type 3] Min → Gradient point → Min paths...")
			for init_x, init_y, init_e, init_d in initial_minima[:3]:
				for final_x, final_y, final_e, final_d in final_minima[:3]:
					for grad_x, grad_y, grad_e, grad_mag in relevant_gradients[:5]:  # Use top 5 gradient points
						
						path_coords_x = [init_x, grad_x, final_x]
						path_coords_y = [init_y, grad_y, final_y]
						path_energies = [init_e, grad_e, final_e]
						
						barrier_height = (grad_e - min(init_e, final_e)) * 0.239006
						
						path_info = {
							'id': path_id,
							'type': 'min_gradient_min',
							'initial': (init_x, init_y, init_e),
							'middle': (grad_x, grad_y, grad_e),
							'final': (final_x, final_y, final_e),
							'barrier': barrier_height,
							'gradient_magnitude': grad_mag,
							'coords_x': path_coords_x,
							'coords_y': path_coords_y,
							'energies': path_energies
						}
						
						generated_paths.append(path_info)
						path_types['min_gradient_min'] += 1
						path_id += 1
		
		# Type 4: min → intermediate → saddle → min (for more complex paths)
		if relevant_intermediates and relevant_saddles:
			print("    [Type 4] Min → Intermediate → Saddle → Min paths...")
			for init_x, init_y, init_e, init_d in initial_minima[:2]:
				for final_x, final_y, final_e, final_d in final_minima[:2]:
					for int_x, int_y, int_e in relevant_intermediates[:3]:
						for sad_x, sad_y, sad_e, sad_eigs in relevant_saddles[:3]:
							
							path_coords_x = [init_x, int_x, sad_x, final_x]
							path_coords_y = [init_y, int_y, sad_y, final_y]
							path_energies = [init_e, int_e, sad_e, final_e]
							
							barrier_height = (sad_e - min(init_e, final_e)) * 0.239006
							
							path_info = {
								'id': path_id,
								'type': 'min_intermediate_saddle_min',
								'initial': (init_x, init_y, init_e),
								'intermediate': (int_x, int_y, int_e),
								'saddle': (sad_x, sad_y, sad_e),
								'final': (final_x, final_y, final_e),
								'barrier': barrier_height,
								'coords_x': path_coords_x,
								'coords_y': path_coords_y,
								'energies': path_energies
							}
							
							generated_paths.append(path_info)
							path_types['min_intermediate_saddle_min'] += 1
							path_id += 1
		
		# Type 5: min → gradient → saddle → min
		if relevant_gradients and relevant_saddles:
			print("    [Type 5] Min → Gradient → Saddle → Min paths...")
			for init_x, init_y, init_e, init_d in initial_minima[:2]:
				for final_x, final_y, final_e, final_d in final_minima[:2]:
					for grad_x, grad_y, grad_e, grad_mag in relevant_gradients[:3]:
						for sad_x, sad_y, sad_e, sad_eigs in relevant_saddles[:3]:
							
							path_coords_x = [init_x, grad_x, sad_x, final_x]
							path_coords_y = [init_y, grad_y, sad_y, final_y]
							path_energies = [init_e, grad_e, sad_e, final_e]
							
							barrier_height = (sad_e - min(init_e, final_e)) * 0.239006
							
							path_info = {
								'id': path_id,
								'type': 'min_gradient_saddle_min',
								'initial': (init_x, init_y, init_e),
								'gradient': (grad_x, grad_y, grad_e),
								'saddle': (sad_x, sad_y, sad_e),
								'final': (final_x, final_y, final_e),
								'barrier': barrier_height,
								'gradient_magnitude': grad_mag,
								'coords_x': path_coords_x,
								'coords_y': path_coords_y,
								'energies': path_energies
							}
							
							generated_paths.append(path_info)
							path_types['min_gradient_saddle_min'] += 1
							path_id += 1
		
		print(f"    Generated {len(generated_paths)} total candidate paths")
		print(f"      - Type 1 (Min→Saddle→Min): {path_types['min_saddle_min']}")
		print(f"      - Type 2 (Min→Intermediate→Min): {path_types['min_intermediate_min']}")
		print(f"      - Type 3 (Min→Gradient→Min): {path_types['min_gradient_min']}")
		print(f"      - Type 4 (Min→Intermediate→Saddle→Min): {path_types['min_intermediate_saddle_min']}")
		print(f"      - Type 5 (Min→Gradient→Saddle→Min): {path_types['min_gradient_saddle_min']}")
		
		# 5. Extract frames along each path and export
		print("\n[5] Exporting paths as NEB starting trajectories...")
		exported_paths = []
		
		for path_info in generated_paths:
			path_id = path_info['id']
			output_path_folder = os.path.join(output_folder, f"NEB_Path_{path_id:03d}.ptGeo")
			
			if not os.path.exists(output_path_folder):
				os.makedirs(output_path_folder)
			
			# Collect frames along the path
			frames_data = []
			frame_count = 0
			
			# Process initial minimum
			init_file = os.path.join(path_folder, f"frame{path_info['coords_x'][0]}_{path_info['coords_y'][0]}.pkl")
			if os.path.exists(init_file):
				try:
					system.coordinates3 = ImportCoordinates3(init_file, log=None)
					output_file = os.path.join(output_path_folder, f"frame{frame_count:04d}.pkl")
					Pickle(output_file, system.coordinates3)
					frames_data.append({
						'frame': frame_count,
						'x': path_info['coords_x'][0],
						'y': path_info['coords_y'][0],
						'energy': path_info['energies'][0]
					})
					frame_count += 1
				except Exception as e:
					print(f"      Warning: Could not process initial point: {e}")
			
			# Process saddle point
			sad_file = os.path.join(path_folder, f"frame{path_info['coords_x'][1]}_{path_info['coords_y'][1]}.pkl")
			if os.path.exists(sad_file):
				try:
					system.coordinates3 = ImportCoordinates3(sad_file, log=None)
					output_file = os.path.join(output_path_folder, f"frame{frame_count:04d}.pkl")
					Pickle(output_file, system.coordinates3)
					frames_data.append({
						'frame': frame_count,
						'x': path_info['coords_x'][1],
						'y': path_info['coords_y'][1],
						'energy': path_info['energies'][1]
					})
					frame_count += 1
				except Exception as e:
					print(f"      Warning: Could not process saddle point: {e}")
			
			# Process final minimum
			fin_file = os.path.join(path_folder, f"frame{path_info['coords_x'][2]}_{path_info['coords_y'][2]}.pkl")
			if os.path.exists(fin_file):
				try:
					system.coordinates3 = ImportCoordinates3(fin_file, log=None)
					output_file = os.path.join(output_path_folder, f"frame{frame_count:04d}.pkl")
					Pickle(output_file, system.coordinates3)
					frames_data.append({
						'frame': frame_count,
						'x': path_info['coords_x'][2],
						'y': path_info['coords_y'][2],
						'energy': path_info['energies'][2]
					})
					frame_count += 1
				except Exception as e:
					print(f"      Warning: Could not process final point: {e}")
			
			# Write log file for this path
			log_file = os.path.join(output_path_folder, "NEB_Path.log")
			with open(log_file, 'w') as f:
				f.write("Frame  X_Index  Y_Index  Energy(eV)  Barrier(kcal/mol)\n")
				ref_energy = min(path_info['energies'][0], path_info['energies'][2])
				for data in frames_data:
					barrier = (data['energy'] - ref_energy) * 0.239006
					f.write(f"{data['frame']:5d}  {data['x']:7d}  {data['y']:7d}  {data['energy']:10.6f}  {barrier:10.6f}\n")
			
			# Write path info file
			info_file = os.path.join(output_path_folder, "Path_Info.txt")
			with open(info_file, 'w') as f:
				f.write("NEB PATH INFORMATION\n")
				f.write("="*60 + "\n\n")
				f.write(f"Path ID: {path_id}\n")
				f.write(f"Total frames: {frame_count}\n")
				f.write(f"Barrier height: {path_info['barrier']:.4f} kcal/mol\n\n")
				f.write("Path segments:\n")
				f.write(f"  Initial minimum:  ({path_info['initial'][0]}, {path_info['initial'][1]}), E = {path_info['initial'][2]:.6f} eV\n")
				f.write(f"  Saddle point:     ({path_info['saddle'][0]}, {path_info['saddle'][1]}), E = {path_info['saddle'][2]:.6f} eV\n")
				f.write(f"  Final minimum:    ({path_info['final'][0]}, {path_info['final'][1]}), E = {path_info['final'][2]:.6f} eV\n")
				f.write(f"\nRenumbered frames (0 to {frame_count-1}) ready for NEB calculation\n")
			
			if frame_count >= 3:  # Only count as valid if we got all 3 key points
				exported_paths.append({
					'id': path_id,
					'folder': output_path_folder,
					'frames': frame_count,
					'barrier': path_info['barrier']
				})
				print(f"    ✓ Path {path_id:03d}: {frame_count} frames, Barrier = {path_info['barrier']:.4f} kcal/mol")
			else:
				print(f"    ✗ Path {path_id:03d}: Insufficient frames ({frame_count})")
		
		# 6. Generate summary report
		print("\n[6] Generating summary report...")
		summary_report = self._GenerateNEBPathReport(exported_paths, analysis_results)
		
		summary_file = os.path.join(output_folder, "NEB_Paths_Summary.txt")
		with open(summary_file, 'w') as f:
			f.write(summary_report)
		
		print(f"    Summary saved to: {summary_file}")
		
		print("\n" + "="*70)
		print(f"SUCCESS: Generated {len(exported_paths)} NEB starting paths")
		print(f"Output folder: {output_folder}")
		print("="*70 + "\n")
		
		return {
			'paths': exported_paths,
			'summary_report': summary_report,
			'total_generated': len(exported_paths)
		}
	
	#----------------------------------------------------------------------------------------
	def _GenerateNEBPathReport(self, exported_paths, analysis_results):
		"""Generate formatted report of generated NEB paths."""
		report = ""
		report += "NEB STARTING PATHS - GENERATION REPORT\n"
		report += "="*70 + "\n\n"
		
		report += "SURFACE STATISTICS\n"
		report += "-"*70 + "\n"
		report += f"Total local minima found: {len(analysis_results['minima'])}\n"
		report += f"Total saddle points found: {len(analysis_results['saddle_points'])}\n"
		report += f"Total NEB paths generated: {len(exported_paths)}\n\n"
		
		report += "GENERATED NEB PATHS (sorted by barrier height)\n"
		report += "-"*70 + "\n"
		
		sorted_paths = sorted(exported_paths, key=lambda x: x['barrier'])
		
		for i, path in enumerate(sorted_paths, 1):
			report += f"\n{i}. Path ID: {path['id']:03d}\n"
			report += f"   Frames: {path['frames']}\n"
			report += f"   Barrier: {path['barrier']:.4f} kcal/mol\n"
			report += f"   Folder: NEB_Path_{path['id']:03d}.ptGeo\n"
			report += f"   Status: Ready for NEB calculation\n"
		
		report += "\n" + "-"*70 + "\n"
		report += "USAGE INSTRUCTIONS\n"
		report += "-"*70 + "\n"
		report += "Each path folder contains:\n"
		report += "  • frame0000.pkl, frame0001.pkl, ... - Renumbered coordinate files\n"
		report += "  • NEB_Path.log - Energy profile along path\n"
		report += "  • Path_Info.txt - Path details and statistics\n\n"
		report += "To use with NEB calculations:\n"
		report += "  1. Load the frame files as your initial guess\n"
		report += "  2. Renumbered frames ensure correct ordering (0 to n)\n"
		report += "  3. Use the barrier height to estimate convergence criteria\n"
		report += "  4. Check lowest-barrier paths first for efficiency\n"
		
		report += "\n" + "="*70 + "\n"
		return report
	
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
				energy_kcal = self.energies1D[i] * 0.239006
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
					energy_kcal = self.multiple1Dplot[k][i] * 0.239006
					log_text += "{} {} pickPath {} {}\n".format(i,self.multiple1Dplot[k][i],energy_kcal,kcats[i])
			new_log.write(log_text)

#================================================================================================#
#======================================END OF THE FILE===========================================#
#================================================================================================#



