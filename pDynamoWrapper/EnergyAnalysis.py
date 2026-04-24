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
			self.multiple1Dplot.append(energyTmp)
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
		#----------------------------------
		elif self.Type == "2D":
			for line in reading:
				if i > 0:
					lns = line.split()
					self.RC1.append( float(lns[2] ) )
					self.RC2.append( float(lns[3] ) )
					m = int(lns[0])				
					n = int(lns[1])	
					self.energiesMatrix[n][m] = float(lns[4]) 
				i += 1		
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
		#----------------------------------
		elif self.Type == "FE1D":
			energyTmp = np.zeros( (self.xlen), dtype=float )
			for line in reading:
				lns = line.split()
				m = int( lns[0] )
				energyTmp[m]  = float(lns[1]) 
			self.energies1D = energyTmp		
		#----------------------------------
		elif self.Type == "FE2D":
			for line in reading:
				lns = line.split()
				m = int( lns[0])				
				n = int( lns[1])
				self.energiesMatrix[n][m] = float(lns[2])
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
		plt.savefig(self.baseName+"1d.png",dpi=1000)
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
	#===============================================
	def Plot2D(self,contourlines,crd1label,crd2label,_xlim=None,_ylim=None,SHOW=False,_figS=[7,5],_reverserc1=False,_reverserc2=False):
		'''
		Plot contour plot for potential, free energy and potential of mean field
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
		plt.savefig(plotName+".png",dpi=1000)
		if SHOW: plt.show()
		plt.clf()
		plt.close()
		fig.clf()
	#----------------------------------------------------------------------------------------
	def MultPlot2D(self,contourlines,crd1label,crd2label,_xlim=None,_ylim=None,SHOW=False,_reverserc1=False,_reverserc2=False):
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
		plt.savefig(self.baseName+".png",dpi=1000)
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
	def Path_From_PES(self, in_point,fin_point,_path,_folder_dst,_system):
		'''
		'''
		#setting current point as initial
		cp = in_point

		z = self.energiesMatrix

		pathx = [in_point[0]] 
		pathy = [in_point[1]]
		
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
			try: self.energies1D.append( z[ cp[1], cp[0] ] )
			except: break 
			pathx.append(cp[0])
			pathy.append(cp[1])
		
		#---------------------------------------------------------------
		if not os.path.exists( os.path.join(_folder_dst,"traj1d.ptGeo") ): 
			os.makedirs( os.path.join(_folder_dst,"traj1d.ptGeo") )
		kcats = []
		new_idx = 0
		min_energy = self.energies1D[0]
		for indx in range(len(pathx)):
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
	
		return(pathx, pathy, self.energies1D, kcats)


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


		# 3. If still too few, add evenly spaced points from the original path
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
				add_indices = [candidates[int(i*step)] for i in range(1, needed+1)]
				keep.update(add_indices)
				kept = sorted(keep)
    
    	# 3. If still too many, subsample evenly
		if len(kept) > max_points:
			step = len(kept) / max_points
			kept = [kept[int(i*step)] for i in range(max_points-1)] + [kept[-1]]

		

    
		return kept

	def Analyze_Path(self, path, in_point=None, fin_point=None):
		'''
		Analyze the reaction path and identify potential energy barriers,
		  intermediates, and transition states.
		  Use the identified features to estimate reaction rates, calculate kcat values, and generate a report.
		  calculate kcat values, and export resumed trajectory.
		'''
		if self.Type == "2D" or self.Type == "WHAM2D" or self.Type == "FE2D" or self.Type == "2DRef":
			if in_point == None:
				in_point = [0,0]
			if fin_point == None:
				fin_point = [self.xlen-1,self.ylen-1]
			self.Path_From_PES(in_point,fin_point,path,"path_analysis",_system=None)

		pass

#================================================================================================#
#======================================END OF THE FILE===========================================#
#================================================================================================#



