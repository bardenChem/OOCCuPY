#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = Analysis.py



#=============================================================
import os, sys, glob,shutil
import numpy as np
#--------------------------------------------------------------
import matplotlib.pyplot as plt
#--------------------------------------------------------------
from .commonFunctions 			import *
from pBabel                     import *                                     
from pCore                      import *                                     
from pMolecule                  import *                   
from pScientific                import * 
from pScientific.Symmetry       import *                                        
from pScientific.Statistics     import *
from pScientific.Arrays         import *
from pSimulation                import *


#=====================================================================
class TrajectoryAnalysis:
	'''
	Functions to perform analysis from molecular dynamics trajectories
	'''
	#-----------------------------------------
	def __init__(self,_trajFolder,_system,t_time):
		'''
		Default constructor. Initializa the atributes. 
		'''
		self.trajFolder = _trajFolder
		self.molecule   = _system
		self.RG         = []
		self.RMS        = []
		self.total_time = t_time
		self.rc1_MF     = 0.0
		self.rc2_MF     = 0.0
		self.rg_MF      = 0.0
		self.rms_MF     = 0.0
		self.energies   = []
		self.label_RC1  = None 
		self.label_RC2  = None 
		self.RCs        = {}

		self.trajectory = None

		if self.trajFolder[-4:] == ".dcd":			
			Duplicate(self.trajFolder,self.trajFolder[:-4]+".ptGeo",self.molecule)
			self.trajFolder = self.trajFolder[:-4]+".ptGeo"
			self.trajectory = ImportTrajectory( self.trajFolder , self.molecule )
		elif self.trajFolder[-6:] == ".ptGeo":
			self.trajectory = ImportTrajectory( self.trajFolder , self.molecule )
		
		self.trajectory.ReadHeader()		

	#=================================================
	def Split_Traj(self,_bp):
		'''
		'''
		first_part_trj  = self.trajFolder[:-4] + "_1st_trj.ptGeo"
		second_part_trj = self.trajFolder[:-4] + "_2nd_trj.ptGeo"

		print(first_part_trj,second_part_trj)

		if not os.path.exists(first_part_trj):  os.makedirs(first_part_trj)
		if not os.path.exists(second_part_trj): os.makedirs(second_part_trj)

		pkls = glob.glob(self.trajFolder+"/*.pkl")

		cnt = 0 
		for pkl in pkls:
			newName = "framenone.pkl"
			if cnt < _bp:
				newName = os.path.join(first_part_trj,"frame{}.pkl".format( cnt ) )
			elif cnt >= _bp:
				newName = os.path.join(second_part_trj,"frame{}.pkl".format( cnt-_bp ) )
			cnt +=1
			shutil.copy(pkl,newName)

    #=================================================
	def CalculateRG_RMSD(self,qc_mm=False,protein=False):
		'''
		Get Radius of Gyration and Root Mean Square distance for the trajectory
		'''
		masses    = Array.FromIterable ( [ atom.mass for atom in self.molecule.atoms ] )
		self.crd3 = Unpickle(os.path.join(self.trajFolder,"frame0.pkl"))[0]
		system  = None 
		rg0     = None
		try:
			self.molecule.coordinates3 = self.crd3
		except:
			self.crd3 = Unpickle(os.path.join(self.trajFolder,"frame0.pkl"))
			self.molecule.coordinates3 = self.crd3
		energy0 = self.molecule.Energy(log=None)
		if qc_mm:
			system= Selection( list(self.molecule.qcState.pureQCAtoms) )
			# . Calculate the radius of gyration.
			rg0 = self.crd3.RadiusOfGyration(selection = system, weights = masses)
		elif protein:
			system  = AtomSelection.FromAtomPattern ( self.molecule, "*:*:CA" )		
			# . Calculate the radius of gyration.
			rg0 = self.crd3.RadiusOfGyration(selection = system, weights = masses)
		else:
			system  = Selection.FromIterable ( range ( len ( self.molecule.atoms ) ) )
			# . Calculate the radius of gyration.
			rg0 = self.crd3.RadiusOfGyration(selection = system, weights = masses)
		#------------------------------------------------------------------------------
		# . Save the starting coordinates.
		reference3 = Clone(self.crd3)  
		#------------------------------------------------------------------------------
		n = []
		m = 0   
		#-------------------------------------------------------------------------------
		while self.trajectory.RestoreOwnerData():
			self.energies.append( self.molecule.Energy(log=None) - energy0 )
			self.molecule.coordinates3.Superimpose ( reference3, selection = system, weights = masses )
			self.RG.append  ( self.molecule.coordinates3.RadiusOfGyration( selection = system, weights = masses ) )
			self.RMS.append ( self.molecule.coordinates3.RootMeanSquareDeviation( reference3, selection = system, weights = masses ) )
			n.append(m)
			m+=1
		# . Set up the statistics calculations.        
		rgStatistics  = Statistics(self.RG)
		rmsStatistics = Statistics(self.RMS)
		#-------------------------------------------------------------------------------
		# . Save the results.        
		textLog = open( self.trajFolder+"_MDanalysis", "w" ) 
		#-------------------------------------------------------------------------------
		_Text = "rg0 rgMean rgSD rgMax rgMin\n"
		_Text += "{} {} {} {} {}\n".format(rg0,rgStatistics.mean,rgStatistics.standardDeviation,rgStatistics.maximum,rgStatistics.minimum )
		_Text += "rmsMean rmsSD rmsMax rmsMin\n"
		_Text += "{} {} {} {}\n".format(rmsStatistics.mean,rmsStatistics.standardDeviation,rmsStatistics.maximum,rmsStatistics.minimum )
		#-------------------------------------------------------------------------------
		_Text += "Frame RG RMS\n"
		for i in range(len(self.RG)):
			_Text += "{} {} {}\n".format(i,self.RG[i],self.RMS[i])
		#--------------------------------------------------------------------------------
		textLog.write(_Text)
		textLog.close()
	#===================================================================================================
	def Calculate_RDF(self,_selection_1,_selection_2=None,_selection_name="",_box_size= 25.0):
		'''
		'''

		self.molecule.symmetry =  PeriodicBoundaryConditions.WithCrystalSystem ( CrystalSystemCubic ( ) )
		self.molecule.symmetryParameters = self.molecule.symmetry.MakeSymmetryParameters ( a = _box_size)
		print(self.molecule.symmetryParameters)
		rdf_dat = RadialDistributionFunction ( self.trajFolder, self.molecule, selection1 = _selection_1, selection2=_selection_2, upper=10.0 )

		textLog = open( self.trajFolder+"_"+_selection_name+"_"+"_rdf.log", "w" )
		_text = "Distance(A) G(r) \n"

		for i in range( len(rdf_dat[0]) ):	_text += "{} {}\n".format(rdf_dat[0][i],rdf_dat[1][i])
		textLog.write(_text)
		textLog.close()

		fig1, (ax1) = plt.subplots(nrows=1)
		plt.plot(rdf_dat[0], rdf_dat[1])
		ax1.set_xlabel("Distance $\AA$")
		ax1.set_ylabel("G(r)")
		plt.savefig(self.trajFolder+"_"+_selection_name+"_rdf.png")
		plt.close()
		fig1.clf()

	#-----------------------------------------------------------------------------------------------------
	def Calculate_SD(self,_selection,_selection_name=""):
		'''
		'''		
		rdf_dat = SelfDiffusionFunction( self.trajFolder, self.molecule, selection  = _selection )
		
		textLog = open( self.trajFolder+"_"+_selection_name+"_"+"_sdf.log", "w" )
		_text = "Time(ps) SDF \n"

		for i in range( len(rdf_dat[0]) ):	_text += "{} {}\n".format(rdf_dat[0][i],rdf_dat[1][i])
		textLog.write(_text)
		textLog.close()

		fig1, (ax1) = plt.subplots(nrows=1)
		plt.plot(rdf_dat[0], rdf_dat[1])
		ax1.set_xlabel("Time (ps)")
		ax1.set_ylabel("Dself")
		plt.savefig(self.trajFolder+"_"+_selection_name+"_sdf.png")
		plt.close()
		fig1.clf()
	#--------------------------------------------------
	def ExtractFrames(self):
		'''			
        '''
		try: 	from sklearn.neighbors import KernelDensity
		except:	pass
		kde = KernelDensity(bandwidth=1.0, kernel='gaussian')
		
		self.RMS        = np.array(self.RMS, dtype=np.float32)
		self.RG         = np.array(self.RG, dtype=np.float32)
		self.RMS.reshape(-1,1)
		self.RG.reshape(-1,1)
		
		try:
			kde.fit(self.RMS[:,None])
			density_rms = np.exp(kde.score_samples(self.RMS[:,None]))
			self.rms_MF = max(density_rms[:,None])
			kde.fit(self.RG[:,None])
			density_rg  = np.exp(kde.score_samples(self.RG[:,None]))
			self.rg_MF  = max(density_rg[:,None])
		
		#------------------------------------------------------------------------------
			distold = abs(density_rms[0] - self.rms_MF)
			distnew = 0.0
			fn      = 0 
			for i in range( len(density_rms) ):
				distnew = abs(density_rms[i] - self.rms_MF)
				if distnew < distold:
					distold = distnew
					fn = i
		#------------------------------------------------
			print( os.path.join(self.trajFolder,"frame{}.pkl".format(fn) ) )
			self.molecule.coordinates3 = ImportSystem( os.path.join(self.trajFolder,"frame{}.pkl".format(fn) ) )
			ExportSystem( os.path.join( self.trajFolder, "mostFrequentRMS.pdb" ),self.molecule,log=None )
			ExportSystem( os.path.join( self.trajFolder, "mostFrequentRMS.pkl" ),self.molecule,log=None )
		#------------------------------------------------------------------------------
		except:	pass	
		self.molecule.coordinates3 = AveragePositions(self.trajFolder,self.molecule)
		ExportSystem( os.path.join( self.trajFolder,"Average.pdb"), self.molecule,log=None  )
		ExportSystem( os.path.join( self.trajFolder,"Average.pkl"), self.molecule,log=None )
		self.molecule.coordinates3 = self.crd3

	#=================================================
	def ExtractFrames_biplot(self,rc_1,rc_2):
		'''
		'''
		try: 	from sklearn.neighbors import KernelDensity
		except:	
			print("Error loading KernelDensity library!")
			exit()
		kde = KernelDensity(bandwidth=1.0, kernel='gaussian')

		distances1 = np.array(rc_1[1], dtype=np.float32)
		distances2 = np.array(rc_2[1], dtype=np.float32)

		kde.fit(distances1[:, None])
		density_rc1 = kde.score_samples(distances1[:,None])
		density_rc1 = np.exp(density_rc1)
		rc1_MF = max(density_rc1)
		distances2.reshape(-1,1)
		kde.fit(distances2[:,None])
		density_rc2 = np.exp(kde.score_samples(distances2[:,None]))
		rc2_MF = max(density_rc2)

		distoldRC1 = abs(density_rc1[0] - rc1_MF)
		distoldRC2 = abs(density_rc2[0] - rc2_MF)
		distold    = abs(distoldRC1-distoldRC2)
		distnew    = 0.0
		fn         = 0 
		for i in range( len(distances1) ):
			distnew = abs( abs(density_rc1[i] - rc1_MF) -  abs(density_rc2[i] - rc2_MF) )
			if distnew < distold:
				distold = distnew
				fn = i		
		a  = Unpickle( os.path.join(self.trajFolder,"frame{}.pkl".format(fn) ) )
		b  = Unpickle( os.path.join(self.trajFolder,"frame{}.pkl".format( len(distances2)-1 ) ) )

		self.molecule.coordinates3 = a[0]
		ExportSystem( os.path.join( self.trajFolder,"mostFrequentRC1RC2.pdb"), self.molecule,log=None  )
		ExportSystem( os.path.join( self.trajFolder,"mostFrequentRC1RC2.pkl"), self.molecule,log=None )

		try:
			import seaborn as sns
			g=sns.jointplot(x=distances1,y=distances2,kind="kde",cmap="plasma",shade=True,height=6,widht=8)
			g.set_axis_labels(rc_1[0].label,rc_2[0].label)
			plt.savefig( os.path.join( self.trajFolder,label_text+"_Biplot.png"),dpi=1000 )
			if SHOW: plt.show()
			plt.close()
		except:
			print("Error in importing seaborn package!\nSkipping biplot distribution plot!")
			pass
		
	#=================================================
	def PlotRG_RMS(self,SHOW=False):
		'''
		Plot graphs for the variation of Radius of Gyration and RMSD
		'''
		fig1, (ax1) = plt.subplots(nrows=1)

		n = np.linspace( 0, self.total_time, len(self.RG) )
		plt.plot(n, self.RG)
		ax1.set_xlabel("Time (ps)")
		ax1.set_ylabel("Radius of Gyration $\AA$")
		plt.savefig( os.path.join( self.trajFolder,"analysis_mdRG.png") )
		if SHOW: plt.show()
		plt.clf()
		plt.close()
		fig1.clf()

		fig2, (ax2) = plt.subplots(nrows=1)
		#--------------------------------------------------------------------------
		plt.plot(n, self.RMS)
		ax2.set_xlabel("Time (ps)")
		ax2.set_ylabel("RMSD $\AA$")
		plt.savefig( os.path.join( self.trajFolder,"analysis_mdRMSD.png") )
		if SHOW: plt.show() 
		plt.clf()  
		plt.close() 
		fig2.clf()    
		#---------------------------------------------------------------------------
		try:
			import seaborn as sns
			std_rg  = np.var(self.RG)
			std_rms = np.var(self.RMS)
			self.RG  = self.RG/std_rg
			self.RMS = self.RMS/std_rms

			g = sns.jointplot(x=self.RG,y=self.RMS,kind="kde",cmap="plasma",shade=True)
			g.set_axis_labels("Radius of Gyration $\AA$","RMSD $\AA$")
			plt.savefig( os.path.join( self.trajFolder,"rg_rmsd_biplot.png") )
			if SHOW: plt.show()
			plt.close()

		except:
			print("Error in importing seaborn package!\nSkipping biplot distribution plot!")
			pass
		#---------------------------------------------------------------------------
		fig, (ax1) = plt.subplots(nrows=1)
		plt.plot(n, self.energies)
		ax1.set_xlabel("Time (ps)")
		ax1.set_ylabel("Energy kJ/mol")
		plt.savefig(self.trajFolder+"_MDenergy.png")
		if SHOW: plt.show()
		plt.clf()
		plt.close()
		fig.clf()

	#=========================================================================
	def DistancePlots(self,RCs,SHOW=False):
		'''
		Calculate distances for the indicated reaction coordinates.
		'''	

		cnt = 0 
		for rc in RCs:
			self.RCs[str(cnt)] = [ rc, [] ]
			cnt +=1
		cnt -= 1

		print(self.RCs)

		frames = 0 
		while self.trajectory.RestoreOwnerData():
			for key in self.RCs:
				atom1 = self.RCs[key][0].atoms[0]			
				atom2 = self.RCs[key][0].atoms[1]			
				self.RCs[key][1].append( self.molecule.coordinates3.Distance(atom1, atom2) )
			frames +=1
								
		#------------------------------------------------------------------------
		# . Save the results. 
		textLog = open( self.trajFolder+"_DA.log", "w" )         
		_Text = "Frame "
		for key in self.RCs: _Text += "{} ".format(self.RCs[key][0].label[:-6])
		_Text += "\n"
		for j in range(frames):
			_Text += "{} ".format(j)
			for key in self.RCs:
				_Text += "{} ".format(self.RCs[key][1][j])			
			_Text += "\n"			
		#.------------------	
		textLog.write(_Text)
		textLog.close()
		#-------------------------------------------------------------------------
		n = np.linspace( 0, self.total_time, frames )
		#-------------------------------------------------------------------------		
		fig2, (ax2) = plt.subplots(nrows=1)		
		for key in self.RCs: 
			plt.plot( n, self.RCs[key][1], label=self.RCs[key][0].label[:-6] )
		#---------------------------------------------
		plt.xlabel("Time (ps)")
		plt.ylabel("Distances $\AA$")
		plt.legend()
		plt.savefig(self.trajFolder+"_DA.png",dpi=1000)
		plt.clf()
		plt.close()
		fig2.clf()

		self.ExtractFrames_biplot(self.RCs["0"],self.RCs["1"])

	#=========================================================================
	def Save_DCD(self):
		'''
		'''
		traj_save = os.path.join(self.trajFolder[:-4] + ".dcd")
		Duplicate(self.trajFolder,traj_save,self.molecule)
	#=========================================================================
	def Print(self):
		'''
		'''
		print("Claas printing information for Debug!")
		print( "Printing trajectory folder path: {}".format(self.trajFolder) )
		print( "RG array lenght: {}".format( len(self.RG) ) )
		print( "RMS array lenght: {}".format( len(self.RMS) ) )
		print( "RC1 most frequent:{}".format( self.rc1_MF) ) 
		print( "RC2 most frequent:{}".format( self.rc2_MF) ) 
		print( "RG most frequent:{}".format( self.rg_MF) ) 
		print( "RMS most frequent:{}".format( self.rg_MF) ) 
		
#=================================================================================

