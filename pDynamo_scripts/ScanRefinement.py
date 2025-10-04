#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = ScanRefinemnt.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

#================================
import pymp
from commonFunctions import *
from pMolecule import *
from pMolecule.QCModel import *

from GeometrySearcher import * 

import os, glob, sys, shutil
import numpy as np 

from pSimulation import *
from QuantumMethods import *
#================================
#**********************************************************

class ScanRefinement:
	'''
	Fazer refinamento relaxado de Scans semiempirico usando métodos de DFT.
	* Por hora, para funcionar com pySCF
	* Por hora, funciona só para caminhos unidimensionais com duas restrições de distância
	'''

	def __init__(self, _parameters):
		self.system 	   = _parameters["active_system"].system
		self.traj   	   = _parameters["source_folder"]
		self.outFolder	   = _parameters["folder"]
		self.sigma_a1_a3   = []
		self.sigma_a3_a1   = []
		self.massConstraint= False
		self.fileLists     = []
		self.forceC        = []
		self.xsize         = 0
		self.atoms         = []
		self.RCTypes       = []
		self.energies1D    = []
		self.RC1           = []
		self.RC2           = []
		self.optmizer      = _parameters["optmizer"]
		self.pureQCAtoms   = list(self.system.qcState.pureQCAtoms)
		self.logname       = os.path.join(self.outFolder,"ScanRefined.log")
		self.forceC        = [2500.0,2500.0]

		self.baseName = self.outFolder
		if not os.path.exists(self.baseName): os.makedirs(self.baseName)

		_path = os.path.join( self.traj,"")
		self.fileLists  = glob.glob(_path + "frame*.pkl")
		self.xsize      = len(self.fileLists)

		self.energies1D = pymp.shared.array( (self.xsize) , dtype='float')
		self.indexArrayX  = pymp.shared.array( (self.xsize) , dtype='uint8')
		self.GeoOptPars =   { "maxIterations":_parameters["maxIterations"]  ,
                              "rmsGradient"  :_parameters["rmsGradient"]   }

		if "force_constants"  in _parameters:
			cnt=0
			for fc in _parameters["force_constants"]:
				self.forceC[cnt] = fc
				cnt +=1
	#-------------------------------------------------------
	def SetReactionCoord(self,_RC):
		'''
		Set reaction coordinate, determining initial parameters from the atoms information
		'''
		#------------------------------------------------------------
		self.atoms.append(_RC.atoms)
		self.sigma_a1_a3.append(_RC.weight13)
		self.sigma_a3_a1.append(_RC.weight31)
		self.massConstraint         = _RC.massConstraint

		if len(_RC.atoms)   == 3: self.RCTypes.append("multipleDistance")
		elif len(_RC.atoms) == 2: self.RCTypes.append("Distance")
		else: print("Wrong Number of atoms in RC!!")

	#-----------------------------------------------------
	def RunRelaxedRefinement(self,_method,_base,_SCF_type):

		
		pySCF_pars = {"functional":_method          ,
					  "pySCF_method":_SCF_type      ,
					  "active_system":self.system   ,
					  "region":self.pureQCAtoms     , 
					  "QCcharge":0        ,
					  "method_class":"pySCF"        ,
					  "multiplicity":1              ,
					  "basis":_base                 }


		distance1 = 0.0
		distance2 = 0.0
		rmodel1   = None
		rmodel2   = None
		en0       = 0.0 

		restraints = RestraintModel( )

		for i in range(self.xsize):

			lsFrames= GetFrameIndex(self.fileLists[i][:-4])
			pySCF_pars["molden_name"] = os.path.join( self.outFolder, os.path.basename(self.fileLists[i])[:-4] + ".molden") 
			qcmol = QuantumMethods(pySCF_pars)
			qcmol.system.coordinates3 = ImportCoordinates3(self.fileLists[i])
			lsFrames= GetFrameIndex(self.fileLists[i][:-4])	
			qcmol.Set_QC_System()
			qcmol.system.DefineRestraintModel( restraints )

			if self.RCTypes[0] == "Distance":
				distance1 = qcmol.system.coordinates3.Distance( self.atoms[0][0], self.atoms[0][1])
				rmodel1   = RestraintEnergyModel.Harmonic( distance1, self.forceC[0] )
				restraint1 = RestraintDistance.WithOptions( energyModel = rmodel1, point1= self.atoms[0][0], point2= self.atoms[0][1] )
				restraints["RC1"] = restraint1
			elif self.RCTypes[0] == "multipleDistance":
				distance1  = ( qcmol.system.coordinates3.Distance( self.atoms[0][1], self.atoms[0][2]) ) *self.sigma_a3_a1[0]
				distance1 -= ( qcmol.system.coordinates3.Distance( self.atoms[0][0], self.atoms[0][1]) ) *self.sigma_a1_a3[0]
				rmodel1 = RestraintEnergyModel.Harmonic( distance1, self.forceC[0] )
				restraint1 = RestraintMultipleDistance.WithOptions( energyModel = rmodel1, distances= [ [ self.atoms[0][1],self.atoms[0][0] , self.sigma_a1_a3[0] ], [ self.atoms[0][1], self.atoms[0][2], self.sigma_a3_a1[0] ] ] )
				restraints["RC1"] = restraint1
			if self.RCTypes[1] == "Distance":
				distance2 = qcmol.system.coordinates3.Distance( self.atoms[1][0], self.atoms[1][1])
				rmodel2   = RestraintEnergyModel.Harmonic( distance2, self.forceC[1] )
				restraint2 = RestraintDistance.WithOptions( energyModel = rmodel2, point1= self.atoms[1][0], point2= self.atoms[1][1] )
				restraints["RC2"] = restraint2
			elif self.RCTypes[1] == "multipleDistance":
				distance2  = ( qcmol.system.coordinates3.Distance( self.atoms[1][1], self.atoms[1][2]) ) *self.sigma_a3_a1[1]
				distance2 -= ( qcmol.system.coordinates3.Distance( self.atoms[1][0], self.atoms[1][1]) ) *self.sigma_a1_a3[1]
				rmodel2   = RestraintEnergyModel.Harmonic( distance2, self.forceC[1] )
				restraint2 = RestraintMultipleDistance.WithOptions( energyModel = rmodel2, distances= [ [ self.atoms[1][1],self.atoms[1][0] , self.sigma_a1_a3[1] ], [ self.atoms[1][1], self.atoms[1][2], self.sigma_a3_a1[1] ] ] )
				restraints["RC2"] = restraint2		         
            #--------------------------------------------------------------------
			relaxRun = GeometrySearcher(qcmol.system, self.baseName)
			relaxRun.ChangeDefaultParameters(self.GeoOptPars)
			relaxRun.Minimization(self.optmizer)
            #--------------------------------------------------------------------
			if lsFrames[0] == 0:
				en0 = qcmol.system.Energy(log=None)
				self.energies1D[0] = 0.0
			else: self.energies1D[ lsFrames[0] ] = self.molecule.Energy(log=None) - en0 
			self.indexArrayX[ lsFrames[0] ] = lsFrames[0]
            #--------------------------------------------------------------------
			if self.RCTypes[0] == "Distance":
				self.RC1.append( qcmol.system.coordinates3.Distance( self.atoms[0][0], self.atoms[0][1]) )
			elif self.RCTypes[0] == "multipleDistance":
				self.RC1.append( qcmol.system.coordinates3.Distance( self.atoms[0][1], self.atoms[0][2]) -  qcmol.system.coordinates3.Distance( self.atoms[0][0], self.atoms[0][1]) )
			if self.RCTypes[1] == "Distance":	
				self.RC2.append( qcmol.system.coordinates3.Distance( self.atoms[1][0], self.atoms[1][1]) ) 
			elif self.RCTypes[1] == "multipleDistance":
				self.RC2.append( qcmol.system.coordinates3.Distance( self.atoms[1][1], self.atoms[1][2]) -  qcmol.system.coordinates3.Distance( self.atoms[1][0], self.atoms[1][1]) )

			Pickle( os.path.join( self.outFolder,"ScanRefined.ptGeo", "frame{}.pkl".format(lsFrames[0]) ), self.molecule.coordinates3 )

			qcmol.system.DefineRestraintModel(None)

			trajName = os.path.join( self.outFolder, "ScanRefined.dcd")
			trajpath = os.path.join( self.outFolder, "ScanRefined.ptGeo" )
			Duplicate( trajpath, trajName, self.system )  

 	#----------------------------------------------------------------------------
	def WriteLog(self):

		text = "x energy\n"
		textfile = open(self.logname)
		for i in range(self.xsize):              
			#text+= "{} {} {} {} {}\n".format( i,"0",self.RC1[i], self.RC2[i], self.energiesMatrix[i])
			text+= "{}  {}\n".format( i, self.energiesMatrix[i])
		textfile.write(textfile)
		textfile.close()
		return (self.logname)











