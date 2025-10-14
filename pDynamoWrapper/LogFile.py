#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = LofFile.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

from pCore import *
from datetime import datetime
from timeit import default_timer as timer
import os

#---------------------------------------------------------
class LogFile:
	'''
	'''
	def __init__(self):
		'''
		'''
		self.data = {}
		self.data["Log_Type"]   = "None"
		self.data["dimensions"] = 0
		self.data["QM_methods"] = []

		pass 

	#-----------------------------------------------------
	@classmethod
	def From_Reaction_Scan(self,file_name,data_object):
		'''
		'''
		pass 

	@classmethod
	def From_Energy_Refinement(self,file_name,data_object):
		'''
		'''
		pass 
	#------------------------------------------------------
	@classmethod
	def From_Molecular_Dynamics(self,file_name,data_object):
		'''
		'''
		pass
	#------------------------------------------------------
	@classmethod
	def From_Umbrella_Sampling(self,file_name,data_object):
		'''
		'''
		pass 

	#------------------------------------------------------
	@classmethod 
	def From_PMF_Analysis(self,file_name,data_object):
		'''
		'''
		pass 

#============================================================






