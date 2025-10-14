#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = pDynamoWrapper.py

#--------------------------------------------------------------
import os, glob, sys

from commonFunctions import *
from SimulationSystem import SimulationSystem
from Simulation import *
from Analysis import *

#==============================================================
class Wrapper:
	'''
	Set Wrapper and control the library of the presetted simulations
	'''

	def __init__(self, _projectFolder=None):
		'''
		Initialize project folder and paths.
		'''
		self.activeSystem    = None 
		self.system_historic = []
		self.projectFolder   = os.getcwd()

		if not _projectFolder: self.projectFolder = os.getcwd()
		else: self.projectFolder = os.path.join( os.getcwd(), _projectFolder )		
		if not os.path.exists(self.projectFolder): os.makedirs(self.projectFolder)
		

	#-----------------------------------------
	@classmethod
	def From_Input(selfClass,_inputFile,_projectFolder_=None):
		'''
		Reads input and automate system setting, simulation and analysis.
		'''
		self            = selfClass(_projectFolder=_projectFolder_)

		RUN_SIMULATIONS = 0
		RUN_ANALYSIS    = 0 
		SET_CRD_NMB     = 0

		Spherical_Pruning = False
		Fixed_atoms = False
		reactions_crds = []
		_parameters = {}
		simulations = []
		_parameters["mass_constraints"] = []
		_parameters["crd_labels"] = []
		save_name = None

		inpFile = open(_inputFile,"r")
		for line in inpFile:
			lines = line.split()
			if len(lines) > 0:
				if  lines[0] == "#INPUT_TYPE":
					_parameters["Input_Type"] = lines[1]
					if lines[1] == "geometry":
						_parameters["crd_file"] = lines[2]
					elif lines[1] == "amber" or lines[1] == "gromacs":
						_parameters["top_file"] = lines[2]
						_parameters["crd_file"] = lines[3]
					elif lines[1] == "pkl":
						_parameters["pkl_file"] = lines[2]
					elif lines[1] == "protein":
						_parameters["pdb_file"] = lines[2]
				elif lines[0] == "#SPHERICAL_PRUNNING":
					if not lines[1] == "no":
						_parameters["spherical_prune"] = lines[1] # put the patterns to center the selection
						_parameters["spherical_prune_radius"] = float(lines[2])
				elif lines[0] == "#FIXED_ATOMS":
					if not lines[1] == "no":
						_parameters["set_fixed_atoms"] = lines[1] # put the patterns to centert the selection
						_parameters["free_atoms_radius"] = float(lines[2])
				elif lines[0] == "#SET_REACTION_CRD":
					if not lines[1] == "no":
						SET_CRD_NMB+=1
						_parameters["type_rc"+str(SET_CRD_NMB)]     = lines[1]
						_parameters["mass_constraints"].append(lines[2])
						_parameters["atoms_rc"+str(SET_CRD_NMB)] = []
						for i in range(0, int(lines[3])):
							_parameters["atoms_rc"+str(SET_CRD_NMB)].append(lines[i+4])
				elif lines[0] == "#SET_INITIAL_CRD":
					if not lines[1] == "no": _parameters["set_initial_crd"] = lines[1]
				elif lines[0] == "#SET_QMMM_REGION":
					if not lines[1] == "no":
						_parameters["set_qc_region"] = "yes"
						if lines[1] == "#residue_patterns":
							list_res = []
							for i in range( 0,int(lines[2]) ):
								list_res.append(lines[i+3]) 
							_parameters["residue_patterns"] = list_res						
						elif lines[1] == "#center_atom":
							_parameters["center_atom"] = lines[2]
							_parameters["radius"] = float(lines[3])
				elif lines[0] == "#SET_ENERGY_MODEL_QM":
					if not lines[1] == "no": 
						_parameters["method_class"] = lines[1]
						_parameters["set_energy_model"] = "QM"
				elif lines[0] == "#HAMILTONIAN": _parameters["Hamiltonian"] = lines[1]
				elif lines[0] == "#RUN_SIMULATION":	_parameters["simulation_type"] = simulations.append(lines[1])
				elif lines[0] == "#QC_CHARGE": _parameters["QCcharge"] = int(lines[1])
				elif lines[0] == "#MULTIPLICITY": _parameters["multiplicity"] = int(lines[1])
				elif lines[0] == "#ORCA_METHOD": _parameters["orca_method"] = lines[1] 
				elif lines[0] == "#FUNCTIONAL": _parameters["functional"] = lines[1] 
				elif lines[0] == "#SOFTWARE": _parameters["Software"] = lines[1]
				elif lines[0] == "#BASIS": _parameters["basis"] = lines[1]
				elif lines[0] == "#PYSCF_METHOD": _parameters["pySCF_method"] = lines[1]
				elif lines[0] == "#METHODS_SMO":
					_parameters["methods_lists"] = []
					for idx in range( 1,len(lines) ):
						_parameters["methods_lists"].append(lines[idx])
				elif lines[0] == "#MOPAC_KEYS": 
					_parameters["mopac_keywords"] = []
					for idx in range( 1,len(lines) ):
						_parameters["mopac_keywords"].append(lines[idx])
				elif lines[0] == "#XNBINS": _parameters["xnbins"] = int(lines[1])
				elif lines[0] == "#YNBINS": _parameters["ynbins"] = int(lines[1])
				elif lines[0] == "#MAX_ITERATIONS": _parameters["maxIterations"] = int(lines[1])
				elif lines[0] == "#RMS_GRADIENT": _parameters["rmsGradient"] = float(lines[1])
				elif lines[0] == "#OPT_ALG":  _parameters["optmizer"] = lines[1]
				elif lines[0] == "#DINCRE_RC1": _parameters["dincre_rc1"] = float(lines[1])
				elif lines[0] == "#DINCRE_RC2": _parameters["dincre_rc2"] = float(lines[1])
				elif lines[0] == "#NDIMS":  _parameters["ndim"] = int(lines[1])
				elif lines[0] == "#X_NSTEPS": _parameters["nsteps_rc1"] = int(lines[1])
				elif lines[0] == "#Y_NSTEPS": _parameters["nsteps_rc2"] = int(lines[1])
				elif lines[0] == "#MD_METHOD": _parameters["MD_method"] = lines[1]
				elif lines[0] == "#PRESSURE": _parameters["pressure"] = float(lines[1])
				elif lines[0] == "#PRESSURE_COUPLING": _parameters["pressure_coupling"] = int(lines[1])
				elif lines[0] == "#PRESSURE_CONTROL": _parameters["pressure_control"] = lines[1]
				elif lines[0] == "#TEMPERATURE_SCALE_OPTION": _parameters["temperature_scale_option"] = lines[1]
				elif lines[0] == "#TEMPERATURE_SCALE": _parameters["temperature_scale"] = float(lines[1])
				elif lines[0] == "#TIME_STEP": _parameters["timeStep"] = float(lines[1])
				elif lines[0] == "#COLLISION_FREQ": _parameters["coll_freq"] = float(lines[1])
				elif lines[0] == "#INITIAL_TEMPERATURE": _parameters["start_temperature"] = float(lines[1])
				elif lines[0] == "#EQUILIBRATION_STEPS": _parameters["equilibration_nsteps"] = int(lines[1])
				elif lines[0] == "#HEATING_STEPS": _parameters["heating_nsteps"] = int(lines[1])
				elif lines[0] == "#PRODUCTION_STEPS": _parameters["production_nsteps"] = int(lines[1])
				elif lines[0] == "#SAMPLING_HEATING": _parameters["heating_sampling"] = int(lines[1])
				elif lines[0] == "#SAMPLING_EQUILIBRATION": _parameters["sampling_equilibration"] = int(lines[1])
				elif lines[0] == "#SAMPLING_PRODUCTION": _parameters["sampling_production"] = int(lines[1])
				elif lines[0] == "#FORCE_CONSTANTS":
					_parameters["force_constants"] = [] 
					for idx in range(1,len(lines)):
						_parameters["force_constants"].append( float(lines[idx]) )
				elif lines[0] == "#CRD_INPUT": _parameters["crd_format"] = lines[1]
				elif lines[0] == "#XNWINDOWS": _parameters["xwindows"] = int(lines[1])
				elif lines[0] == "#YNWINDOWS": _parameters["ywindows"] = int(lines[1])
				elif lines[0] == "#CYCLES": _parameters["cycles"] = int(lines[1])
				elif lines[0] == "#MODE": _parameters["mode"] = int(lines[1])
				elif lines[0] == "#FRAMES": _parameters["frame"] = int(lines[1])
				elif lines[0] == "#BINS": _parameters["traj_bins"] = int(lines[1])
				elif lines[0] == "#INITIAL_CRD": _parameters["init_coord"] = lines[1]
				elif lines[0] == "#FINAL_CRD": _parameters["final_coord"] = lines[1]
				elif lines[0] == "#SPRING_FORCE": _parameters["spring_force_constant"] = float(lines[1])
				elif lines[0] == "#FIXED_NEB": _parameters["fixed_terminal_images"] = lines[1]
				elif lines[0] == "#RMS_GROWING_STRING": _parameters["RMS_growing_intial_string"] = lines[1]
				elif lines[0] == "#RESTART": _parameters["restart"] = lines[1]
				elif lines[0] == "#ANALYSIS_ONLY": _parameters["analysis_only"] = lines[1]
				elif lines[0] == "#MAX_NUM_OF_THREADS": _parameters["NmaxThreads"] = int(lines[1])
				elif lines[0] == "#TEMPERATURE": _parameters["temperature"] = float(lines[1])
				elif lines[0] == "#LOG_FREQUENCY": _parameters["log_frequency"] = int(lines[1])
				elif lines[0] == "#SAVE_FORMAT": _parameters["save_format"] = lines[1]
				elif lines[0] == "#TRAJECTORY_NAME": _parameters["trajectory_name"] = lines[1]
				elif lines[0] == "#SEED": _parameters["seed"] = int(lines[1])
				elif lines[0] == "#RELAX": _parameters["relax"] = lines[1]
				elif lines[0] == "#REVERSE_RC1": _parameters["reverse_rc1"] = lines[1]
				elif lines[0] == "#REVERSE_RC2": _parameters["reverse_rc2"] = lines[1]
				elif lines[0] == "#CORRECT_QMMM_CHARGES": _parameters["correct_QMMM_charge"] = lines[1]
				elif lines[0] == "#SOURCE_FOLDER": _parameters["source_folder"] = lines[1]
				elif lines[0] == "#FOLDER": _parameters["folder"] = lines[1]
				elif lines[0] == "#FIGSIZE": _parameters["fig_size"] = [ int(lines[1]), int(lines[2]) ]
				elif lines[0] == "#CONTOUR_LINES": _parameters["contour_lines"] = int(lines[1])
				elif lines[0] == "#XLIM": _parameters["xlim"] = [ float(lines[1]), float(lines[2]) ]
				elif lines[0] == "#YLIM": _parameters["ylim"] = [ float(lines[1]), float(lines[2]) ]
				elif lines[0] == "#CRD_LABEL_1": _parameters["crd_labels"].append(lines[1])
				elif lines[0] == "#CRD_LABEL_2": _parameters["crd_labels"].append(lines[1])
				elif lines[0] == "#SAVE_NAME": save_name = lines[1]
				elif lines[0] == "#SCRATCH": _parameters["scratch"] = lines[1] 
				elif lines[0] == "#REFINE_METHODS": _parameters["refine_methods"] = lines[1] 

		
		_parameters["set_reaction_crd"] = SET_CRD_NMB			
		self.Set_System(_parameters)

		for sim in simulations:
			_parameters["simulation_type"] = sim
			self.Run_Simulation(_parameters)

		self.SaveSystem(save_name)
		inpFile.close()		

	#-----------------------------------------
	def Set_System(self,_parameters):
		'''
		Sets up the simulation system based on the provided parameters.
		'''
		_system4load =  SimulationSystem()

		mass_constraints = []

		if "mass_constraints" in _parameters:
			for mc in _parameters["mass_constraints"]:
				if mc == "yes": mass_constraints.append(True)
				elif mc == "no": mass_constraints.append(False)

		_debug_ok = False
		if "DEBUG" in _parameters: _debug_ok = True 
		if "Input_Type" not in _parameters:	raise KeyError("Missing required parameter: Input_Type")

		#Load system based on Input_Type using a dictionary for clarity
		input_methods = {
			"geometry": SimulationSystem.From_Coordinates,
			"amber": SimulationSystem.From_AMBER,
			"gromacs": SimulationSystem.From_Gromacs,
			"pkl": SimulationSystem.From_PKL,
			"protein":SimulationSystem.Protein_From_Coordinates,
		}

		try:
			input_type = _parameters["Input_Type"]
			load_function = input_methods[input_type]
			if input_type == "geometry":
				_system4load = load_function(_parameters["crd_file"])
			elif input_type == "amber":
				_system4load = load_function(_parameters["top_file"], _parameters["crd_file"])
			elif input_type == "gromacs":
				_system4load = load_function(_parameters["top_file"], _parameters["crd_file"])
			elif input_type == "pkl":
				_system4load = load_function(_parameters["pkl_file"])
			elif input_type == "protein":
				_system4load = load_function(_parameters["pdb_file"])
			else: raise ValueError(f"Unsupported Input_Type: {input_type}")

		except (KeyError, TypeError) as e:
		# Handle potential issues with missing keys or invalid values
			raise ValueError(f"Error loading system: {e}") from e  # Chain the original exception
	
		self.activeSystem = _system4load
		#------------------------------------
		if "spherical_prune" in _parameters:			
			self.activeSystem.Spherical_Pruning(_parameters["spherical_prune"],float(_parameters["spherical_prune_radius"]))
		if "set_fixed_atoms" in _parameters:
			self.activeSystem.Setting_Free_Atoms(_parameters["set_fixed_atoms"],float(_parameters["free_atoms_radius"]))
		if "set_reaction_crd" in _parameters:
			for rc in range(0,_parameters["set_reaction_crd"]):				
				self.activeSystem.Set_Reaction_crd( _parameters["atoms_rc"+str(rc+1)],_parameters["type_rc"+str(rc+1)],mass_constraints[rc])
		if "set_initial_crd" in _parameters:
			if ( _parameters["set_initial_crd"][-4:] ) == ".pkl":
				try:				
					self.activeSystem.system.coordinates3 = Unpickle(_parameters["set_initial_crd"])[0]
				except:
					self.activeSystem.system.coordinates3 = Unpickle(_parameters["set_initial_crd"])
			else: self.activeSystem.system.coordinates3 = ImportCoordinates3(_parameters["set_initial_crd"])
		if "set_qc_region" in _parameters:
			_residue_list = []
			_centerAtom = None
			_radius = None
			if "residue_patterns" in _parameters: _residue_list = _parameters["residue_patterns"]
			if "center_atom"      in _parameters: _centerAtom   = _parameters["center_atom"]
			if "radius"           in _parameters: _radius       = _parameters["radius"]
			if _parameters["set_qc_region"] == "yes":
				self.activeSystem.Set_QCMM_Region(_residue_list,_centerAtom,_radius)
		if "set_energy_model" in _parameters:
			if _parameters["set_energy_model"] == "QM":
				self.activeSystem.Set_QC_Method(_parameters,_DEBUG=_debug_ok)

		#-----------------------------------
		self.system_historic.append(_system4load)
		self.activeSystem.Check()


	#-----------------------------------------
	def Run_Simulation(self,_parameters):
		'''
		'''
		_parameters["active_system"] = self.activeSystem
		_parameters["project_folder"]= self.projectFolder
		self.system_historic.append(self.activeSystem) 
		#----------------------------------------------- 
		_Run_Simulation   = Simulation(_parameters)
		self.activeSystem.system = _Run_Simulation.Execute()

	#-----------------------------------------
	def Run_Analysis(self,_parameters):
		'''
		'''
		_parameters["active_system"] = self.activeSystem
		self.system_historic.append(self.activeSystem) 
		#--------------------------------------------
		_Analysis = Analysis(_parameters)
		_Analysis.Execute()
		
	#-----------------------------------------
	def PrintSystems(self):
		'''
		Method to print the summary of the loaded systems 
		'''
		print("There are {} loaded systems".format( len(self.system_historic) ) )
		ctn = input("Type any key to print the Summary of the Systems, or 'N' to cancel this")
		if not ctn == "N":
			if len(self.system_historic) > 0:
				for system in self.system_historic:
					system.system.Summary()
					print("***************************************************")
				print("Now, printing the current system Summary:")
				self.activeSystem.system.Summary()
		
			elif len(self.system_historic) == 1: print("There is only the current System loaded!\n Printing its information below!")
			else:                                print( "There are no loaded systems!")

	#.-------------------------------------------------------------------------
	def SaveProject(self):
		'''
		The complete version of this function intends to save in pkl and another coordinate format
		the systems and trajectories worked in this simulations
		Though, in the current state only will save the current system to a pkl file
		'''
		_baseName = self.activeSystem.basename
	
		savePathPkl = os.path.join(self.projectFolder,_baseName+".pkl")
		savePathPdb = os.path.join(self.projectFolder,_baseName+".pdb")
	
		for _System in self.system_historic:
			_label = _System.label
			savePathPkl = os.path.join(savePathPkl, _baseName+"_"+_label+".pkl")
			Pickle( savePathPkl, _System )        
    #.-------------------------------------------------------------------------
	def SaveSystem( self, _cname=None):
		'''
		'''
	
		_baseName = self.activeSystem.baseName
		if _cname:
			savePathPkl = os.path.join(self.projectFolder,_cname+".pkl")
			savePathPdb = os.path.join(self.projectFolder,_cname+".pdb")
		else:
			savePathPkl = os.path.join(self.projectFolder,_baseName+".pkl")
			savePathPdb = os.path.join(self.projectFolder,_baseName+".pdb")
	
			i = 0;
			while os.path.exists(savePathPdb):
				i += 1
				savePathPdb = os.path.join(self.projectFolder,_baseName+".pdb")
				savePathPdb = savePathPdb[:-4] + "_#{}.pdb".format(i)
			while os.path.exists(savePathPkl):
				i += 1
				savePathPkl = os.path.join(self.projectFolder,_baseName+".pkl")
				savePathPkl = savePathPkl[:-4] + "_#{}.pkl".format(i)
		#----------------------------------------------------------------
		self.activeSystem.system.Summary()
		ExportSystem( savePathPkl,self.activeSystem.system )
		ExportSystem( savePathPdb,self.activeSystem.system )
#==============================================================
if __name__=="__main__":
	'''
	'''
	run_script = Wrapper.From_Input(sys.argv[1],sys.argv[2])
	
#======================================================================================================#
#======================================END OF FILE=====================================================#
#======================================================================================================#