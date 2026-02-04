#/usr/bin/env python
# -*- coding: utf-8 -*-
#OOCCuPY.py

import argparse,sys,os
import subprocess 
from pDynamoWrapper.pDynamoWrapper import Wrapper
from pathlib import Path

# Add config utilities
try:
    from config import get_config, find_data_file
except ImportError:
    # For direct execution
    sys.path.insert(0, str(Path(__file__).parent))
    from config import get_config, find_data_file


#=======================================================
class Tee:
    def __init__(self, filename, mode='a'):
        self.file = open(filename, mode)
        self.stdout = sys.stdout  # Save original stdout

    def write(self, text):
        self.file.write(text)
        self.stdout.write(text)  # Print to screen

    def flush(self):  # Required for file flushing
        self.file.flush()
        self.stdout.flush()

    def close(self):
        self.file.close()

#========================================================

class Interface:
	'''
	Class to handle the functionalities of OOCCuPy
	'''

	def __init__(self,args):
		'''
		'''
		self.args  = args
		self.tee   = Tee(self.args.output, 'w')
		sys.stdout = self.tee

	def __del__(self):
		'''Clean up when object is destroyed'''
		if hasattr(self, 'tee'):
			sys.stdout = self.tee.stdout  # Restore original stdout
			self.tee.close()

	def pDynamoWrapper_handler(self):
		'''
		'''
		if not self.args.proj_folder: 
			self.args.proj_folder = None


		# Make sure we print something through the Tee first to verify it's working
		print(f"Starting pDynamoWrapper_handler - output will be saved to {self.args.output}")

		config = get_config()

		if self.args.tests:  
        	# Find test directory
        	test_dir = config.get_test_data_path("pDynamoWrapper")
        	print(f"Looking for tests in: {test_dir}")
        
        	# Try to find Run_All.py
        	run_all_path = find_data_file("Run_All.py", "pDynamoWrapper")
        	if run_all_path:
            	self._run_subprocess_real_time(["python3", str(run_all_path)])
        	else:
        	    print(f"⚠ Run_All.py not found")
            	print(f"Searching in: {test_dir}")
            
            	# Run individual tests
            	test_files = list(test_dir.glob("test_*.py"))
            	if test_files:
                	print(f"Found {len(test_files)} test files:")
                	for tf in sorted(test_files):
                	    print(f"  - {tf.name}")
                
                	# Option: Run all tests sequentially
                	for test_file in sorted(test_files):
                	    print(f"\n{'='*50}")
                	    print(f"Running: {test_file.name}")
                	    print('='*50)
                	    self._run_subprocess_real_time(["python3", str(test_file)])
            	else:
                	print("No test files found!")
                
    	elif self.args.test_number: 
        	# Find specific test
        	test_file = find_data_file(f"test_{self.args.test_number:02d}.py", "pDynamoWrapper")
        	if test_file:
            	self._run_subprocess_real_time(["python3", str(test_file)])
        	else:
         	   print(f"⚠ Test {self.args.test_number} not found")
            
    	else:
        	if self.args.inp_file:
            	run_input = Wrapper.From_Input(self.args.inp_file, self.args.proj_folder)



	def MD_prep_handler(self):
		'''
		'''

	def QM_input_handler(self):
		'''
		'''

	def run_testes(self,_number):
		'''
		'''

#================================================================================
def main():
	parser = argparse.ArgumentParser(
						prog="OOCCuPy",
						description="Object Oriented Computational Chemistry ")

	parser.add_argument("-o", "--output", default="output.txt", help="Output file (default: output.txt)")
	parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode")
	
	subparsers 		 = parser.add_subparsers(dest="command", help="Comandos disponíveis")
	pdynamo_parser 	 = subparsers.add_parser("pDynamo", help="Module pDynamo Wrapper ")
	mdtools_parser 	 = subparsers.add_parser("mdtools", help="Module to prep and handle molecular dynamics ")
	qm_inputs_parser = subparsers.add_parser("QM_inputs", help="Module to create input for QM softwares")
	config_parser    = subparsers.add_parser("config", help="Configuration management")
	config_group     = config_parser.add_mutually_exclusive_group()

	pdynamo_parser.add_argument("--tests",action="store_true",help="Run All Tests!!")
	pdynamo_parser.add_argument("--test",type=int,dest="test_number",help="Run specif test!!")
	pdynamo_parser.add_argument("--input",type=str,dest="inp_file",help="Run pDynamo from input file!")
	pdynamo_parser.add_argument("--proj_folder",type=str,dest="proj_folder")
	
	config_group.add_argument("--show", action="store_true", help="Show configuration")
	config_group.add_argument("--paths", action="store_true", help="Show paths")
	config_group.add_argument("--setup", action="store_true", help="Setup environment")

	args = parser.parse_args()
	if args.verbose: print("Verbose mode enabled!")

	RunOOCCuPy = Interface(args)
	if 		args.command == "pDynamo": RunOOCCuPy.pDynamoWrapper_handler()
	elif 	args.command == "mdtools": RunOOCCuPy.MD_prep_handler()
	elif 	args.command == "QM_inputs": RunOOCCuPy.QM_input_handler() 
	elif args.command == "config":
    	if args.show:
    	    config = get_config()
    	    config.show()
    	elif args.paths:
        	config = get_config()
        	print("\nOOCCuPY Paths:")
        	print("="*40)
        	print(f"Config: {config.config_dir}")
        	print(f"Tests: {config.get_test_data_path()}")
        	print(f"Examples: {config.get_examples_path()}")
        
        	# Show package paths if available
        	pkg_tests = config.get('paths.package_tests')
        	pkg_examples = config.get('paths.package_examples')
        	if pkg_tests:
        	    print(f"Package Tests: {pkg_tests}")
        	if pkg_examples:
        	    print(f"Package Examples: {pkg_examples}")
    	elif args.setup:
        	from config import setup_environment
        	setup_environment()


#=============================================================
if __name__=="__main__":
	main()







