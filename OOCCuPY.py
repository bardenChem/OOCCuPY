#/usr/bin/env python
# -*- coding: utf-8 -*-
#OOCCuPY.py

import argparse,sys,os
import subprocess 

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
		if self.args.tests:  
			# Use subprocess to capture output
			result = subprocess.run(
				["python3", "Tests/pDynamoWrapper/Run_All.py"], 
				stdout=subprocess.PIPE, 
				stderr=subprocess.STDOUT,  # Merge stderr with stdout
				text=True
			)
			print(result.stdout)  # This goes through our Tee
            
		elif self.args.test_number: 
			test_file = "Tests/pDynamoWrapper/test_{:02d}.py".format(self.args.test_number)
			# Use subprocess to capture output
			result = subprocess.run(
				["python3", test_file], 
				stdout=subprocess.PIPE, 
				stderr=subprocess.STDOUT,
				text=True
			)
			print(result.stdout)


	def MD_prep_handler(self):
		'''
		'''

	def QM_input_handler(self):
		'''
		'''

	def run_testes(self,_number):
		'''
		'''
#=============================================================
if __name__=="__main__":
	parser = argparse.ArgumentParser(
						prog="OOCCuPy",
						description="Object Oriented Computational Chemistry ")

	parser.add_argument("-o", "--output", default="output.txt", help="Output file (default: output.txt)")
	parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode")
	
	subparsers 		 = parser.add_subparsers(dest="command", help="Comandos disponíveis")
	pdynamo_parser 	 = subparsers.add_parser("pDynamo", help="Module pDynamo Wrapper ")
	mdtools_parser 	 = subparsers.add_parser("mdtools", help="Module to prep and handle molecular dynamics ")
	qm_inputs_parser = subparsers.add_parser("QM_inputs", help="Module to create input for QM softwares")

	pdynamo_parser.add_argument("--tests",action="store_true",help="Run All Tests!!")
	pdynamo_parser.add_argument("--test",type=int,dest="test_number",help="Run specif test!!")
	pdynamo_parser.add_argument("--input",action="store_true",help="Run specif test!!")
	
	args = parser.parse_args()
	print(args)

	if args.verbose: print("Verbose mode enabled!")

	RunOOCCuPy = Interface(args)
	if 		args.command == "pDynamo": RunOOCCuPy.pDynamoWrapper_handler()
	elif 	args.command == "mdtools": RunOOCCuPy.MD_prep_handler()
	elif 	args.command == "QM_inputs": RunOOCCuPy.QM_input_handler() 







