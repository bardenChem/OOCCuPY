#/usr/bin/env python
# -*- coding: utf-8 -*-
#OOCCuPY.py

import argparse

class Interface:
	'''
	Class to handle the functionalities of OOCCuPy
	'''

	def __init__(self):
		'''
		'''

	def pDynamoWrapper_handler(self):
		'''
		'''


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
	parser = argparse.ArgumentParser(description="OOCCuPy Software")
	parser.add_argument("input", help="Input file")
	parser.add_argument("-o", "--output", default="output.txt", help="Output file (default: output.txt)")
	parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode")
	
	subparsers 		 = parser.add_subparsers(dest="command", help="Comandos disponíveis")
	pdynamo_parser 	 = subparsers.add_parser("pDynamoWrapper", help="Module pDynamo scripts ")
	mdtools_parser 	 = subparsers.add_parser("mdtools", help="Module to prep and handle molecular dynamics ")
	qm_inputs_parser = subparsers.add_parser("QM_inputs", help="Module to create input for QM softwares")
	
	if args.verbose: print("Verbose mode enabled!")





