#!/usr/bin/env python
# -*- coding: utf-8 -*-
# mopac_module.py

#=======================================================================

#load modules

import os
from xyz_class import*
from pdb_class import*

#========================================================================

class mopac_inp:

	def __init__(self  ,
				xyzfile,
				charge ,
				multi  ,
				inpnam ,
				method):

		self.name    = xyzfile
		self.inpnam  = inpnam
		self.charge  = charge
		self.multi   = multi
		self.mult    = "Singlet"
		self.solvent = True
		self.hamilt  = method
		self.xyz     = None

		if not self.multi == 1:
			self.mult = "Doublet"

		self.xyz = xyz_parser(self.name)
		self.xyz.parse_xyz()

	def write_mop(self):

		mop_inp = open(self.inpnam,'w')
		mop_text = ''
		mop_text += '{0} 1SCF PL T=1D TIMES charge={1} {2} GEO-OK\n\n\n'.format(self.hamilt,self.charge,self.mult)

		for i in range(self.xyz.Natoms):
			mop_text +="{0}  {1}  1 {2}  1 {3} \n".format(self.xyz.AtomLabels[i],self.xyz.xCoord[i],self.xyz.yCoord[i],self.xyz.zCoord[i])

		mop_inp.write(mop_text)
		mop_inp.close()


def run_all(met):
	lists = glob.glob('*.pdb')
	for pdb in lists:
		a = protein(name=pdb,amber=True)
		a.pdb_parse(pdb)
		a.write_xyz()
	listxyz = glob.glob('*.xyz')
	for xyz in listxyz:
		a = mopac_inp(xyz,1,1,xyz[:-4]+"neutro.mop",met)
		a.write_mop()
		b = mopac_inp(xyz,2,2,xyz[:-4]+"cation.mop",met)
		b.write_mop()
		c = mopac_inp(xyz,0,2,xyz[:-4]+"anion.mop",met)
		c.write_mop()				
	script = open("run_mopac.sh",'w')
	script_text ="#!/bin/sh \n"
	listsmop = glob.glob('*.mop')
	for mop in listsmop:
		script_text += "/opt/mopac/MOPAC2016.exe {0}\n".format(mop)
	script.write(script_text)
	script.close()
