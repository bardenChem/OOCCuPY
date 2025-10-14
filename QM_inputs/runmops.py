#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = runmops.py
import pymp
import os, glob, sys
#--------------------------------
def Run_mops(nthread):
	file_list = glob.glob("*.mop")
	with pymp.Parallel(nthread) as p:
		for i in p.range( len(file_list) ):
			os.system( "/home/venom/igor/build/mopac {}".format(file_list[i]) )
			os.system( "gzip {}.aux {}.gz".format( file_list[i][:-4],file_list[i][:-4] ) )
			try: os.system( "rm {}.aux".format(file_list[i][:-4]) )
			except: pass
#================================
if __name__ == '__main__': Run_mops( sys.argv[1] )
