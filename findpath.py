import sys
import os

def pathfinder():
	try:
		abspath=os.getcwd()
		rootpath=os.path.join(abspath[:abspath.index('N2D')],'N2D')
	except:
		ppath=sys.path
		for p in ppath:
			if 'N2D' in p:
				rootpath=os.path.join(p[:p.index('N2D')],'N2D')
	return rootpath