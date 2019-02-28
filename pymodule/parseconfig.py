#!/usr/bin/python3.6

import os

class directoryPaths:
	'''
	Parse path configure file
	USAGE : env = directoryPaths(config = './config.txt')
	'''
  
	def __init__(self, config = './config.txt'):
		self.file = config
		env = self.getConfigPaths(config)
		for key, val in zip(env.keys(), env.values()):
			if key is not None and val is not None:
				setattr(self, key, val)
	
	def getConfigPaths(self, config = './config.txt'):
		
		# check if file exists
		if(not os.path.isfile(config)):
			raise ValueError('NO SUCH FILE!')
		
		env = dict()
		with open(config, 'rt') as f:
			for line in f:
				line = line.strip()
				if line.startswith('#'):
					continue
				if not line:
					continue
				key, path = line.split('=')
				env[key] = path
		for key in env.keys():
			env.update((k, v.replace('$'+key, env[key])) for k, v in env.items())
		
		return(env)



