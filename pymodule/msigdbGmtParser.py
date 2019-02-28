## mSigDB gmt gene pathway file parser

class msigdbGmtParser:
	"""
	Read msigdb pathway gmt files into class
	- fname : gmt file from msigDB (ex 'c5.bp.v6.1.entrez.gmt')
	- idf : gene pathway simplified id format (ex 'go{:04}')
	"""
	
	def __init__(self, fname, idf='id{:04}'):
		self.fname = fname
		self.idf = idf
		self.load_gmt()
	
	def __len__(self):
		return len(self.gmt)
	
	def load_gmt(self):
		with open(self.fname, 'r') as f:
			gmt = f.readlines()
		self.gmt = gmt
	
	def makeIDs(self):
		gsids = dict()
		for i, line in enumerate(self.gmt):
			ID = self.idf.format(i)
			gsids[ID] = line.strip().split('\t')[0]
		return gsids
	
	def pathway2gene(self):
		gsets = dict()
		for i, line in enumerate(self.gmt):
			ID = self.idf.format(i)
			gsets[ID] = line.strip().split('\t')[2:]
		return gsets
	
	def gene2pathway(self):
		gsetrev = dict()
		gsets = self.pathway2gene()
		for k,v in gsets.items():
			for x in v:
				gsetrev.setdefault(x,[]).append(k)
		return gsetrev


"""	
Usage :
obj = msigdbGmtParser(fname, idf=idf)
gsids = obj.makeIDs()
gsets = obj.pathway2gene()
gsetrev = obj.gene2pathway()
"""
