#E MultipleSequence Alignment
#Alejandro Otero Bravo
#May 14, 2018
#v0.1
#Separates taxa in an alignment based on multiple alignment metrics
#Including distance from an outgroup, base composition, and kmer frequency.
#The main point will be to separate taxa that are likely affected by LBA.
#Bug fixes: changes the alphabet function to not iterate once per record. 

from Bio import Phylo
from Bio import Align
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator as distcalc
from collections import Counter
from copy import copy
import numpy as np
from sklearn.decomposition import PCA
import logging

class emsa(Align.MultipleSeqAlignment):
	'''Modification of the class MultipleSeqAlignment to include a separate alphabet
	determination and several functions.'''
	def __init__(self, records):
		super().__init__(records)
		
	def Alphabet(self):
		#Determine the alphabet of the alignment
		nucl = set("ATCGUN-")
		dna = set("ATCG")
		rna = set("AUCG")
		prot = set("ABCDEFGHIKLMNPQRSTVWYZX-")
		self.ambig = False
		self.gapped = False
		self.alphabet = None
		letterset = set()
		for i in [x.seq for x in self]:
			for j in i:
				letterset.update(j)
		logging.debug("alphabet is: %s" % letterset)
		self.letters = list(letterset)
		if "-" in letterset:
			self.gapped = True
			letterset.remove("-")
			logging.debug("Gaps are present in the alignment.")
			logging.debug(letterset)
		if len(letterset - nucl) == 0:
			if "N" in letterset:
				self.ambig = True
				letterset.remove("N")
				logging.debug("Ns are present in the alignment.")
				logging.debug(letterset)
			if len(letterset - dna) == 0:
				self.alphabet = "DNA"
			else:
				self.alphabet = "RNA"
		elif len(letterset - prot) == 0:
			self.alphabet = "PROT"
			if 'X' in letterset:
				self.ambig = True
		#return(self.alphabet)

	#alphabet = self.alphabet()
	
	def is_named(self, id):
		'''finds the index of the sequence named id'''
		for i in range(len(self)):
			if self[i].id == id:
				return(i)
		raise Exception("Sequence with id %s not found" % id)
	
	def strip_col(self, index):
		'''removes a column from the alignment'''
		if index == 0:
			return(alignment2(self[:,1:]))
		else:
			return(alignment2(self[:,:index]+self[:,index+1:]))
	
	def strip_ambig(self, alphabet = ['A', 'C', 'T', 'G']):
		'''strips the alignment of any columns with characters not in the alphabet list'''
		rows_w_ambig = set()
		for i in self:
			for j in range(len(i.seq)):
				if i.seq[j] not in alphabet:
					rows_w_ambig.update([j])
		out = copy(self)
		count = 0
		for i in sorted(rows_w_ambig):
			out = out.strip_col(i - count)
			count += 1
		return(out)
	
	def clean(self):
		'''Removes gaps and ambiguous bases from the alignment'''
		#self.Alphabet()
		if self.alphabet is None:
			logging.info("Alphabet contains the following letters: %s" % self.letters)
			raise Exception("Alphabet was not recognized.")
		out = copy(self)
		if self.gapped:
			out.gapped = False
			out = out.strip_ambig(alphabet = [x for x in self.letters if x not in ['-']])
		if self.ambig:
			if self.alphabet == "PROT":
				out.ambig = False
				out = out.strip_ambig(alphabet = [x for x in self.letters if x not in ['X', '?']])
			elif self.alphabet in ["DNA", "RNA"]:
				out.ambig = False
				out = out.strip_ambig(alphabet = [x for x in self.letters if x not in ['N']])
		return(out)
	
	def without_index(self, index):
		'''returns the same alignment excluding the sequence of the index provided
		if a list of indices is provided it will exclude all of them.'''
		if type(index) == int:
			return([x for i,x in enumerate(self) if i != index])
		elif type(index) == str:
			return([x for i,x in enumerate(self) if i != self.is_named(index)])
		elif type(index) == list:
			exclude = []
			for i in index:
				if type(i) == str:
					exclude.append(self.is_named(i))
				elif type(i) == int:
					exclude.append(i)
			return([x for i,x in enumerate(self) if i not in index])
		else:
			raise Exception("index not an integer, string or list")

	def distance_from(self, sequence, model = 'identity'):
		'''returns a dictionary with the keys being sequence ids and the distance from the
		sequence provided'''
		calculator = distcalc(model)
		distance = {}
		for i in range(len(self)):
			distance[self[i].id] = calculator._pairwise(sequence, self[i])
		return(distance)
	
	def composition(self):
		'''returns a dictionary with keys being sequence ids and a set of the sequence
		alphabet with their frequencies.'''
		comp_dict = {}
		for i in self:
			compos = Counter()
			for j in i.seq:
				compos[j] += 1
			comp_dict[i.id] = compos
		return(comp_dict)
	
	def kmer_composition(self, kmer = 3):
		'''returns a dictionary with keys being sequence ids and a set of the sequence
		kmers with their frequencies'''
		composition_dict = {}
		for i in self:
			compos = Counter()
			for j in range(len(i.seq) - kmer):
				compos[str(i.seq[j:j + kmer])] += 1
			composition_dict[i.id] = compos
		return(composition_dict)
	
	def stat_array(self, kmer = 2, outgroup = 1):
		'''returns an array of statistics on each sequence including the distance, base 
		composition, and kmer frequency.'''
		dist = self.distance_from(self[outgroup])
		comp = self.composition()
		kmer = self.kmer_composition(kmer = kmer)
		distance_list = []
		composition_dict = {}
		for item in self.letters:
			composition_dict[item] = []
		kmer2_set = set() #We identify the number of total kmers in the sample
		for key in [x.id for x in self]:
			kmer2_set.update(kmer[key])
		kmer2_dict = {} #And now for each we create a list in a dictionary
		for kmer_val in kmer2_set:
			kmer2_dict[kmer_val] = []
		for key in [x.id for x in self]:
			distance_list.append(dist[key])
			for item in self.letters:
				composition_dict[item].append(comp[key][item])
			for kmer_val in kmer2_dict:
				kmer2_dict[kmer_val].append(kmer[key][kmer_val])
		output_tuple = tuple([distance_list]) + tuple([composition_dict[x] for x in self.letters]) + tuple([kmer2_dict[x] for x in kmer2_set]) 
		end_array = np.column_stack(output_tuple)
		return(end_array)
	
	def get_pca(self, array = None, groups = 1):
		'''Returns a list with the first item being a dictionary with Principal Component 
		values with each id as a key and the second item being the proportion of total 
		variance explained by this component. The array from stat_array is used as default
		but a separate numpy array can be provided given that rows match the order of ids.
		'''
		if array is None:
			stripped = self.strip_ambig()
			array = stripped.stat_array()
		pca_object = PCA(n_components = groups)
		pca_object.fit(array)
		pca_results = pca_object.fit_transform(array)
		ids = [x.id for x in self]
		values = {}
		for i in range(len(ids)):
			values[ids[i]] = float(pca_results[i])
		return([values, pca_object.explained_variance_ratio_] )
	
	def get_different(self, pca = None, cutoff = 0):
		'''Separates taxa into two lists of similar groups. To avoid recalculation of the 
		PCA, values can be given as a dictionary (see get_pca).'''
		if pca is None:
			pca = self.get_pca()[0]
		else:
			if type(pca) == list:
				pca = pca[0]
		values = pca
		group1 = []
		group2 = []
		for key in values:
			if values[key] >= cutoff:
				group1.append(key)
			else:
				group2.append(key)
		return([group1, group2])
		



