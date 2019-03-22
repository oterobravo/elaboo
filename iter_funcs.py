#Functions for the iteration of the algorithm
#Alejandro Otero Bravo
#v0.4
#Modules are defined here.

#Disclaimer: The following is a preliminary version that has not been thoroughly tested.

from treefuns import *
from emsa import *
from Bio.Phylo.Applications import RaxmlCommandline
from Bio.Phylo.Applications import _Fasttree
import glob
import os
import argparse
import random
import string 
import logging
import matplotlib.pyplot as hs

def elaboo_parser():
	parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter)                                               

	parser.add_argument("-a", "--alignment", type = str, required = True, help = "Input alignment")
	parser.add_argument("-f", "--input_format", type = str, default = "fasta", help = "Input alignment format, default: fasta")
	parser.add_argument("-b", "--bootstraps", type = int, default = 10, help = "Number of bootstrap replicates per taxon in question")
	parser.add_argument("-o", "--outgroup", type = str, help = "Outgroup from the alignment")
	parser.add_argument("-s", "--taxa_stats_file", default = None, type=str, help = "File including values for PCA\n(must conform to numpy array standards).")
	parser.add_argument("-p", "--pca_threshold", type=int, default = 0, help = "Threshold for separating groups using the principal component")
	parser.add_argument("-m", "--model", type = str, default = None, help = "Sequence evolution model for tree building.\n(default GTRCAT for nucleotide, PROTCATDAYHOFF for protein)")
	parser.add_argument("-g", "--taxa_list", type = str, help = "Text file including the list of taxa to be analyzed separately.")
	parser.add_argument("-t", "--trees", type = str, default = 'resulting_trees.tre', help = "Newick file of trees to be used for consolidation.")
	parser.add_argument("-d", "--distribution", type = str, default = None, help = "Plot distribution of PCA values on given file.")
	parser.add_argument("-c", "--clean", type = bool, default = False, help = "File provided has no gaps nor ambiguities.")
	parser.add_argument("-y", "--algorithm", type = str, default = 'a', help = '''Which method to take.
		a: (default) do all from just alignment.
		b: Only calculate statistics for the alignment.
		c: List potentially problematic taxa from alignment statistics.
		d: Generate trees for each problem taxa.
		e: Split and generate trees from a separate array.
		f: Generate trees for my list of problem taxa.
		g: Generate consensus tree from trees provided in file resulting_trees.tre
		h: Split taxa based on precalculated summary statistics
		i: Build trees for a given set of taxa''')
	parser.add_argument("-e", "--tree_build", type = str, default = "fast", help = '''Method for tree building.
		fast: FastTree
		epa: RAxML EPA
		rax: RAxML''')
	parser.add_argument("--kmer", type = int, default = 2, help = "Length of k-mer to use in sequence.")
	parser.add_argument("--histbreaks", type = int, default = 20, help = "Number of histogram breaks for the distribution file.")
	parser.add_argument("-l", "--log_file", type = str, default = 'treeiter.log', help = "File to save the log.")
	parser.add_argument("-w", "--log_level", type = int, default = 3, help = "logging level between 0 (none) and 3")
	parser.add_argument("-k", "--keep_alignments", action = 'store_true', help = "Do not remove individual alignment files")
	parsed_args = parser.parse_args()
	return(parsed_args)


def calculate(alignment, array_file, kmer_num, outgroup):
	logging.info("Begin CALCULATE")
	taxa_stats = alignment.stat_array(kmer = kmer_num, outgroup = outgroup)
	np.set_printoptions(threshold=1000000, linewidth=1000, precision = 4)
	#logging.info("Statistic array: \n\t%s" % "\n\t".join([", ".join([str(y) for y in x]) for x in taxa_stats.tolist()]))
	taxa_stats_list = taxa_stats.tolist()
	logging.info("Statistic array:")
	if len(alignment) == len(taxa_stats_list):
		logging.info("\n" + "\n".join([alignment[index].id + "\t" + "\t".join(str(value) for value in taxa_stats_list[index]) for index in range(len(alignment))]))
	np.savetxt(array_file, taxa_stats)
	logging.info("CALCULATE finalized.")
	return(taxa_stats)

def split_taxa(alignment, statistics_array, outgroup, threshold = 0, hist_file = None, breaks = 20):
	PCA_results = alignment.get_pca(array = statistics_array)
	logging.info("PCA values:")
	logging.info("\n%s" % "\n".join([str(x) + "\t" + str(round(PCA_results[0][x], 5)) for x in PCA_results[0]]))
	logging.info("Proportion of variance explained by PCA %s" % PCA_results[1])
	get_different = alignment.get_different(pca = PCA_results, cutoff = threshold) 
	logging.debug('Two classes identified:')
	logging.debug('First group (to be evaluated individually): %s' % get_different[0])
	logging.debug('Second group (to be used as the base tree): %s' % get_different[1])
	if hist_file is not None:
		logging.info("Saving histogram to file %s" % hist_file)
		logging.disable(logging.CRITICAL) ##TO STOP MATPLOTLIB LOGGING
		hs.hist([PCA_results[0][x] for x in PCA_results[0]], breaks)
		hs.axvline(x = threshold, color = 'b')
		hs.savefig(hist_file)
		logging.disable(logging.NOTSET) ##TO ENABLE LOGGING AGAIN
	if outgroup in get_different[0]:
		raise Exception("Outgroup was identified as an outlier. It is recommended to use a sequence that is not as diverged.")
	with open("problematic_taxa.txt","w") as problematic_taxa:
		for taxon in get_different[0]:
			problematic_taxa.write("%s\n" % taxon)
	return(get_different[0], PCA_results)

def generate_alignments(alignment, directory, prob_taxa):
	count=0
	alignment_list = []
	
	for problematic_taxon in prob_taxa:
		logging.debug('Creating alignment for %s' % problematic_taxon)
		count += 1
		others_index = []
		for other_ones in [x for x in prob_taxa if x != problematic_taxon]:
			others_index.append([x.id for x in alignment].index(other_ones))
		new_alignment = emsa([alignment[x] for x in range(len(alignment)) if x not in others_index])
		namepath = os.getcwd()+'/'+directory+'/'+'iteration_'+str(count)
		new_alignment_path = namepath+'.temp.phy'
		alignment_list.append(new_alignment_path)
		aln_handle = open(new_alignment_path, "w")
		AlignIO.write(new_alignment, aln_handle, "phylip-relaxed")
		aln_handle.close()
	return(alignment_list)


def tree_iterate_all(alignment_list, directory, boots = 100, model = None, data_type = "PROT", method = "fast"):
	#Run tree iteration to create multiple trees.
	#Run a tree search for each alignment using the method specified
	count = 0
	namestr = 'iteration_'+str(count)
	if (method == "rax" and model is None):
		#logging.debug("Default model selected for RAxML")
		if data_type == "PROT":
			evo_model = "PROTCATDAYHOFF"
			#logging.debug("Protein data: PROTCATDAYHOFF selected as model")
		else:
			evo_model = "GTRCAT"
			#logging.debug("Nucleotide data: GTRCAT selected as model")
	for aln in alignment_list:
		count += 1
		namestr = 'iteration_'+str(count)
		logging.debug('Begin run for alignment %s' % aln)
		if method == "fast":
			fasttreecmd = _Fasttree.FastTreeCommandline(input = aln, out = directory+"/ft_"+str(count)+".tre")
			logging.info("FastTree called as: %s" % fasttreecmd)
			fasttreecmd()
			logging.info("FastTree run finalized")
		elif method == "rax":
			seed = random.randint(10000,99999)
			raxml_obj = RaxmlCommandline(sequences = aln, model = evo_model, name=namestr, num_replicates = boots, working_dir = os.getcwd()+'/'+directory+'/', bootstrap_seed = seed)
			logging.info("RAxML called as: %s" % raxml_obj)
			raxml_obj()
			logging.info("RAxML run finalized")
		else:
			raise Exception("Incorrect method specified")
	#Clean up
	if method == "fast":
		iterated_trees = glob.glob(os.getcwd()+'/'+directory+"/ft_*")
	else:
		iterated_trees = glob.glob(os.getcwd()+'/'+directory+"/RAxML_bootstrap."+namestr+"*")
	with open(os.getcwd()+'/resulting_trees.tre','a+') as resulting_file:
		for file in iterated_trees:
			with open(file, 'r') as curr_result:
				resulting_file.write(curr_result.read())
			logging.info('Appended tree(s) from %s into resulting tree file' % file)
			os.remove(file)
	##Clean up RAxML files
	if method == "rax":
		raxml_output = glob.glob(directory+"/RAxML*")
		for file in raxml_output:
			os.remove(file)


def consolidate_trees(treesfile, outgroup):
	out_trees = []
	with open(os.getcwd() + treesfile, 'r') as treefile:
		for line in treefile.readlines():
			out_trees.append(Tree(line))
	logging.info("%i trees loaded." % len(out_trees))
	for i in out_trees:
		try:
			i.set_outgroup(i&outgroup)
		except:
			raise ValueError("Outgroup not found in tree")
	logging.info("Outgroup set on all trees. Begin consensus algorithm.")
	ufc = unequal_full_cons(out_trees)

	#print output
	ufc.ladderize()
	logging.info("Finished. Resulting tree:\n%s" % ufc.write(format=2))
	return(ufc.write(format=2))

def tree_placement(alignment, prob_taxa, model = None, reps = 5, keep = False, temporal_dir = None):
	logging.info("Begin Placement Method")
	if temporal_dir is None:
		temporal_dir = '.elaboo_epa_temp'
		os.mkdir(temporal_dir)
		logging.info('Created temporal directory with name %s' % temporal_dir)
	namepath = os.getcwd()+'/'+temporal_dir+'/'
	base_tree_aln = emsa([x for x in alignment if x.id not in prob_taxa])
	logging.info("Creating alignment excluding problematic taxa.")
	aln_handle = open(namepath+"base_alignment.phy", "w")
	AlignIO.write(base_tree_aln, aln_handle, "phylip-relaxed")
	aln_handle.close()
	#Determine sequence data for default model
	if model is None:
		if alignment.alphabet == "PROT":
			evo_model = "PROTCATDAYHOFF"
		else:
			evo_model = "GTRCAT"
	else:
		evo_model = model
	seed = random.randint(10000,99999)
	raxml_obj = RaxmlCommandline(sequences = namepath+"base_alignment.phy", model = evo_model, name="basetree", num_replicates = reps, working_dir = namepath, parsimony_seed = seed)
	logging.info("RAxML called as: %s" % raxml_obj)
	raxml_obj()
	logging.info("RAxML run finalized")
	if not keep:
		for f in glob.glob(temporal_dir+"/RAxML_[!b]*"):
			os.remove(f)
	full_aln_handle = open(namepath+"clean.phy", "w")
	AlignIO.write(alignment, full_aln_handle, "phylip-relaxed")
	full_aln_handle.close()
	seed = random.randint(10000,99999)
	raxml_epa_obj = RaxmlCommandline(sequences = namepath+"clean.phy", model = evo_model, name="EPA.tre", working_dir = namepath, algorithm = 'v', starting_tree = namepath+"RAxML_bestTree.basetree")
	logging.info("RAxML called as: %s" % raxml_epa_obj)
	raxml_epa_obj()
	logging.debug("RAxML run finalized")
	for f in glob.glob(temporal_dir+"/*"):
		if ("RAxML" in f):
			os.rename(f, f.replace(temporal_dir+'/', "", 1)) #Move files up.
		else:
			os.remove(f)
	os.rmdir(temporal_dir)

