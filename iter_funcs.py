#Functions for the iteration of the algorithm
#Alejandro Otero Bravo
#v0.1
#Modules are defined here.

#Disclaimer: The following is a preliminary version that has not been thoroughly tested.

from treefuns import *
from emsa import *
from Bio.Phylo.Applications import RaxmlCommandline
import glob
import os
import random
import string 
import logging
import matplotlib.pyplot as hs

def calculate(alignment, array_file, kmer_num, outgroup):
	taxa_stats = alignment.stat_array(kmer = kmer_num, outgroup = outgroup)
	np.set_printoptions(threshold=1000000, linewidth=1000, precision = 4)
	logging.info("Statistic array: \n\t%s" % "\n\t".join([", ".join([str(y) for y in x]) for x in taxa_stats.tolist()]))
	np.savetxt(array_file, taxa_stats)
	return(taxa_stats)

def split_taxa(alignment, statistics_array, outgroup, threshold = 0, hist_file = None, breaks = 20):
	PCA_results = alignment.get_pca(array = statistics_array)
	logging.info("PCA values:\n%s" % "\n\t".join([str(round(PCA_results[0][x], 5))+": "+str(x) for x in PCA_results[0]]))
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
		alignment_list.append(namepath+'.phy')
		aln_handle = open(namepath+'.phy', "w")
		AlignIO.write(new_alignment, aln_handle, "phylip-relaxed")
		aln_handle.close()
	return(alignment_list)

def tree_iterate(alignment_list, model, boots, directory, data_type):
	#Run tree iteration to create multiple trees. For each problematic taxon, create an alignment
	#excluding all other problematic taxa and run RAxML.
	count = 0
	for aln in alignment_list:
		count += 1
		logging.debug('Begin RAxML run for alignment %s' % aln)
		seed = random.randint(10000,99999)
		namestr = 'iteration_'+str(count)
		#Determine sequence data for default model
		if model is None:
			if data_type == "PROT":
				evo_model = "PROTCATDAYHOFF"
			else:
				evo_model = "GTRCAT"
		else:
			evo_model = model		
		raxml_obj = RaxmlCommandline(sequences = aln, model = evo_model, name=namestr, num_replicates = boots, working_dir = os.getcwd()+'/'+directory+'/', bootstrap_seed = seed)
		logging.info("RAxML called as: %s" % raxml_obj)
		raxml_obj()
		logging.info("RAxML run finalized")
		
		#Save desired output and delete the rest
		raxml_trees = glob.glob(os.getcwd()+'/'+directory+"/RAxML_bootstrap."+namestr+"*")
		with open(os.getcwd()+'/resulting_trees.tre','a+') as resulting_file:
			for file in raxml_trees:
				with open(file, 'r') as curr_result:
					resulting_file.write(curr_result.read())
				logging.info('Appended tree(s) from %s into resulting tree file' % file)
				os.remove(file)
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
	for f in glob.glob(temporal_dir+"/RAxML_*EPA*"):
		os.rename(f, f.replace(temporal_dir+'/', "", 1)) #Move files up.
	os.rmdir(temporal_dir)
	logging.info("RAxML run finalized")
