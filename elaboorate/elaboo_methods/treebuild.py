import logging
import tempfile
import os
import glob
import random
import emsa #Replace with Bio.MultipleSeqAlignment later
from Bio.Phylo import Applications
from Bio import Align, AlignIO

def check_for_problem_taxa(problem_taxa, taxa_list):
	#requires logging
	if problem_taxa is None and args.taxa_list is None:
		raise Exception("Taxa to be used have not been calculated nor provided. Switch mode or provide taxa_list via -g.")
	elif problem_taxa is None:
		problem_taxa = []
		with open(args.taxa_list, 'r') as list_from_file:
			for line in list_from_file: 
				problem_taxa.append(line.strip())
		logging.info("Taxa read from file %s" % args.taxa_list)
			
	return(problem_taxa)

def treebuilder(alignment, problem_taxa, args):
	#Requires tempfile, logging, 
	logging.info("Begin TREE_BUILD")
	with tempfile.TemporaryDirectory() as temporal_dir:
		logging.info('Created temporal directory with name %s' % temporal_dir)
		outdir = os.getcwd()+"/elaboo_"+args.tree_build+"_results"
		try:
			os.mkdir(outdir)
		except:
			raise Exception("Could not create output directory. Please remove results from previous run.")
		if args.tree_build == "fast":
			result_file = fast_method(alignment, temporal_dir, problem_taxa, outdir)####
		else:
			if args.model is None:
				if alignment.alphabet == "PROT":
					evo_model = "PROTCATDAYHOFF"
				else:
					evo_model = "GTRCAT"
			else:
				evo_model = args.model
			if args.tree_build == "epa":
				result_file = epa_method(alignment, problem_taxa, model = evo_model, keep = args.keep_alignments, temporal_dir = temporal_dir) ###
			elif args.tree_build == "rax":
				result_file = rax_method(alignment, problem_taxa, evo_model, temporal_dir, boots = args.bootstraps, outdir = outdir) #####
			else:
				raise Exception("Invalid Tree Building Method. Use 'rax', 'epa' or 'fast'.")
			for f in glob.glob(temporal_dir+"/*"):
				if ("RAxML" in f):
					logging.info("Moving %s to %s" % (f, f.replace(temporal_dir+'/', outdir+"/", 1)))
					os.rename(f, f.replace(temporal_dir+'/', outdir+"/", 1)) #Save files
				else:
					logging.info("Removing %s" % f)
					os.remove(f)
		return(result_file)
			

def generate_alignments(alignment, directory, prob_taxa):
	count=0
	alignment_list = []
	
	for problematic_taxon in prob_taxa:
		logging.debug('Creating alignment for %s' % problematic_taxon)
		count += 1
		others_index = []
		for other_ones in [x for x in prob_taxa if x != problematic_taxon]:
			others_index.append([x.id for x in alignment].index(other_ones))
		new_alignment = emsa.Emsa([alignment[x] for x in range(len(alignment)) if x not in others_index])
		namepath = directory+'/'+'iteration_'+str(count)
		new_alignment_path = namepath+'.temp.phy'
		alignment_list.append(new_alignment_path)
		with open(new_alignment_path, "w") as aln_handle:
			AlignIO.write(new_alignment, aln_handle, "phylip-relaxed")

	return(alignment_list)

def epa_method(alignment, prob_taxa, model = None, reps = 5, keep = False, temporal_dir = None):
	#Requires os, emsa, Bio.Phylo.Applications, AlignIO, glob
	logging.info("Begin Placement Method")
	namepath = temporal_dir+'/'
	base_tree_aln = emsa.Emsa([x for x in alignment if x.id not in prob_taxa])
	
	logging.info("Creating alignment excluding problematic taxa.")
	with open(namepath+"base_alignment.phy", "w") as aln_handle:
		AlignIO.write(base_tree_aln, aln_handle, "phylip-relaxed")
	
	#Base tree without problem taxa
	seed = random.randint(10000,99999)
	raxml_obj = Applications.RaxmlCommandline(sequences = namepath+"base_alignment.phy", model = model, name="basetree", num_replicates = reps, working_dir = temporal_dir, parsimony_seed = seed)
	logging.info("RAxML called as: %s" % raxml_obj)
	raxml_obj()

	with open(namepath+"total.phy", "w") as full_aln_handle:
		AlignIO.write(alignment, full_aln_handle, "phylip-relaxed")
	
	#Run EPA
	seed = random.randint(10000,99999)
	raxml_epa_obj = Applications.RaxmlCommandline(sequences = namepath+"total.phy", model = model, name="EPA.tre", working_dir = namepath, algorithm = 'v', starting_tree = namepath+"RAxML_bestTree.basetree")
	logging.info("RAxML called as: %s" % raxml_epa_obj)
	raxml_epa_obj()
	logging.debug("RAxML-EPA run finalized")
	return(None)

def fast_method(alignment, tempdir, problem_taxa, outdir):
	logging.info('Generating alignments leaving all but one taxa out.')
	alignment_list = generate_alignments(alignment, tempdir, prob_taxa = problem_taxa)
	logging.info('Iterating tree reconstruction using FastTree')
	count = 0
	namestr = 'iteration_'+str(count)
	for aln in alignment_list:
		count += 1
		namestr = 'iteration_'+str(count)
		logging.debug('Begin run for alignment %s' % aln)
		fasttreecmd = Applications._Fasttree.FastTreeCommandline(input = aln, out = outdir+"/ft_"+str(count)+".tre")
		logging.info("FastTree called as: %s" % fasttreecmd)
		fasttreecmd()
		logging.info("FastTree run finalized")
	logging.info("Tree reconstruction Finalized.")
	iterated_trees = glob.glob(outdir+"/ft_*.tre")
	result_file = outdir+"/elaboo_resulting_trees.tre"
	with open(result_file, 'a+') as resulting_file:
		for file in iterated_trees:
			with open(file, 'r') as curr_result:
				resulting_file.write(curr_result.read())
			logging.info('Appended tree(s) from %s into resulting tree file' % file)
	return(result_file)

def rax_method(alignment, problem_taxa, model, directory = ".", boots = 10, outdir = "."):
	logging.info('Generating alignments leaving all but one taxa out.')
	alignment_list = generate_alignments(alignment, directory, prob_taxa = problem_taxa)
	logging.info('Iterating tree reconstruction using RAxML')
	
	count = 0
	namestr = 'iteration_'+str(count)
	for aln in alignment_list:
		count += 1
		namestr = 'iteration_'+str(count)
		logging.debug('Begin run for alignment %s' % aln)
		seed = random.randint(10000,99999)
		raxml_obj = Applications.RaxmlCommandline(sequences = aln, model = model, name=namestr, num_replicates = boots, working_dir = directory+'/', bootstrap_seed = seed)
		logging.info("RAxML called as: %s" % raxml_obj)
		raxml_obj()
		logging.info("RAxML run finalized")

	#Clean up
	iterated_trees = glob.glob(directory+"/RAxML_bootstrap."+namestr+"*")
	result_file = outdir+"/elaboo_resulting_trees.tre"
	with open(result_file, 'a+') as resulting_file:
		for file in iterated_trees:
			with open(file, 'r') as curr_result:
				resulting_file.write(curr_result.read())
			logging.info('Appended tree(s) from %s into resulting tree file' % file)
	logging.info("Iteration of RAxML tree construction finished.")
	return(result_file)
