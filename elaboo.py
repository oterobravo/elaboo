# ELABOORATE
#Extract & Leave-All-But-One-Out Reconstruction.
#Alejandro Otero Bravo
#v0.3
#Disclaimer: The following is a preliminary version that has not been thoroughly tested.

import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)                                               

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
args = parser.parse_args()


import glob
import os
import random
import string 
import logging
from iter_funcs import *

if args.log_level == 3:
	log_level = logging.DEBUG
elif args.log_level == 2:
	log_level = logging.INFO
elif args.log_level == 1:
	log_level = logging.WARNING
elif args.log_level == 0:
	log_level = logging.CRITICAL
else:
	raise ValueError("log level (-w) is invalid.")

logging.getLogger(__name__)
logging.basicConfig(filename = args.log_file, filemode = 'w', level = log_level, format = '%(asctime)s - %(funcName)s - %(message)s', datefmt='%I:%M:%S')
logging.info('Log started')
logging.info('Elaboo called with parameters: %s' % " ".join(["\n\t"+x+": "+str(getattr(args, x)) for x in vars(args)]))


CALCULATE = True if (args.algorithm in "abcd") else False
SPLIT = True if (args.algorithm in "acdeh") else False
TREE_BUILD = True if (args.algorithm in "adefi") else False
CONSOLIDATE = True if (args.algorithm in "aefg") else False

logging.info("Algorithm selected: %s\n\tCalculate: %s\n\tSplit: %s\n\tTree Building: %s\n\tConsolidate: %s" % (args.algorithm, CALCULATE, SPLIT, TREE_BUILD, CONSOLIDATE))
##

#Always require alignment
input_alignment = AlignIO.read(open(args.alignment), args.input_format)
logging.info("input file %s successfully read" % args.alignment)
logging.info("Taxa: %i, length: %i" % (len(input_alignment), len(input_alignment[0])))

alignment_prepared = emsa(input_alignment._records)
alignment_prepared.Alphabet()
if args.clean:
	alignment_prepared = alignment_prepared.clean()

logging.debug('Alignment transformed and ambiguities have been stripped.')
logging.info("Taxa: %i, length: %i" % (len(alignment_prepared), len(alignment_prepared[0])))

datatype = alignment_prepared.Alphabet()

outgroup = args.outgroup
if outgroup is None:
	outgroup = alignment_prepared[0].id
	logging.info('Output not provided, %s was chosen.' % outgroup)
else:
	logging.info('Output provided: %s' % outgroup)
	if outgroup not in [x.id for x in alignment_prepared]:
		raise ValueError('Output provided not found in input alignment.')


if (CALCULATE & (args.taxa_stats_file is None)):
	logging.info("Begin CALCULATE")
	taxa_stats = calculate(alignment_prepared, array_file = "alignment_stats.txt", kmer_num = args.kmer, outgroup = alignment_prepared.is_named(outgroup))
	logging.info("CALCULATE finalized.")

if SPLIT:
	logging.info("Begin SPLIT")
	if args.taxa_stats_file is not None:
		logging.info("Reading from taxa_stats file")
		try:
			taxa_stats = np.genfromtxt(args.taxa_stats_file)
		except:
			raise Exception("Statistics for each taxon were not calculated and no file was given. Add a file using -s or set CALCULATE to True.")
	problem_taxa, PCA_results = split_taxa(alignment_prepared, taxa_stats, outgroup, threshold = args.pca_threshold, hist_file = args.distribution, breaks = args.histbreaks )
	logging.info("SPLIT finalized.")

if TREE_BUILD:
	logging.info("Begin TREE_BUILD")
	temporal_dir = '.elaboo_temp'
	os.mkdir(temporal_dir)
	logging.info('Created temporal directory with name %s' % temporal_dir)
	try:
		problem_taxa
	except NameError:
		try:
			problem_taxa = []
			with open(args.taxa_list, 'r') as list_from_file:
				for line in list_from_file: 
					problem_taxa.append(line.strip())
		except:
			raise Exception("Taxa to be used have not been calculated nor provided. Switch mode or provide taxa_list.")
		else:
			logging.info("Taxa read from file %s" % args.taxa_list)
	if args.tree_build == "epa":
		CONSOLIDATE = False
		tree_placement(alignment_prepared, problem_taxa, model = args.model, temporal_dir = temporal_dir, keep = args.keep_alignments)
		logging.info("EPA Finalized.")
		#os.rmdir(temporal_dir)
	elif args.tree_build == "rax":
		logging.info('Generating alignments leaving all but one taxa out.')
		laboo_alignments = generate_alignments(alignment_prepared, temporal_dir, prob_taxa = problem_taxa)
		logging.info('Iterating tree reconstruction using RAxML')
		tree_iterate_all(laboo_alignments, model = args.model, boots = args.bootstraps, data_type = datatype, directory = temporal_dir, method = "rax")
		#Remove intermediate files, if desired
		if args.keep_alignments == False:
			for f in glob.glob(temporal_dir+"/*"):
				os.remove(f)
				logging.debug('Removed %s.' % f)
			logging.info('Temporary alignments removed.')
		else:
			for f in glob.glob(temporal_dir+"/*"):
				os.rename(f, f.replace(temporal_dir+"/", "", 1))
				logging.debug("Saved %s." % f)
		os.rmdir(temporal_dir)
	elif args.tree_build == "fast":
		logging.info('Generating alignments leaving all but one taxa out.')
		laboo_alignments = generate_alignments(alignment_prepared, temporal_dir, prob_taxa = problem_taxa)
		logging.info('Iterating tree reconstruction using FastTree')
		tree_iterate_all(laboo_alignments, temporal_dir, method = "fast")
		logging.info("Tree reconstruction Finalized.")
		if args.keep_alignments == False:
			for f in glob.glob(temporal_dir+"/*"):
				os.remove(f)
			logging.info('Temporary alignments removed.')
		else:
			for f in glob.glob(temporal_dir+"/*"):
				os.rename(f, f.replace(temporal_dir+"/", "", 1))
		os.rmdir(temporal_dir)
		##Troubleshoot this and remove redundancy
	else:
		raise Exception("Invalid Tree Building Method. Use 'rax', 'epa' or 'fast'.")
	logging.info("TREE_BUILD finalized.")

if CONSOLIDATE:
	try:
		final_tree = consolidate_trees('/'+args.trees, outgroup)
	except IOError:
		raise Exception('Problem opening file of trees for consolidation')
	logging.info("Final topology:\n %s" % final_tree)
	print(final_tree)
