#Elaboorate Parser
#Alejandro Otero Bravo

import argparse
import logging

def elaboo_parser():
	parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter)                                               
	parser.add_argument("-a", "--alignment", type = str, required = True, help = "Input alignment")
	parser.add_argument("-f", "--input_format", type = str, default = "fasta", help = "Input alignment format, default: fasta")
	parser.add_argument("-b", "--bootstraps", type = int, default = 10, help = "Number of bootstrap replicates per taxon in question")
	parser.add_argument("-o", "--outgroup", type = str, help = "Outgroup from the alignment")
	parser.add_argument("-s", "--taxa_stats_file", default = None, type=str, help = "File including values for PCA\n(must conform to numpy array standards).")
	parser.add_argument("-p", "--pca_threshold", type=int, default = 0, help = "Threshold for separating groups using the principal component")
	parser.add_argument("-m", "--model", type = str, default = None, help = "Sequence evolution model for tree building.\n(default GTRCAT for nucleotide, PROTCATDAYHOFF for protein)")
	parser.add_argument("-g", "--taxa_list", type = str, default = None, help = "Text file including the list of taxa to be analyzed separately.")
	parser.add_argument("-t", "--trees", type = str, default = None, help = "Newick file of trees to be used for consolidation.")
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
	parser.add_argument("-l", "--log_file", type = str, default = 'elaboo.log', help = "File to save the log.")
	parser.add_argument("-w", "--log_level", type = int, default = 3, help = "logging level between 0 (none) and 3")
	parser.add_argument("-k", "--keep_alignments", action = 'store_true', help = "Do not remove individual alignment files")
	parsed_args = parser.parse_args()
	return(parsed_args)

def start_logging(args):
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