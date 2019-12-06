#!/usr/bin/env python3
#ELABOORATE
#Extract & Leave-All-But-One-Out Reconstruction.
#Alejandro Otero Bravo
#v0.5

import glob
import os
import logging
import emsa 
from elaboo_parser import elaboo_parser, start_logging
import elaboo_methods as em


def main(args):
	start_logging(args)
	
	CALCULATE = True if (args.algorithm in "abcd") else False
	SPLIT = True if (args.algorithm in "acdeh") else False
	TREE_BUILD = True if (args.algorithm in "adefi") else False
	CONSOLIDATE = True if (args.algorithm in "aefg" and args.tree_build != "epa") else False

	logging.info("Algorithm selected: %s\n\tCalculate: %s\n\tSplit: %s\n\tTree Building: %s\n\tConsolidate: %s" % (args.algorithm, CALCULATE, SPLIT, TREE_BUILD, CONSOLIDATE))
	##
	#Always require alignment
	alignment_prepared = emsa.read_alignment(args.alignment, args.input_format, args.clean)
	datatype = alignment_prepared.Alphabet()

	outgroup = args.outgroup
	if outgroup is None:
		outgroup = alignment_prepared[0].id
		logging.info('Output not provided, %s was chosen.' % outgroup)
	else:
		logging.info('Output provided: %s' % outgroup)
		if outgroup not in [x.id for x in alignment_prepared]:
			raise ValueError('Output provided not found in input alignment.')

	if CALCULATE:
		taxa_stats = em.calculate(alignment_prepared, kmer_num = args.kmer, outgroup = alignment_prepared.is_named(outgroup), taxa_stats_file = args.taxa_stats_file)

	if SPLIT:
		PCA_results = em.get_PCAs(alignment_prepared, taxa_stats)
		problem_taxa = em.split_taxa(alignment_prepared, outgroup, threshold = args.pca_threshold, hist_file = args.distribution, breaks = args.histbreaks, PCA_results = PCA_results )
	else:
		problem_taxa = None
	
	if TREE_BUILD:
		problem_taxa = em.check_for_problem_taxa(problem_taxa, args.taxa_list)
		result_trees = em.treebuilder(alignment_prepared, problem_taxa, args)
	else:
		result_trees = None

	if CONSOLIDATE:
		if result_trees is None and args.trees is None:
			raise Exception("No trees to consolidate. Use TREE_BUILD or provide a tree file via -t")
		elif result_trees is None:
			result_trees = "/"+args.trees
		try:
			final_tree = em.consolidate_trees(result_trees, outgroup)
		except IOError:
			raise Exception('Problem opening file of trees for consolidation')
		logging.info("Final topology:\n %s" % final_tree)
	return(1)

if __name__ == "__main__":
	parsed_args = elaboo_parser()
	main(parsed_args)
