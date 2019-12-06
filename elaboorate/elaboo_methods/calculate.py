import logging
import numpy as np


def calculate(alignment, kmer_num, outgroup, array_file = "elaboo_alignment_stats.txt", taxa_stats_file = None):
	if taxa_stats_file is not None:
		logging.info("Reading from taxa_stats file")
		try:
			taxa_stats = np.genfromtxt(taxa_stats_file)
		except:
			raise Exception("File passed via -s cannot be read correctly.")
		return(taxa_stats)
	logging.info("Begin CALCULATE")
	taxa_stats = alignment.stat_array(kmer = kmer_num, outgroup = outgroup)
	np.set_printoptions(threshold=1000000, linewidth=1000, precision = 4)
	taxa_stats_list = taxa_stats.tolist()
	logging.info("Statistic array:")
	if len(alignment) == len(taxa_stats_list):
		logging.info("\n" + "\n".join([alignment[index].id + "\t" + "\t".join(str(value) for value in taxa_stats_list[index]) for index in range(len(alignment))]))
	np.savetxt(array_file, taxa_stats)
	logging.info("CALCULATE finalized.")
	return(taxa_stats)