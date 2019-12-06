import logging

def draw_histogram(values, hist_file, threshold, breaks):
	try: 
		logging.getLogger('matplotlib').setLevel(logging.ERROR) #Disable matplotlib logging
		import matplotlib.pyplot as hs
	except ImportError:
		logging.info("Error importing pyplot. Cannot draw histogram.")
		return
	logging.info("Saving histogram to file %s" % hist_file)
	hs.hist(values, breaks)
	hs.axvline(x = threshold, color = 'b')
	hs.savefig(hist_file)
	return

def get_PCAs(alignment, statistics_array):
	PCA_results = alignment.get_pca(array = statistics_array)
	logging.info("PCA values:")
	logging.info("\n%s" % "\n".join([str(x) + "\t" + str(round(PCA_results[0][x], 5)) for x in PCA_results[0]]))
	logging.info("Proportion of variance explained by PCA %s" % PCA_results[1])
	return(PCA_results)

def split_taxa(alignment, outgroup, threshold = 0, hist_file = None, breaks = 20, PCA_results = None):
	logging.info("Begin SPLIT")
	get_different = alignment.get_different(pca = PCA_results, cutoff = threshold) 
	logging.debug('Two classes identified:')
	logging.debug('First group (to be evaluated individually): %s' % get_different[0])
	logging.debug('Second group (to be used as the base tree): %s' % get_different[1])
	if hist_file is not None:
		draw_histogram([PCA_results[0][x] for x in PCA_results[0]], hist_file, threshold, breaks)
	if outgroup in get_different[0]:
		raise Exception("Outgroup %s was identified as an outlier. It is recommended to use a sequence that is not as diverged." % outgroup)
	with open("elaboo_problematic_taxa.txt", "w") as problematic_taxa:
		for taxon in get_different[0]:
			problematic_taxa.write("%s\n" % taxon)
	logging.info("SPLIT finalized.")
	return(get_different[0])