import logging
import os
import treefunctions as tf
from ete3 import Tree

def consolidate_trees(treesfile, outgroup):
	out_trees = []
	with open(treesfile, 'r') as treefile:
		for line in treefile.readlines():
			out_trees.append(Tree(line))
	logging.info("%i trees loaded." % len(out_trees))
	for i in out_trees:
		try:
			i.set_outgroup(i&outgroup)
		except:
			raise ValueError("Outgroup not found in tree")
	logging.info("Outgroup set on all trees. Begin consensus algorithm.")
	ufc = tf.unequal_full_cons(out_trees)

	#print output
	ufc.ladderize()
	logging.info("Finished. Resulting tree:\n%s" % ufc.write(format=2))
	return(ufc.write(format=2))