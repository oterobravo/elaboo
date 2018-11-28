#Treefuns
#Functions for working with trees
#includes the majority rule plus consensus algorithm
##from Jansson, et al., 2013. https://arxiv.org/pdf/1307.7821.pdf
#Alejandro Otero Bravo
#v0.1 
#Disclaimer: The following is a preliminary version that has not been thoroughly tested.

from Bio import Phylo
from Bio.Phylo.Consensus import *
from io import StringIO
from ete3 import Tree
from collections import Counter


def copy_tree_list(treelist):
	'''copies a list of Tree objects from ete3
	'''
	newlist = []
	for tree in treelist:
		newlist.append(tree.copy())
	return(newlist)

def get_leaves_obj(node, string = False):
	'''returns either a list or a string of the names of the leaves
	sorted for easy comparison. 
	Requires a Tree object from ete3
	'''
	if string == False:
		return(sorted([x.name for x in node]))
	else:
		return(','.join(sorted(x.name for x in node)))

def check_cluster_compatibility(clus1, clus2):
	'''evaluates cluster compatibility based on whether
	C1 is a subset of C2, other way around, or if C1 intersect C2 is 0.
	clusters are list of taxa names
	'''
	leaves1 = set(clus1)
	leaves2 = set(clus2)
	if leaves1 == leaves2:
		return(True)
	elif (leaves1.issubset(leaves2) | leaves2.issubset(leaves1)):
		return(True)
	elif (len(leaves1.intersection(leaves2)) == 0):
		return(True)
	else: 
		return(False)

def check_node_compatibility(node1, node2, exclude = []):
	'''checks if node1 is compatible with node 2 at the exclusion 
	of taxa listed in exclude, if any.
	Requires two Tree objects from ete3 and a list of taxa names in strings
	'''
	clusters1 = []
	for node in node1.traverse():
		clusters1.append(get_leaves_obj(node))
	clusters2 = []
	for node in node2.traverse():
		clusters2.append(get_leaves_obj(node))
	if len(exclude) > 0:
		for item in exclude:
			clusters1 = [c for c in clusters1 if c not in exclude]
			clusters2 = [c for c in clusters2 if c not in exclude]
	compatible = True
	for cluster1 in clusters1:
		for cluster2 in clusters2:
			if not check_cluster_compatibility(cluster1, cluster2):
				compatible = False
				break
	return(compatible)

def majority_rule_plus(treelist):
	'''implements the Majority Rule (+) from Jansson, et al., 2013. 
	Reference at: https://arxiv.org/pdf/1307.7821.pdf
	Requires a list of Tree items from ete3.
	Returns a list including:
	T: tree topology
	N: dictionary with clusters as keys and number of compatible nodes.
	K: dictionary with clusters as keys and number of trees supporting that node.
	Q: dictionary with clusters as keys and number of trees incompatible with that node.
	'''
	t = copy_tree_list(treelist)
	#Phase 1
	#Step 1
	T = t[0].copy()
	#print("Starting tree is:")
	#print(T)
	#Step 2
	N = {}
	#print("First tree is:")
	#print(T)
	for node in T.traverse():
		if node.is_leaf():
			continue
		N.setdefault(get_leaves_obj(node, string = True), 1) 
	
	#Step 3
	for j in range(1, len(t)):
	#Step 3.1
		#print("Starting the evaluation of tree number %s" %j)
		#print(t[j])
		for node in T.traverse():
			if node.is_leaf() | node.is_root():
				continue
			nodeleaves = get_leaves_obj(node)
			#print("Looking for nodes with leaves %s" %nodeleaves)
			matched = False
			for node2 in t[j].traverse():
				if node2.is_leaf() | node2.is_root():
					continue
				#print("node %s from tree T %d has leaves %s" % (node2, j, get_leaves_obj(node2)))
				if nodeleaves == get_leaves_obj(node2):
					#print("Node from T[%d] found" %j)
					#print(node)
					#print(node2)
					#print(T)
					#print(N)
					N[get_leaves_obj(node, string = True)] += 1
					matched = True
					break
			if matched == False:
				#print("No nodes found with leaves %s" %nodeleaves)
				for node2 in t[j].traverse():
					if check_node_compatibility(node, node2) == False:
						#print("Node incompatibility between %s and %s" % (node, node2))
						#print(N)
						N[get_leaves_obj(node, string = T)] += -1
						break
	#Step 3.2
		for key in N:
			if N[key] < 1 & len(key) > 1:
				incompatible_node = T.get_common_ancestor(key.split(","))
				incompatible_node.delete()
			
	
	#Step 3.3
		for node2 in t[j].traverse():
			if node2.is_leaf() | node2.is_root(): break
			node2leaves = get_leaves_obj(node2, string = True)
			for node in T.traverse():
				if node.is_leaf() | node2.is_root(): break
				matched = False
				compatible = False
				if node2leaves == get_leaves_obj(node, string = True):
					matched = True
					#print("%s matched with %s" % (node2, node))
					break
				if check_node_compatibility(node, node2):
					compatible = True
					#print("%s is compatible with %s" % (node2, node))
				if (not matched) & compatible:
					#print("These nodes are compatible and do not match %s %s" % (node2, node))
					#print("leaves are %s" %node2leaves)
					#print(node2)
					try:
						ancestor = T.get_common_ancestor(get_leaves_obj(node2))
					except:
						print("Error, nodes are not connected")
						print("Current tree:")
						return([T, node2])
					#print(ancestor)
					T.prune(set(T.get_leaf_names()).difference(node2.get_leaf_names()))
					ancestor.add_child(node2.copy())
					#add count for this node in the dictionary
					N[','.join([x.name for x in node2])] = 1
		#print(N)
		#print(T)
	#Phase 2
	#Step 4
	K = {}
	Q = {}
	for node in T.traverse():
		if node.is_leaf() | node.is_root() :
			continue
		K.setdefault(get_leaves_obj(node, string = True), 0) 
		Q.setdefault(get_leaves_obj(node, string = True), 0)
	#Step 5
	for j in range(1, len(t)):
		for node in T.traverse():
			if node.is_leaf() | node.is_root() :
				continue
			nodeleaves = get_leaves_obj(node, string = True)
			match = False
			for node2 in t[j].traverse():
				if nodeleaves == get_leaves_obj(node2, string = True):
					K[nodeleaves] += 1
					match = True
					break
			if not match:
			 	compatible = True
			 	for node2 in t[j].traverse():
			 		if not check_node_compatibility(node, node2):
			 			compatible = False
			 			break
			 	if not compatible:
			 		Q[nodeleaves] += 1
	#Step 6
	for node in T.traverse():
		if node.is_leaf() | node.is_root():
			continue
		nodeleaves = get_leaves_obj(node, string = True)
		if K[nodeleaves] <= Q[nodeleaves]:
			try:
				incompatible_node = T.get_common_ancestor(nodeleaves.split(","))
			except:
				print("Error, nodes are not connected")
				#print("Current tree:")
				#print(T)
				#print("current node:")
				#print(node)
				return([T, node])
			incompatible_node.delete()
		
	return([T, N, K, Q])
	
def taxa_set(t):
	'''Requires a list of tree objects
	Returns a list with three sets:
	all unique taxa in all trees
	taxa common to all trees
	and taxa not present in all trees
	'''
	common = set(get_leaves_obj(t[0]))
	all = set(get_leaves_obj(t[0]))
	for trees in t[1:]:
		new_leaves = set(get_leaves_obj(trees))
		common.intersection_update(new_leaves)
		all.update(new_leaves)
	uncommon = all - common
	return([all, common, uncommon])

def unequal_taxa_cons(treelist):
	'''In the presence of unequal taxa sets for all trees, consensus building algorithms
	won't work.
	This creates a consensus tree by ignoring all taxa not present in all trees.
	Uses the majority consensus (+) method.
	'''
	t = copy_tree_list(treelist)
	taxa = taxa_set(t)
	for label in taxa[2]:
		for tree in t:
			try: 
				uncommon_taxon = tree&label
			except:
				continue
			uncommon_taxon.delete()
	return(majority_rule_plus(t))

def find_node_placement(treelist, nodename):
	'''determines the most common node placement of a node named nodename.
	Requires a string 'nodename' with the leaf name and a list of ete3 Tree objects.
	Returns a counter object with a string of taxa names as the keys and the number
	of trees that node appears in.
	'''
	t = copy_tree_list(treelist)
	trees_w_node = []
	for tree in t:
		if nodename in get_leaves_obj(tree):
			trees_w_node.append(tree)
	sisters = []
	for tree in trees_w_node:
		n = tree.search_nodes(name=nodename)[0]
		sisters.append(get_leaves_obj(n.up, string = True))
	sister_count = Counter()
	for word in sisters:
		sister_count[word] += 1
	return(sister_count)

def find_node_placement2(treelist, nodename, exclude  = []):
	def find_parent_excluding(node, exclude = []):
		if node.is_root():
			raise ValueError("%s has no parent that excludes all taxa in question" % node)
		parent = node.up
		parent_leaves = get_leaves_obj(parent)
		if len(exclude) > 1:
			#unexcluded_leaves = [x for x in parent_leaves if x not in exclude.append(nodename)]
			unexcluded_leaves = [x for x in parent_leaves if x not in exclude]
		else:
			#unexcluded_leaves = [x for x in parent_leaves if x not in nodename]
			unexcluded_leaves = parent_leaves
		if len(unexcluded_leaves) > 0:
			return(parent)
		else:
			return(find_parent_excluding(parent, exclude))
	
	t = copy_tree_list(treelist)
	trees_w_node = []
	for tree in t:
		if nodename in get_leaves_obj(tree):
			trees_w_node.append(tree)
	sisters = []
	for tree in trees_w_node:
		n = tree.search_nodes(name = nodename)[0]
		n.sis =  find_parent_excluding(n, exclude = exclude)
		sisters.append(get_leaves_obj(n.sis, string = True))
	sister_count = Counter()
	for word in sisters:
		sister_count[word] += 1
	return(sister_count)
			

def add_most_common_one(nodename, treelist):
	'''This function requires a list of ete3 Tree objects with unequal taxa names.
	By providing a taxa name that is not present in all trees (nodename), it will generate
	a consensus tree with only taxa present in all trees and place taxa nodename on the 
	node that most commonly appears in the trees where it is present.'''
	t = copy_tree_list(treelist)
	place = find_node_placement(t, nodename)
	consensus = unequal_taxa_cons(t)[0]
	try:
		MLnode = place.most_common()[0][0]
	except:
		print("Error")
		return(place)
	MLnode_leaves = MLnode.split(",")
	MLnode_leaves.remove(nodename)
	if len(MLnode_leaves) > 1: #if there is just one sister, it will choose the whole node
		most_common = consensus.get_common_ancestor(MLnode_leaves)
		if len(most_common.get_leaves()) != len(MLnode_leaves):
			raise Exception("sister node is not found in consensus")
		most_common.up.add_child(Tree('new;'))
		new = consensus.search_nodes(name='new')[0]
		new.add_child(most_common.copy())
		most_common.detach()
		new.add_child(Tree('%s;' % nodename))
	elif len(MLnode_leaves) == 1:
		most_common = consensus.search_nodes(name = MLnode_leaves[0])[0]
		sister = most_common.copy()
		most_common.add_child(Tree('%s;' % nodename))
		most_common.add_child(sister)
	return([MLnode_leaves, consensus, most_common])


def add_most_common(nodenames, treelist):
	'''Similar to add_most_common_one. The only difference is this takes a list of node names
	and attempts to find the consensus for all of them, even if taxa are not present in
	the same tree.'''
	t = copy_tree_list(treelist)
	outlist = []
	nodes_to_add = {}
	consensus = unequal_taxa_cons(t)[0]
	outlist.append(consensus.copy())
	for nodename in nodenames:
		place = find_node_placement2(t, nodename, exclude = nodenames)#[x for x in nodenames if x not in nodename])
		try:
			MLnode = place.most_common()[0][0]
		except:
			raise Exception("Error finding most common placement")
		new_node  = ",".join([l for l in MLnode.split(",") if l not in nodenames])
		if len(new_node) < 1:
			print("new_node is %s" % new_node)
			print("nodename is %s" % nodename)
			raise Exception("Two nodes in question are closely related") ###FIX THIS
		nodes_to_add.setdefault(new_node, [])
		nodes_to_add[new_node].append(MLnode.split(","))
	for key in nodes_to_add:
		if len([x for x in key.split(",") if x not in nodenames ]) == 1:
			sibling = consensus.search_nodes(name = key)[0]
		else:
			sibling = consensus.get_common_ancestor(key.split(","))
		try:
			parent = sibling.up
			parent.add_child(name=key+'new')
		except AttributeError:
			print("parent could not be found\nkey currently in play:\n%s" % key)
			break
		except:
			print("unknown error\nkey currently in play: %s" % key)
			break
		new = consensus.search_nodes(name=key+'new')[0]
		free_sibling = sibling.detach()
		new.add_child(free_sibling)
		appendlist = set()
		for newtaxa in nodes_to_add[key]:
			appendlist.update([n for n in newtaxa if n in nodenames])
		for taxa in appendlist:
			new.add_child(name=taxa)
			outlist.append(consensus.copy())
	return([consensus, outlist, nodes_to_add])

def is_node_in(node, tree):
	'''function to find if a certain node is contained in a tree'''
	leaves = sorted([x.name for x in node.get_leaves()])
	for node2 in tree.traverse():
		if leaves == sorted([x.name for x in node2.get_leaves()]):
			return(True)
	return(False)

def is_node_in_excluding(node, tree, exclude):
	'''This will find if a certain node is contained in a tree if a certain number of 
	leaves are ignored or excluded'''
	leaves = sorted([x.name for x in node.get_leaves()])
	for node2 in tree.traverse():
		if [x for x in leaves if x not in exclude] == sorted([x.name for x in node2.get_leaves() if x.name not in exclude]):
			return(True)
	return(False)

def getboot(basetree, trees, percent = True, exclude_nodes = []):
	'''Returns a dictionary with a list of nodes of basetree with values being the number
	of trees that node appears in trees. Percent determines if it is reported as a 
	percentage instead of absolute numbers.'''
	outtree = basetree.copy()
	support = {}
	for node in basetree.traverse():
		leaves = sorted([x.name for x in node.get_leaves()])
		support.setdefault(",".join(leaves), 0)
		for tree in trees:
			if len(exclude_nodes) == 0:
				if is_node_in(node, tree):
					support[",".join(leaves)] += 1
					continue
			else:
				if is_node_in_excluding(node, tree, exclude_nodes):
					support[",".join(leaves)] += 1
					continue
	if percent:
		support.update((k, v / len(trees)) for k,v in support.items())
	return(support)

def mapboot(basetree, trees, percent = True, exclude_nodes = []):
	'''adds bootstrap support values to a tree object based on a list of trees'''
	outtree = basetree.copy()
	sup = getboot(basetree, trees, percent = percent, exclude_nodes = exclude_nodes)
	for key in sup:
		node = outtree.get_common_ancestor(key.split(","))
		node.support = sup[key]
	return(outtree)

def unequal_full_cons(treelist, include_boot = True):
	'''Automatically determines the taxa that are not present in all trees and creates
	a consensus tree based on the common taxa and adds the uncommon ones to the most
	frequent node'''
	taxaset = taxa_set(treelist)
	cons = add_most_common(taxaset[2], treelist)[0]
	if include_boot:
		cons_boot = mapboot(cons, treelist, exclude_nodes = list(taxaset[2]))
		return(cons_boot)
	else:
		return(cons)

	
