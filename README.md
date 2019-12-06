# ELABOORATE
Extract & Leave-All-But-One-Out Reconstruction. 

For taxa with high mutation rates, LBA (long branch attraction) artifacts can generate erroneous relationships. Elaboorate attempts to identify and extract taxa that are likely to be erroneously placed, and use a leave-all-but-one-out method to indicate the correct placement of each taxa. 

## Getting Started

### Requirements  
ELABOORATE requires python 3+ and the packages listed in the requirements.txt file.  
```
pip install -r requirements
```

## Quick Start

For simplicity, you can use a virtual environment. If using `venv`:
```
python3 -m venv elaboo
cd elaboo/
source bin/activate
pip install -r requirements.txt
```

(Also for simplicity,do `alias elaboo='/path/to/elaboo.py'`. )

The minimal requirements for elaboorate are an alignment and the name of the outgroup to use. This will run the full pipeline with the default parameters.  
```
python elaboo.py -s alignment.fas -o outgroup_name
```
For a complete list of options and the default parameters use the help option `-h`.

To run individual modules, see the Algorithms section below.

## Modules
Elaboo contains 4 modules that perform the different steps of the process. They can be used together or separately to evaluate intermediate steps of the process or to use external files to proceed. 
It's recommended to run each step individually with a new dataset to see what's happening and ensure the taxa selected by SPLIT make sense. Check `elaboo.log` in between runs.
```
elaboo -s alignment.fas -o outgroup_name -y b 
```
Creates `elaboo_alignment_stats.txt`.

```
elaboo -s alignment.fas -o outgroup_name -y h -s elaboo_alignment_stats.txt 
```
Generates principal component and splits taxa, list on `elaboo_problematic_taxa.txt`.
Recommended: `-d histogram.png` to view the distribution of parameters and cutoff (requires matplotlib).

```
elaboo -s alignment.fas -o outgroup_name -y i -g elaboo_problematic_taxa.txt -fast 
```
Iterate tree reconstruction. Produce all the trees and save in `elaboo_resulting_trees.txt`.

```
elaboo -s alignment.fas -o outgroup_name -y g -t elaboo_resulting_trees.txt 
```
Creates consensus tree from resulting trees.


The four modules are:

##### CALCULATE
Generates summary statistics for each taxa to indicate possible problematic taxa. 
The default summary statistics used are distance from an outgroup, base composition, and 2-kmer frequency.

Required input: alignment file (`-a`), outgroup(`-o`)
Suggested input: distribution (`-d file`) to create a figure file with the distribution of the first principal component. This will help in the next step.
Output: taxa stats file (text file with statistic matrix).

##### SPLIT
Uses PCA on the summary statistics to identify potentially problematic taxa.

Required input: alignment file (`-a`), outgroup(`-o`), taxa stats file (`-s`).
Suggested input: PCA threshold (`-p`) depending on your data.
Output: taxa list file (text file with list of taxa to be analyzed separately).

##### BUILD TREES
Creates individual alignments and infers trees from each alignment using RAxML.

Required input: alignment file (`-a`), outgroup(`-o`), taxa list file (`-g`).
Suggested input: EPA (`-e`) placement algorithm for speed in exploratory runs.
Output: tree file (list of trees from individual reconstructions).

##### CONSOLIDATE
Generates a consensus tree from the trees built.

Required input: alignment file (`-a`), outgroup(`-o`), tree files (`t`).
Output: tree file with consensus tree. 

## Algorithms
Elaboorate allows every module to be called individually and to use external data at every step of the way. The different algorithms control which modules are used as shown here:

|mode|CALCULATE|SPLIT|BUILD_TREES|CONSOLIDATE|
|:---:|:---:|:---:|:---:|:---:|
|a|*|*|*|*|
|b|*||||
|c|*|*|||
|d|*|*|*||
|e||*|*|*|
|f|||*|*|
|g||||*|
|h||*|||
|i|||*||
