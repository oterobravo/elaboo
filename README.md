# ELABOORATE
Extract & Leave-All-But-One-Out Reconstruction. 

For taxa with high mutation rates, LBA (long branch attraction) artifacts can generate erroneous relationships. Elaboorate attempts to identify and extract taxa that are likely to be erroneously placed, and use a leave-all-but-one-out method to indicate the correct placement of each taxa. 

## Getting Started

### Requirements  
ELABOORATE requires python 3+ and the following extra packages (maybe less):

```  
biopython==1.70
cycler==0.10.0
ete3==3.1.1
kiwisolver==1.0.1
matplotlib==2.2.2
numpy==1.14.1
p4===1.2.0.-2016-05-03-
pyparsing==2.2.0
python-dateutil==2.7.2
pytz==2018.4
scikit-learn==0.19.1
scipy==1.0.1
six==1.11.0
sklearn==0.0
```  
Additionally, the scripts emsa.py, treefuns.py and iter_funcs.py must be in the same folder or in the library path. 


## Quick Start

For simplicity, you can use a virtual environment.
```
python3 -m venv elaboo
cd elaboo/
source bin/activate
pip install -r requirements.txt
```

The easiest way to run Elaboorate is to use all default parameters. Simply attach an alignment and indicate an outgroup.

```
python elaboorate.py -s alignment.fas -o outgroup_name
```

For a complete list of options use the help option `-h`.

## Modules
Elaboo contains 4 modules that perform the different steps of the process. They can be used together or separately to evaluate intermediate steps of the process or to use external files to proceed. 

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
