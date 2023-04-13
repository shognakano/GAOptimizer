# GAOptimizer
GAOptimizer is a script that uses a genetic algorithm to optimize the combination of consensus mutations that stabilize a protein structure. 
It introduces a consensus mutation into the template protein structure (PDB) by random mutation (specified by mutation_rate, default is 30%) or by recombination mutation (70%).
Total 100 of children would be generated per one generation, and selection of the children would be performed as tournament selection based on the scores implemented in PyRosetta.

[Required software]
Python 3.7.6 or higher, PyRosetta-4, Biopython 1.77, MAFFT (should be able to run from command line)

[Required files]
3D structure of the template protein to introduce the mutation (PDB data)
Sequence data of the template protein homologs used to identify consensus mutations (Fasta format, up to 50 sequences are recommended)


