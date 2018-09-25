# TE-invasion-simulations

This program simulates TE invasion into a Drosophila melanogaster like genome, and allows for the evolution of host repression via small RNA mediated silencing. 





Choosing Parameters. Parameters for each simulation are established in the input file (see "Input.txt"). Users can set any of the following parameters below. Parameters are described fully in the associated manuscript.

a. Population size;
b. Basic transposition rate (new insertion per generation per genome per existing copy);
c. Inactivation rate (probability that a copy lose its ability to produce transposase every generation);
d. Max_dysgenesis (maximum fertility cost of TE activity)
e. Number of generations; 
f. Ectopic recombination rate (a modifier to calculate the ectopic recombination cost to fitness);
g. Snapshot (generations in the simulation for which you'd like full genotypes reported for every individual)




Establishing fitness effects. The composition of the genome is described in the associated manuscript. Briefly, each genome is comprised of pisites (p), pseudo-pi sites (q), neutral sites (n), and functional sites (a-j). The input file contains 3 chromosomes (X,2,3) with a single letter indicating the class for each site. For functional sites, the following fitness values are assigned to homozygotes:

b = 0
d = 0.25
f = 0.05
h = 0.99
j = 1.01

To implement dominant fitness effects, letters a,c,e,g,i can be used. Their fitness values must be assigned directly in the simulator.pl code. Homozygous effects of classes b,d,f,h,j can also be altered in the simulator.pl

Usage.

./simulator.pl Input.txt <output.file> <chrom> <site>
  
  where, chromosome and site indicate the chromosomes (0,1,2 correspond to X, 2nd and 3rd chromosomes repsctively) and the particular site that is occupied by an active TE in 5% of the population at the start of the simulation.
  
Output.

For each run, 3 output files are generated:

1) Snapshot files. Provide full genotypes for every individual in the population are reported, for every generation for which a snap shot was requested in the input file. Unoccupied sites are indicated as 0, occupied sites with active TEs are indicated as 1, occupied sites with inactive TEs are indicated as 2.  Each indvidual is reported as follows:
Individual # 0..n-1
Genotype of chromosome 1, copy 1
Genotype of chromosome 1, copy 2
Genotype of chromosome 2, copy 1
Genotype of chromosome 2, copy 2
Genotype of chromosome 3, copy 1
Genotype of chromosome 3, copy 2
1 1 <number of occupied small RNA sites in maternal genome> <number of occupied non-small RNA sites in maternal genome> 
  
2) Genotype files contain the consolidate genotypes for the full population every 20 generations across the simulation. For each generation the following are reported
Generation #
X chromosome number
autosome number
A tab deliminted table indicating:
<chromosome> <site> <# of chromosomes with active TEs> <# of chromosomes with inactive TEs>
For every site that is occupied in at least one individual. Sites not occupied in any individuals are omitted.
  
3) Population stats files reports summary statistics for every generation.
Number of TE: indicates mean and standard deviation of individual TE copy number
Proportion infected: proportion of the population bearing one more more TE copies
Proportion able to transpose: proportion of the population bearing at least one active TE
Number of piTE: indicates mean and standard deviation of insertions in small RNA sites
Number of pseudo-piTE: indicates mean and standard deviation of insertions in pseudo small RNA sites
Total fitness: mean and standard deviation of individual fitness
Transposition rate: mean and standard deviation in transposition rate per individual genome
Maternal repression: mean and standard deviation in maternal small RNA mediated repression



