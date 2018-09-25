#!/usr/bin/perl -w
use strict;


#Parameters that remain constant throughout the simulation

my $popsize; #Population Size
my @chroms; #Size and number of chromosomes CHROM NO.1 IS SEX CHROMOSOME!
my $gsize_male; #Diploid genome size of male
my $gsize_female; #Diploid genome size of female
my @recom_rate; #Recombination rate
my $ect_basic; #Basic ectopic recombination rate
my $u_basic; #Basic "u" for ONE existing copy of TE. u is number of new insertions.
my $x; #Probability of a TE become dead
my @pi; #Locations of PiRNA regions
my @func; #Locations of functional regions
my $numchr; #Number of chromosomes
my $gener_tot; #Total Generations
my @u_gene_list; #List of u-deciding gene alleles
my @morgan; #The total mean recombination number for each chromosome
my $outfile; #Name of output file
my $evo_gens; #Generations evolved in the whole program
my $chr_start; #Chromosome of starting TE site
my $site_start; #Coord of starting TE site
my $model; #Selection/repression model
my $max_dys; #maximum fitness cost imposed by hybrid dysgenesis

#Parameters pertaining to continued evolutionn
my $pickup; #Whether this run is picking up from another
my $pickfile; #File name of picked up file
my $pickup_lines; #Lines of each individual in the file
my @startings; #File read of picked up file
my @snapshot; #Number of generations when snapshots are taken

#Parameters that change each generation
my @geno; #Genotype of each individual
my @sex; #Sex of each individual; 0-M, 1-F
my @fit; #Fitness of each individual
my @fer; #Fertility of each individual (only consider the sterility caused by high-TE/low-TE mating)
my @fit_tot; #Total fitness; product of the two values above
my @nte; #Number of TE in each individual
my @nte_inpi; #Number of TE in each individual in PiRNA regions
my @nte_inpiq; #Number of TE in each individual in PiRNA and pseudo-PiRNA regions
my @func_disr; #Factor of fitness due to function disruption
my @nte_act; #Whether there are at least one active copy
my @u_gene1;
my @u_gene2; #1st and 2nd copy of u-deciding gene
my @u_gene_freq; #Frequency of u-deciding gene genotypes
my @u; #Transposition rate of each individual, u is number of new insertions.
my $gener; #Generation
my @mum_pi; #Number of TE in mother in pi
my @mum_nonpi; #Number of TE in mother not in pi
my @rep; #strength of piRNA mediated repression for each individuals
my @dysg; #Whether there is dysgenesis effect that negate pi effect
my @telist; #List of all TEs, used for selecting parent TE, format: "3 10" is chromosome 2 copy 1 site 10.
my @telist_temp; #List of all TEs, used for selecting parent TE

my @geno_new; #Genotype of the progeny
my @sex_new; #Sex of the progeny
my @mum_pi_new;
my @mum_nonpi_new;
my @u_gene1_new;
my @u_gene2_new; #1st and 2nd copy of u-deciding gene

my $sum_fit_male;
my $sum_fit_female;

my $total_nte;

#Placeholder, these change usually each individual
my $indi;
my $chr;
my $new;
my $mating_m;
my $mating_f;
my @sperm;
my @egg;
my $allele;
my $random; my $random1; my $random2;

#Output
my @result; #population stats
my @result2; #genotype stats
my @result3; #snapshot/final full genotype
my @props;
#die;#Just to test for syntax errors


#Add an independent gene controlling base u.
#Add the "really dead" TEs that do not jump.

#Initialization
inputing();
#die; #Check input

#Setting up population
if ($pickup == 1) {#Picking up a previous run
 open (PICKUP, $pickfile);
 @startings = <PICKUP>; foreach (@startings) {chomp;}
 #make a vector where each line in the snapshot file corresponds to a specific value
 $pickup_lines = ($#chroms+2)*2;
 #the number of lines per individual should be the last index of the chromosome array which is defined by the index file (2, if there are 3 chromosomes) + 2 (for the chromosome with index zero and the number and maternal genotype of each individual), multiplied by 2. With one pair of sex chromosomes and two autosomes, 8 lines are expected.
 if (($#startings+1)/$pickup_lines != $popsize) {die "Pickup File Size Or Format Error\n";}
 #the population size, as indicated in the input file, should be equal number of lines in the pickup file divided by the expected number of lines per individual 
 foreach $indi (0..($popsize-1)) {
 #loop through each individual in the population according to their number
  unless ($startings[$indi*$pickup_lines] =~ /(\d+)/) {die "Pickup File Size Or Format Error\n";}
  #The 8th line of each individuals genotype information should be comprised of digits, If it isn't, there is a file format error.
  if ($1 != $indi) {die "Pickup File Size Or Format Error\n";}
   #The 8th line of each individuals genotype information should be the genotype at each UDG locus, matpi and mat non-pi. if it matches the individuals number there is a formating error  
  foreach $chr (0..($#chroms)) { #Write in genotype
  #loop through each chromosome and assign the genotype of the indiviual
   $startings[($indi*$pickup_lines)+($chr*2)+1] =~ /(\w+)/;
   $geno[$indi][$chr*2] = $1;
   $startings[($indi*$pickup_lines)+($chr*2)+2] =~ /(\w+)/;
   $geno[$indi][1+($chr*2)] = $1;
  }
  
  if ($geno[$indi][1] =~ /x/) {$sex[$indi] = 0;} #Decide sex; if the second copy chromosome (second homolog of the first chromosome is nothing but x's, individual is male.
  else {$sex[$indi] = 1;}
  
  $startings[(($indi+1)*$pickup_lines)-1] =~ /(\d.+\d)/;
  #snag the last line of genotype info for the individual
  ($u_gene1[$indi], $u_gene2[$indi], $mum_pi[$indi], $mum_nonpi[$indi]) = split "\t", $1;
  #split this line by tabbed values and assign to corresponding values to the udg 1 and 2 genotype, maternal pi and non-pi numbers.
 
 
 }
 push @result2, "Evolution Continued from $pickfile\n";
 #indicate to genotype stats array that this is a continued run


}

else {#Starting a new run
 foreach $indi (0..($popsize-1)) { #foreach individual
  $sex[$indi] = int(rand(2)); #choose a random number O or 1
  foreach $chr (0..($#chroms)) { #Set all genotype to uninfected; second sub is 0,1 for Sexual Chromosomes, 2,3; 4,5; 6,7; etc for Autosomes.
   $geno[$indi][$chr*2] = "0" x ($chroms[$chr]); #assign a value of zero to every site on each chromosome
   $geno[$indi][1+($chr*2)] = "0" x ($chroms[$chr]);
  }
  if ($sex[$indi] == 0) {$geno[$indi][1] = "x" x ($chroms[0])} #Males do not have the second sexual chromosome.
  $mum_pi[$indi] = 0; $mum_nonpi[$indi] = 0; #In the "-1" generation no TEs are present.
 
  #Set u-deciding genotype
  $random = rand(1); #choose a random number from 0to1
  foreach $allele (0..$#u_gene_list) { #for each allele
   $random = $random - $u_gene_freq[$allele]; #if the allele frequency is greater than then random number
   if ($random < 0)  {$u_gene1[$indi] = $u_gene_list[$allele]; last;} #this allele is UDG1
  }
  $random = rand(1);
  foreach $allele (0..$#u_gene_list) { #for each allele
   $random = $random - $u_gene_freq[$allele]; #if the allele frequency is greater than then random number
   if ($random < 0)  {$u_gene2[$indi] = $u_gene_list[$allele]; last;} #this allele is UDG1
  }
 

 }

 #Add the first TE copy here.
 if ((defined $ARGV[2]) && (defined $ARGV[3])) { #the chromosome and position are defined as part of the stdin.
  $chr_start = $ARGV[2]; $site_start = $ARGV[3];
 }
 else {
  $random = random_genome_site(0); # otherwise select a Random male site
  ($chr_start,$site_start) = split " ", $random;
  $chr_start = int($chr_start/2);
 }
 foreach $indi (0..((int($popsize/20))-1)) { #for the first 5% of the population
  substr $geno[$indi][2*$chr_start], $site_start, 1, "1"; #populate this site with an active TE.
 }
 push @result2, "Initial Site: $chr_start $site_start\n";
}
#die; #For testing the initialization;

#Evolving the population
foreach $gener (0..$gener_tot) {
 $sum_fit_male = 0; #starting values for male, fitness, female fitness, and X TEs are zero.
 $sum_fit_female = 0;
 $total_nte = 0;
 @telist = ();
 
 foreach $indi (0..($popsize-1)) { #Calculate fitness, fertility and transposition rate for each individual
  ($nte[$indi],$nte_inpi[$indi],$func_disr[$indi],$nte_act[$indi],$nte_inpiq[$indi]) = count_TE(@{$geno[$indi]}); #count_TE subroutine returnstotal #TEs, piTEs, multiplicative fitness costs of TEs, whether you have any active TEs remaining, & # of piqTEs 
  if ($nte[$indi] > 0) {@{$telist[$indi]} = @telist_temp;} #if there are TEs in the genome, extract their positions indicated by chromosome and site.
  else {@{$telist[$indi]} = ();} #otherwise the telist is empty
  
  $total_nte += $nte[$indi]; #add the TEs in the individual to the population total
  
  if ($model == 2 or $model == 3) {$rep[$indi] = 1;} #in models 2 and 3 there is no piRNA mediated repression
  elsif ($mum_pi[$indi] == 0) {$rep[$indi] = 1;} #if the mother had no insertions in piRNA clusters there is no piRNA mediated repression
  else {$rep[$indi] = exp(.8-($mum_pi[$indi] + $mum_nonpi[$indi]))}; #otherwise repression is given by this function
  if ($nte_act[$indi]) {$u[$indi] = get_u($nte[$indi],$rep[$indi],$u_gene1[$indi],$u_gene2[$indi],$nte_inpiq[$indi]);} #if there are transposase encoded TEs, calculate the transposition rate from # of transposable TEs, repression and UDG
  else {$u[$indi] = 0;} #If no active copy, u=0
  $fit[$indi] = fitness($nte[$indi],$rep[$indi],$func_disr[$indi],$sex[$indi],$nte_inpiq[$indi]); #fitness is determined based on ectopic recombination and TE insertions in functional sites 
  if ($nte_act[$indi]) {$fer[$indi] = 1-$rep[$indi]*fertility($nte[$indi]);} #fertility is determined based on active TEs and repression
  else {$fer[$indi] = 1;} #If no active copies, no hybrid dysgenesis
  if ($model == 1 or $model == 2 or $model == 5) {$fit_tot[$indi] = 1;} #in models 1, 2, and 5 there is no selection
  else {$fit_tot[$indi] = $fit[$indi]*$fer[$indi];} #in all other models, fitness differs between individuals according to the fitness effects above
  
  if ($sex[$indi] == 0) {$sum_fit_male += $fit_tot[$indi];} #add the individual's fitness to the fitness total for males, if male
  else {$sum_fit_female += $fit_tot[$indi];} #add the individual's fitness to the fitness total for female, if female
  
 }
 #End the program if all TE are lost
 if ($total_nte == 0) {print "All TE Copies lost at generation $gener\n"; $evo_gens = $gener; last;} #end program if not individuals remain with TEs
 if ($sum_fit_male == 0) {print "At generation $gener no fertile males are left\n"; $evo_gens = $gener; last;} #end program if no males remain
 if ($sum_fit_female == 0) {print "At generation $gener no fertile females are left\n"; $evo_gens = $gener; last;} #end program if not females remain
 
 #print "Generation $gener\n";
 
 pop_stats($gener); if (($gener > 0) && ($gener%20 == 0)) {geno_stats($gener);} #pop_stats calculates proportion of infected individuals, mean insertion, mean repression for every generation
 if ($snapshot[$gener] == 1) {all_geno($gener);} #saves all genotypes in the population if you have indicated to do so in your input file
 
 if ($gener == $gener_tot) {last;} #end simulation when you read the total input file.
 foreach $new (0..($popsize-1)) { #breeding
  $mating_m = choose_parent(0); #choose a male parents according to relative fitness
  $mating_f = choose_parent(1); $mum_pi_new[$new] = $nte_inpi[$mating_f]; $mum_nonpi_new[$new] = $nte[$mating_f] - $nte_inpi[$mating_f]; #choose a mating female according to relative fitness, save information about TE insertions inside and outside of piRNA clusters to calculate repression
  @sperm = make_sperm($mating_m); #make sperm based on make genotype
  @egg = make_egg($mating_f); #make egg based on female genotype
  foreach $chr (0..($#chroms)) { #create new genotype based on sperm and egg
   $geno_new[$new][$chr*2] = $egg[$chr];
   $geno_new[$new][1+($chr*2)] = $sperm[$chr]; #To make sure the missing X in male is always after the non-missing X
  }
  if ($geno_new[$new][1] =~ "x") {$sex_new[$new] = 0;} #determine sex
  else {$sex_new[$new] = 1;}
  
  $random = int(rand(2)); #Segragation and recombination of u-deciding gene alleles
  if ($random) {$u_gene1_new[$new] = $u_gene1[$mating_m];}
  else {$u_gene1_new[$new] = $u_gene2[$mating_m];}
  $random = int(rand(2));
  if ($random) {$u_gene2_new[$new] = $u_gene1[$mating_f];}
  else {$u_gene2_new[$new] = $u_gene2[$mating_f];}
 
 }

 #Writing the information for the new generation into main storage
 @geno = @geno_new; @geno_new = ();
 @sex = @sex_new; @sex_new = ();
 @mum_pi = @mum_pi_new; @mum_pi_new = ();
 @mum_nonpi = @mum_nonpi_new; @mum_nonpi_new = ();
 @u_gene1 = @u_gene1_new; @u_gene1_new = ();
 @u_gene2 = @u_gene2_new; @u_gene2_new = ();

}

#For testing the distribution of generations for total TE loss
#unless (defined $evo_gens) {$evo_gens = $gener_tot;}
#  open(FILEHANDLE, ">>Test_generation_count.txt");
#  print FILEHANDLE "$evo_gens\n";

#Write all genotypes into file if not lost
if (defined $evo_gens) {
 push @result, "Lost at $evo_gens\n";
}
else {all_geno($gener_tot);}
  
  open(FILEHANDLE, ">$outfile\_Population_stats.txt");
  print FILEHANDLE @result;
  open(FILEHANDLE, ">$outfile\_Genotype_stats.txt");
  print FILEHANDLE @result2;

  
  
  


sub inputing {
 #For initial testing purposes; remove when done
 #$popsize = 20;
 #@chroms = (20, 30);
 #$u_basic = 0.05;
 #$x = 0.0001;
 #$gener_tot = 1000;
 #$ect_basic = 0.0001;
 #$model = 0;
 #$pi[0] = ("1" x 5).("0" x 15);
 #$pi[1] = ("0" x 10).("1" x 5).("0" x 15);
 
 #$func[0] = ("0" x 10).("1" x 5).("0" x 5);
 #$func[1] = ("0" x 25).("1" x 5);
 
 #@u_gene_list = (0.5, 0.6, 0.8, 1, 1.2, 1.4, 1.5);
 #@u_gene_freq = (0.125, 0.125, 0.125, 0.25, 0.125, 0.125, 0.125);
 #foreach (0..10000) {$snapshot[$_] = 0;}
 #$pickup = 0;

 #Real part
 open (INPUTS, "$ARGV[0]") or die;
 my @infile = <INPUTS>;
 if (defined $ARGV[1]) {$outfile = $ARGV[1];}
 else {$outfile = "P_element";}
 #run the program with two arguments, first one is the exact file name, second one is the output file name
 $pickup = 0; #by default its a fresh run in less the input file indicates otherwise.
 my $input_switch = 0; #set stage in inputting process
 my $input_line;
 my $chr_sub = -1; #default value for chromosome number
 my $site_sub = -1; #default value for site
 my $udg_switch = 0; #default value for cycling through reading in the input file
 my $temp_sub;
 my $data_sub;
 my @b_sub; #initialize array.

 foreach $input_line (@infile) {
  chomp $input_line;
  if ($input_switch == 0) {#first step of inputting, identify popsize, max transposition rate, inactivation rate, #generations the simulation will run, ectopic recombination rate, maximum dysgenesis, what model you are using, what generations you'd like genotype snapshots for, and whether you are picking up from a previous simulation. 
   if ($input_line =~ /Popsize\s*=\s*([\d\.]+)/) {$popsize = $1;}
   elsif ($input_line =~ /Transposition_rate\s*=\s*([\d\.]+)/) {$u_basic = $1;}
   elsif ($input_line =~ /Inactivation_rate\s*=\s*([\d\.]+)/) {$x = $1;}
   elsif ($input_line =~ /Generation\s*=\s*([\d\.]+)/) {$gener_tot = $1; foreach (0..($gener_tot)) {$snapshot[$_] = 0;}} #make a snapshot array that is as long as the number of generations you will simulate
   elsif ($input_line =~ /Ectopic_rate\s*=\s*([\d\.]+)/) {$ect_basic = $1;}
   elsif ($input_line =~ /Max_dysgenesis\s*=\s*([\d\.]+)/) {$max_dys = $1;}
   elsif ($input_line =~ /Model\s*=\s*([\d\.]+)/) {$model = $1;}
   elsif ($input_line =~ /Snapshot\s*=\s*(\d.+\d)/) {@b_sub = split " ", $1; foreach (@b_sub) {$snapshot[$_] = 1;}}
   elsif ($input_line =~ /Pickup\s*=\s*(\w.+\w)/) {$pickup = 1; $pickfile = $1;}
   
   elsif ($input_line =~ /UDG/) {$input_switch = 1; next;} # when you reach the UDG line move to the next stage of inputting.
  }
  elsif ($input_switch == 1) {
   if ($udg_switch == 0) {
    @u_gene_list = split " ", $input_line;
    #make an array of the UDG gene values (i.e. multiplicative modifier of transposition)
	$udg_switch = 1; next;
   }
   if ($udg_switch == 1) {
    @u_gene_freq = split " ", $input_line;
     #make an array of the UDG gene frequencies
	$temp_sub = 0;
	foreach (0..($#u_gene_list-1)) {
	#for all except the last UDG gene sum the population frequencies of each gene
	 $temp_sub += $u_gene_freq[$_];
	}
	$u_gene_freq[$#u_gene_list] = 1-$temp_sub;
	#this ensures the frequency of UDG alleles sums to one, in case the user does not automatically specify it that way?
	$input_switch = 2; next;
   }
  }
  elsif ($input_switch == 2) {
  #this part deals with reading in the genome
   $input_line =~ s/\s//g;
   #remove the white space.
   if ($input_line =~ /^([a-q]+)\b/) {
   #for input lines that specify the site class of each position in a chromosome, collect the chromosome, ending with the word boundary
    $data_sub = $1; #assign chromosome to $data_sub
    $chr_sub += 1; #plus one to chromosome number
	$site_sub = -1; #starting site one chromosome is -1
	$pi[$chr_sub] = ""; #make a hash for pisites
	$func[$chr_sub] = ""; #make a hash for functional sites
	$chroms[$chr_sub] = length $data_sub; #assign the length of the chromosome to the $chroms hash, where the hashkey is the chromosome number
	foreach (0..($chroms[$chr_sub]-1)) {
	#foreach position on the chromosome
	 $temp_sub = substr $data_sub, $_, 1; #obtain the letter assigned to the position
	 if ($temp_sub eq "n") {$pi[$chr_sub] .= "0"; $func[$chr_sub] .= "0";} #if its a neutral site, the pisite scalar and functional site scalar for this chromosome are assigned a value of zero
	 elsif ($temp_sub eq "p") {$pi[$chr_sub] .= "1"; $func[$chr_sub] .= "0";} #if its a pi site, the pisite scalar is assigned a value of 1 at this position, and the functional site scalar is a assigned a value of zero
	 elsif ($temp_sub eq "q") {$pi[$chr_sub] .= "2"; $func[$chr_sub] .= "0";} #it is is a pseudo-pi, ths pisite scalar is assigned a value of 2 at this position, and the functional site scalar is assigned a value of zero.
	 else {$pi[$chr_sub] .= "0"; $func[$chr_sub] .= $temp_sub;} #if it is a functional site, the pisite scalar is assigned a value of zero and the functional site scalar is assigned the letter that occurs at this position
	}
   }
   elsif ($input_line =~ /^([\d\.+])/) {
   #this part deals with reading in the recombination rates
    $site_sub += 1; #add one to the site position, first value will be zero 
	$recom_rate[$chr_sub][$site_sub] = $input_line;
	#make an array of recombination rates assigned to the $recom_rate hash that is indexed by chromosome. the assigned recombination rate is 2X what is specific in the infile. why 2X?
   }
  }
 
 }
  
 #Add reading recombination rate from file
 $gsize_male = 0;
 $gsize_female = 0;
 foreach $chr_sub (0..$#chroms) {
 #specify the genome size differences beween males and females due to the absence of Y.
  $gsize_female += 2*($chroms[$chr_sub]);
  #female genome size is two times the length of each chromosome
  if ($chr_sub == 0) {$gsize_male += $chroms[$chr_sub];}
  else {$gsize_male += 2*($chroms[$chr_sub]);}
  #male genome size is 2*the length of each chromosome except for the X.
 
 }
 
 foreach $chr_sub (0..$#chroms) {
 #foreach chromosome
  $morgan[$chr_sub] = 0; #initiate a recombination rate for each chromosome
  foreach $site_sub (0..($chroms[$chr_sub]-2)) {
  #foreach space between sites in the chromosome
   $morgan[$chr_sub] += $recom_rate[$chr_sub][$site_sub];
   #add the calculated recombination rate to the total recombination rate for the chromosome
  }
 
 }
 
 
 $numchr = $#chroms+1; #number of chromosomes is equal to the number of elements in the chromosome array.

}
sub get_u {
#subroutime to calculate transposition rate per individual genome
#[0] nte [1] rep [2] u_gene1 [3] u_gene2 [4] nte_inpiq
 my $u_sub = $u_basic; #start with baseline recombination rate specified in the infile
 $u_sub = $u_sub*($_[2]+$_[3])/2; #multiplicatively modify UDG
 $u_sub = $u_sub*($_[0]-$_[4]); #multiply by number of non-piq elements
 $u_sub = $u_sub*(($_[1]+0.01)/1.01); #reduce based on the presence of piRNA silencing
 return $u_sub;

}



sub fitness {
#subroutine calculates individual fitness
#[0] nte [1] rep [2] func_disr [3] sex [4] nte_inpiq
 my $fit_sub = 1; #start with a fitness of one
 my $ectopic =0; #default ectopic recombination rate of zero
 $fit_sub = $fit_sub * $_[2]; #multiplicatively modify fitness based on insertions in functional sites
 if (($_[0]-$_[4] > 0) && ($_[3] == 1)) { #allow for ectopic recombination in females if there are any insertions outside of piRNA clusters
  $ectopic = $ect_basic * ($_[0]-$_[4]) * ($_[0]-$_[4]-1); #Ectopic recombination rate proportional to selection of two from non-pi TEs e*n*(n-1) ect basic is specific in the input file.
  unless ($model == 4 or $model == 5) {$ectopic = $ectopic*$_[1];} #pi-RNA silencing multiplicatively modifies ectopic recombination except for in models 4 and 5
 
 }
 if ($ectopic <= 1) {$fit_sub = $fit_sub * (1-$ectopic);} #fitness is reduced according to ectopic recombination rate
 else {$fit_sub = 0;} #fitness cost of ectopic recombination cannot be greater than 1.
 
 return $fit_sub;

}

sub fertility {
#calculates fertility fitness cost due to hybrid dysgenesis.
#[0] nte if nte_act > 0
my $fer_sub;
if ($_[0] < 4) {$fer_sub = 0;} #no fertility costs if fewer than 4 TEs
else {$fer_sub = $max_dys*($_[0]-3)/($_[0]-2);} #if 4 or more, fertility costs occurs according to max dysgenesis (D) and #of tes n; D*(n-3)/(n-2) 
return $fer_sub;
}



sub choose_parent {
#[0] sex of parent; 0 if male; 1 if female
#chose parents for each offspring based on their relative fitness?
 my $rand_sub; my $indi_sub; my $temp_sub; my $out;
 if ($_[0] == 0) {$rand_sub = rand(1) * $sum_fit_male;} #for males choose a random number between 0 and 1 and multiply by the summed fitness of all males in the population
 else {$rand_sub = rand(1) * $sum_fit_female;} #same as above for females
 $temp_sub = 0; #default value
 foreach $indi_sub (0..($popsize-1)) {
 #for every individual in the population
  if ($sex[$indi_sub] != $_[0]) {next;}
  #if their sex is not consistent with the sex of the parent you should be choosing, move on
  $temp_sub += $fit_tot[$indi_sub]; #otherwise, add the fitness of that individual to $temp_sub
  if ($temp_sub >= $rand_sub) {$out = $indi_sub; last;} #when temp_sub exceeds the value in $rand_sub, individual is selected, otherwise, the next individual is considered.
 }
 
 unless (defined $out) {die;} #end simulation if no parents left
 
 return $out; #return selected parent.
}
#Subroutine checked and finished

sub make_sperm { #No recombination in Drosophila males
 my @orig_geno = @{$geno[$_[0]]}; #extract the male genotype from the genotype array, one scalar per chromosome.
 my $rand_sub; my $rand2_sub; my $chr_sub; my $site_sub; my $state_sub; my @out_sub = ();
 my $num_transp; my $temp_sub;
 my $rand_site;
 my @list_sub = @{$telist[$_[0]]}; my @parent_te; my @new_te; #List of all non-pi TEs, used for selecting parent TE, format: "3 10" is chromosome 2 copy 1 site 10.
 my $status_sub;
 my $temp_nte = $nte[$_[0]]; #total # of TEs in the genotype

 #Excise and inactivate TEs and make new TE copies
 if ($u[$_[0]] > 0) { #if the mean transposition rate for this individual is greater than zero
  $num_transp = poisson ($u[$_[0]]); #determine the # of new insertions in the gamete, based on the transposition rate, assuming a poisson distribution.
  if ($num_transp > 0) {
   foreach (1..$num_transp) { #foreach new insertion
    $rand_sub = int(rand($#list_sub+1)); #choose a random integer that corresponds to an existing TE in the parent.
    @parent_te = split " ", $list_sub[$rand_sub]; #split parent TE into array [0] chromosome [1] site
	if ($temp_nte < $gsize_male) { #if there are sites available in the genome, select a random site.
	 do {$rand_site = random_genome_site(0); @new_te = split " ", $rand_site;} # select random sites in the male genome; split parent TE into array [0] chromosome [1] site; change to (1) in make egg
	 until ((substr $orig_geno[$new_te[0]], $new_te[1], 1)==0); #until you find an empty site
	 $status_sub = substr $orig_geno[$parent_te[0]], $parent_te[1], 1; #Read parent TE status, 1 is active, 2 is inactive
	 substr $orig_geno[$new_te[0]], $new_te[1], 1, $status_sub; #new TE status is determined by parent TE status
	 if ((substr $pi[int($new_te[0]/2)], $new_te[1], 1) == 0) {push @list_sub, $rand_site;} #Write into new site if not a pi-TE #Write new insertion into active TE list if not piRNA cluster
	 $temp_nte += 1; #add the TE to the total count.
	}
	$rand2_sub = rand(); #select a random number from zero to 1.
	if ($rand2_sub < 0.1) { #Excision in 10% of the time
	 substr $orig_geno[$parent_te[0]], $parent_te[1], 1, "0"; #if it is less than 0.1, remove the parent element.
	 splice @list_sub, $rand_sub, 1; #also remove the parent TE from the list of transposable insertions.
	 unless (defined $list_sub[0]) {last;} #exit if this removes all mobile TEs from the genome? 
	 $temp_nte -= 1; #subtract from the total count.
	}
   }  
  }
 }
 if (defined $list_sub[0]) { #Inactivation process
  foreach $temp_sub (0..$#list_sub) {#for every active TE
   @parent_te = split " ", $list_sub[$temp_sub];
   $status_sub = substr $orig_geno[$parent_te[0]], $parent_te[1], 1;
   if ($status_sub == 1) {#if TE is transposase encoding
    $rand2_sub = rand(); #select a random number zero to 1
	if ($rand2_sub < $x) { #if its less than the specified inactivation rate
	 substr $orig_geno[$parent_te[0]], $parent_te[1], 1, "2"; #Inactivate the TE at the rate specified in the input. #It would be appropriate to have this tied to the transposing element, but this is probably fine.
	}
   }
  }
 }
 
 #Meiosis
 foreach $chr_sub (0..$#chroms) {
  $rand_sub = int(rand(2));#choose a random integer, zero or 1.
  if ($rand_sub) {push @out_sub, $orig_geno[$chr_sub*2];} #if integer is 1 choose second homolog add to gamete
  else {push @out_sub, $orig_geno[($chr_sub*2)+1];} #if integer is 1 choose first homolog, add to gamete
 }
 return @out_sub;
 
}
#Subroutine checked and finished

sub make_egg { #same as make sperm EXCEPT uses a female genome and allows for recombination
 my @orig_geno = @{$geno[$_[0]]};
 my $rand_sub; my $rand2_sub; my $chr_sub; my $site_sub; my $state_sub; my @out_sub = ();
 my $num_transp; my $temp_sub;
 my $rand_site;
 my @list_sub = @{$telist[$_[0]]}; my @parent_te; my @new_te;
 my $status_sub;
 my $recom_num;
 my $temp_nte = $nte[$_[0]];
 
 #Excise and inactivate TEs and make new TE copies
 if ($u[$_[0]] > 0) {
  $num_transp = poisson ($u[$_[0]]);
  if ($num_transp > 0) {
   foreach (1..$num_transp) {
    $rand_sub = int(rand($#list_sub+1));
    @parent_te = split " ", $list_sub[$rand_sub]; #Choose parent TE
	if ($temp_nte < $gsize_female) {
	 do {$rand_site = random_genome_site(1); @new_te = split " ", $rand_site;}
     until ((substr $orig_geno[$new_te[0]], $new_te[1], 1)==0); #Choose new TE site
     $status_sub = substr $orig_geno[$parent_te[0]], $parent_te[1], 1; #Read parent TE status 
     substr $orig_geno[$new_te[0]], $new_te[1], 1, $status_sub;
     if ((substr $pi[int($new_te[0]/2)], $new_te[1], 1) == 0) {push @list_sub, $rand_site;} #Write into new site if not a pi-TE
	 $temp_nte+= 1;
    }
	$rand2_sub = rand();
	if ($rand2_sub < 0.1) { #Excision in 10% of the time
	 substr $orig_geno[$parent_te[0]], $parent_te[1], 1, "0";
	 splice @list_sub, $rand_sub, 1;
	 unless (defined $list_sub[0]) {last;}
	 $temp_nte -= 1;
	}
   }  
  }
 }
 if (defined $list_sub[0]) { #Inactivation process
  foreach $temp_sub (0..$#list_sub) {
   @parent_te = split " ", $list_sub[$temp_sub];
   $status_sub = substr $orig_geno[$parent_te[0]], $parent_te[1], 1;
   if ($status_sub == 1) {
    $rand2_sub = rand();
	if ($rand2_sub < $x) {
	 substr $orig_geno[$parent_te[0]], $parent_te[1], 1, "2"; #Inactivate the TE
	}
   }
  }
 }
 
 #Recombination - unique to females
 my $recom_temp; my $recom_temp2;
 foreach $chr_sub (0..$#chroms) { #for each chromosome
  $recom_num = poisson ($morgan[$chr_sub]); #the number of recombination events is drawn from a poisson distribution with mean equal to total recombination rate for the chromosome.
  if ($recom_num == 0) {next;} #if no crossovers, move on to next chromosome
  foreach (1..$recom_num) { #otherwise for each recombinant
   $rand_sub = rand($morgan[$chr_sub]); #choose a random fraction from zero to the total recombination rate for the chromosome
   $site_sub = -1; #set the default starting site. 
   do {
    $site_sub += 1;
	$rand_sub -= $recom_rate[$chr_sub][$site_sub]; #subtract the recombination rate at each site from the randomly selected value
	if ($site_sub > $chroms[$chr_sub]) {die;} #error if sites exceed chromosome length.
   }
   until ($rand_sub < 0); #until it is less than zero
   $recom_temp = substr $orig_geno[$chr_sub*2], $site_sub; #collect the sites after the recombination site on homolog 1
   $recom_temp2 = substr $orig_geno[($chr_sub*2)+1], $site_sub; #same for homolog 2
   $orig_geno[$chr_sub*2] = (substr $orig_geno[$chr_sub*2], 0, $site_sub).$recom_temp2; #concatenate beginning of homolog 1 with homolog 2
   $orig_geno[($chr_sub*2)+1] = (substr $orig_geno[($chr_sub*2)+1], 0, $site_sub).$recom_temp; #reciprocal
  }
 }
 
 #Meiosis
 foreach $chr_sub (0..$#chroms) {#same as male meiosis
  $rand_sub = int(rand(2));
  if ($rand_sub) {push @out_sub, $orig_geno[$chr_sub*2];}
  else {push @out_sub, $orig_geno[($chr_sub*2)+1];}
 }
 return @out_sub;
 
}
#Subroutine checked and finished

sub count_TE {
#[0] an individuals genotype array, where each scalar indicates the genotype of a particular chromosome
 my $chr_sub; my $site_sub; my $total_te = 0; my $pi_te = 0; my $piq_te = 0; my $live_te = 0; my $fit_te = 1; #default values
 my $zyg_sub; my $func_sub;
 @telist_temp = (); #clear out this array each time
 my $temp_sub;
 
 
 foreach $chr_sub (0..$#chroms) { #foreach chromosome
  foreach $site_sub (0..($chroms[$chr_sub]-1)) { #for each site in the chromosome
   $zyg_sub = 0;
   if (substr $_[$chr_sub*2], $site_sub, 1) { #if site is occupied, as indicated by a value other than zero
    $total_te += 1; #count this TE towards the total # of TEs
	if ((substr $pi[$chr_sub], $site_sub, 1) == 1) {$pi_te += 1; $piq_te += 1;} #if its a pi site it also counts towards piTES
	elsif ((substr $pi[$chr_sub], $site_sub, 1) == 2) {$piq_te += 1;} #if its a pseudo-pi TE, it counts towards pseudo-piTEs
	else {
	 $temp_sub = ($chr_sub*2)." $site_sub";
	 push @telist_temp, $temp_sub; #otherwise it is added to the list of potential parent TEs for transposition
	}
	$zyg_sub += 1; #add to counter indicate zygosity
   }
   if ($chr_sub == 0 && $_[1] =~ /x/) {$zyg_sub *= 2;} #the x-chromosome is hemizygous in males, so fitness effects are revealed.
   elsif (substr $_[($chr_sub*2)+1], $site_sub, 1) { #if site on other homolog is occupied, as indicated by a value other than zero
    $total_te += 1; #count this TE towards the total # of TEs
	if ((substr $pi[$chr_sub], $site_sub, 1) == 1) {$pi_te += 1; $piq_te += 1;} #if its a pi site it also counts towards piTES
	elsif ((substr $pi[$chr_sub], $site_sub, 1) == 2) {$piq_te += 1;} #if its a pseudo-pi TE, it counts towards pseudo-piTEs
	else {
	 $temp_sub = (($chr_sub*2)+1)." $site_sub";
	 push @telist_temp, $temp_sub; #otherwise it is added to the list of potential parent TEs for transposition
	}
	$zyg_sub += 1; #count towards zygosity
   }
   
   if (((substr $func[$chr_sub], $site_sub, 1) ne "0") and ($zyg_sub > 0)) { #if this site is function and it is occupied by a TE, note that dominant sites aren't present in our genome.
    $func_sub = (substr $func[$chr_sub], $site_sub, 1);
	if ($func_sub eq "a") {$fit_te = 0;} #fitness effects for class a sites (dominant)
	if (($func_sub eq "b") and ($zyg_sub == 2)) {$fit_te = 0;} #fitness effects class b sites (recessive)
	if ($func_sub eq "c") {$fit_te *= 0.5;} #fitness effects for class c sites (dominant)
	if (($func_sub eq "d") and ($zyg_sub == 2)) {$fit_te *= 0.75;} #fitness effects class d sites (recessive)
	if ($func_sub eq "e") {$fit_te *= 0.9;} #fitness effects for class e sites (dominant)
	if (($func_sub eq "f") and ($zyg_sub == 2)) {$fit_te *= 0.95;} #fitness effects class f sites (recessive)
	if ($func_sub eq "g") {$fit_te *= 0.95;} #fitness effects for class g sites (dominant)
	if (($func_sub eq "h") and ($zyg_sub == 2)) {$fit_te *= 0.99;} #fitness effects class h sites (recessive)
	if ($func_sub eq "i") {$fit_te *= 1.05;} #fitness effects for class i sites (dominant)
	if (($func_sub eq "j") and ($zyg_sub == 2)) {$fit_te *= 1.01;} #fitness effects class j sites (recessive)
   
   }
   
  }
 }
 
 if ($total_te - $piq_te > 0) { #if there are TEs that aren't in piq sites (i.e TEs that can move)
  foreach $temp_sub (@telist_temp) { #for each insertion that can move
   ($chr_sub, $site_sub) = split " ", $temp_sub; #extract the chromosome position
   if (((substr $pi[int($chr_sub/2)], $site_sub, 1) == 0)&&((substr $_[$chr_sub], $site_sub, 1) == 1)) {$live_te = 1; last;} #if its not in a pi site, you have at least one transposase encoding TE in the genome
  }
 }
 
 
 return ($total_te, $pi_te, $fit_te, $live_te, $piq_te); #return the total #TEs, piTEs, multiplicative fitness costs of TEs, whether you have any active TEs remaining, &E of piqTEs
}
#Subroutine checked and finished

sub random_genome_site {
#[0] sex of parent; 0 if male; 1 if female
 my $rand_sub2; my $chr_sub2;
 if ($_[0] == 0) {$rand_sub2 = int(rand($gsize_male));} #select a random integer indicating a position in the male genome
 else {$rand_sub2 = int(rand($gsize_female));} #select a random integer indicating a position in the female genome
 foreach $chr_sub2 (0..$#chroms) { #foreach chromosome
  if ($rand_sub2 < $chroms[$chr_sub2]) { #if the site is on the current chromosome
   return ($chr_sub2*2)." $rand_sub2"; #return chromosome and position of randomly selected new insertion & exit subroutine
  }
  else {$rand_sub2 -= $chroms[$chr_sub2];} #if the site wasn't on the current chromosome, subtract the number of genomic positions on that chromosome
  
  if (($_[0] == 0)&&($chr_sub2 == 0)) {next;} #for males, skip the second homolog of the X-chromosome 
  if ($rand_sub2 < $chroms[$chr_sub2]) { #if the site is on the current chromosome
   return (($chr_sub2*2)+1)." $rand_sub2"; #return chromosome and position of randomly selected new insertion & exit subroutine
  }
  else {$rand_sub2 -= $chroms[$chr_sub2];} #if the site wasn't on the current chromosome, subtract the number of genomic positions on that chromosome
 }
 print "Random genome site error\n"; die; #give an error if the site specified isn't in the genome.
}
#Subroutine checked and finished

sub pop_stats {
# [0] generation
 push @result, "Generation $_[0]\n";
 my @stats_sub; my $perc_sub; my $temp_sub;
 @stats_sub = n_mean_sd(@nte); #@stats = [0] n [1] mean # TEs [2] sd # of TEs
 push @result, "Number of TE\t$stats_sub[1]\t$stats_sub[2]\n"; #add to results array
 $perc_sub = 0;
 foreach $temp_sub (@nte) {if ($temp_sub >= 1) {$perc_sub += 1;}} #count individuals with one or more TEs
 $perc_sub = $perc_sub/$popsize; #determine what fraction of the population they are
 push @result, "Proportion infected\t$perc_sub\n"; #add proportion infected to result array
 $perc_sub = 0; 
 foreach $temp_sub (@nte_act) {if ($temp_sub >= 1) {$perc_sub += 1;}} #count individuals with one or more active TEs
 $perc_sub = $perc_sub/$popsize; #determine what fraction of the population they are
 push @result, "Proportion able to transpose\t$perc_sub\n"; #add proportion infected to result array
 
 
 @stats_sub = n_mean_sd(@nte_inpi); #@stats = [0] n [1] mean # TEs in pi sites [2] sd # of TEs pisites
 push @result, "Number of piTE\t$stats_sub[1]\t$stats_sub[2]\n";
 @stats_sub = n_mean_sd(@nte_inpiq); #add to results array
 push @result, "Number of piTE and pseudo-piTE\t$stats_sub[1]\t$stats_sub[2]\n"; #@stats = [0] n [1] mean # TEs in piq sites [2] sd # of TEs piq sites
 
 @stats_sub = n_mean_sd(@fit_tot); #@stats = [0] n [1] mean fitness [2] sd fitness
 push @result, "Total fitness\t$stats_sub[1]\t$stats_sub[2]\n"; #add to results
 @stats_sub = n_mean_sd(@u); #@stats = [0] n [1] mean transposition rate (per genome) [2] sd transposition rate (per genome)
 push @result, "Transposition rate\t$stats_sub[1]\t$stats_sub[2]\n";
@stats_sub = n_mean_sd(@rep); #@stats = [0] n [1] mean repression [2] sd repression
 push @result, "Maternal Repression\t$stats_sub[1]\t$stats_sub[2]\n";
 my %udgfreq_sub;
 foreach $temp_sub (@u_gene1) {$udgfreq_sub{$temp_sub} += 1;} 
 foreach $temp_sub (@u_gene2) {$udgfreq_sub{$temp_sub} += 1;} #calculate the total # of each UDG allele
 
 foreach $temp_sub (@u_gene_list) {
  unless (defined ($udgfreq_sub{$temp_sub})) {$udgfreq_sub{$temp_sub} = 0;}
  push @result, "$temp_sub\t"; #add allele to result array
  push @result, ($udgfreq_sub{$temp_sub})/(2*$popsize); #add frequency to result array
  push @result, "\n";
 }
 

}
#Subroutine checked and finished

sub geno_stats { #every 20 generations, [0] generation number
 my @geno_sub;
 my @geno2_sub;
 my $chr_sub;
 my $site_sub;
 my $indi_sub;
 my @gtype_sub;
 my $fems_sub = 0;
 foreach $indi_sub (0..($popsize-1)) {#forevery individual
  if ($sex[$indi_sub] == 1) {$fems_sub += 1;} #count them if they are female
 }
 push @result2, "Generation $_[0]\nX Chromosome Total\t";
 push @result2, $popsize+$fems_sub; #X chromosomes = popsize + # of females
 push @result2, "\nAutosome Total\t";
 push @result2, $popsize*2; #autosomes is twice population size
 push @result2, "\n";

 
 foreach $chr_sub (0..$#chroms) { #foreach chromosome
  foreach $site_sub (0..($chroms[$chr_sub]-1)) { #for each site
   $geno_sub[$chr_sub][$site_sub] = 0; 
   $geno2_sub[$chr_sub][$site_sub] = 0; #make arrays with starting values of zero
  }
 }
 
 foreach $indi_sub (0..($popsize-1)) { #foreach individual
  @gtype_sub = @{$geno[$indi_sub]}; #extract their full genotype
  foreach $chr_sub (0..$#chroms) { #foreach chromosome
   foreach $site_sub (0..($chroms[$chr_sub]-1)) { #foreach site
    if ((substr $gtype_sub[$chr_sub*2], $site_sub, 1) == 1) {
	 $geno_sub[$chr_sub][$site_sub] += 1; #add one to the genotype array if site is occupied by active TE
    }
    elsif ((substr $gtype_sub[$chr_sub*2], $site_sub, 1) == 2) {
	 $geno2_sub[$chr_sub][$site_sub] += 1; #add one to the genotype2 array if site is occupied by inactive TE
    }
    if ($chr_sub == 0 && $sex[$indi_sub] == 0) {next;} #move to next chromosome if this is the second copy of the X-chromosome in males.
    if ((substr $gtype_sub[($chr_sub*2)+1], $site_sub, 1) == 1) {
	 $geno_sub[$chr_sub][$site_sub] += 1; #add one to the genotype array if site is occupied by active TE
    }
    elsif ((substr $gtype_sub[($chr_sub*2)+1], $site_sub, 1) == 2) {
	 $geno2_sub[$chr_sub][$site_sub] += 1; #add one to the genotype2 array if site is occupied by inactive TE
    }
   }
  }
 }
 
 foreach $chr_sub (0..$#chroms) { #foreach chromosome
  foreach $site_sub (0..($chroms[$chr_sub]-1)) { #for each site
   if ($geno_sub[$chr_sub][$site_sub] + $geno2_sub[$chr_sub][$site_sub] == 0) {next;} #if there are no insertions, move to the next site.
   push @result2, "$chr_sub\t$site_sub\t$geno_sub[$chr_sub][$site_sub]\t$geno2_sub[$chr_sub][$site_sub]\n"; #otherwise add the chromosome, site, #of active insertions and # of inactive insertions to the result array
  }
 }
 

}
#Subroutine checked and finished

sub all_geno {
 @result3 = ();
 foreach $indi (0..($popsize-1)) {
 #cycle through every individual in the population
  push @result3, "$indi\n";
  #add their numer to the result vector
  foreach $chr (0..($#chroms)) {# second sub is 0,1 for Sexual Chromosomes, 2,3; 4,5; 6,7; etc for Autosomes.
  #cycle through the each chromosome pair
   push @result3, $geno[$indi][$chr*2]; push @result3, "\n";
   #for first homolog, push the full chromosomal genotype for first homolog onto $result3 vector. Every site is a value of 0,1 or 2 indicating the number of TEs occupying that site.
   push @result3, $geno[$indi][1+($chr*2)]; push @result3, "\n";
   #same for second homolog
 }
  push @result3, "$u_gene1[$indi]\t$u_gene2[$indi]\t$mum_pi[$indi]\t$mum_nonpi[$indi]\n";
  #push onto results vector the genotypes for the multiplicative modifier alleles (always 1 in this study), the number of insertions in piRNA clusters in the maternal genotype, and the number of insertions outside of piRNA clusters in the maternal genotype, which determines repression in offspring.
 }
 #added by erin to ensure 99% infection at generation 500
 @props = grep (/Proportion infected/, @result);
 $props[2000] =~ m/\t(.+?)\n$/g;
 if ($1 >= 0.99)
 	{
 	open(FILEHANDLE, ">$outfile\_Snapshot_Generation_$_[0].txt");
	print FILEHANDLE @result3;
 	print "success";
	}
else
	{
	print "fail";
	}
}	
#Subroutine checked and finished

sub n_mean_sd {
#vector of values from individuals (nte, pite, piqte etc)
 my $n_st = 0; my $mean_st = 0; my $var_st = 0; my $sd_st = 0; my $temp_st;
 my @data_st = @_;
 
 foreach $temp_st (@data_st) {#for reach value
  $n_st += 1; #count of individuals
  $mean_st += $temp_st; #add total sum
 }
 $mean_st = $mean_st/$n_st; #calculate mean

 foreach $temp_st (@data_st) {#for each value
  $var_st += ($temp_st-$mean_st)*($temp_st-$mean_st); #calculate some of squares
 }
 $var_st = $var_st/($n_st-1); #divide by n-1 to calculate variance
 $sd_st = sqrt ($var_st); #calculate sd
 
 return ($n_st, $mean_st, $sd_st); #return popsize, mean and sd.
 
 

}
#Subroutine checked and finished

sub poisson { #determines the number of new insertions in each gamete from a poisson probability distribution based on the transposition rate.
#[0] mean transposition rate; 
 my $rand_pois;
 my $excl_pois = 1; #set starting value for factorials of # of successes
 my $out_pois = -1; #set starting value for # of successes
 my $prob_pois = 0; #set starting value for cumulative poisson probability.
 
 $rand_pois = rand(1); #select a random number from zero to 1
 
 do {
  $out_pois += 1; #add 1
  if ($out_pois > 0) {$excl_pois *= $out_pois;} #calculate factorial for number of successes
  $prob_pois += ($_[0]**$out_pois)*(exp(0-$_[0]))/$excl_pois; #cumulative poisson probability function 
 }
 until ($prob_pois > $rand_pois); #continue to calculate the cumulative poisson probability until it exceeds the random number that you selected.
 
 return $out_pois; #return the #of successes required to exceed the random number. Because the probability "window" of each number of successes depends on the transposition rate, the wider the window, the greater the probability that x falls in that window. 
}
#Subroutine checked and finished

sub debug {
   open(FILEHANDLE, ">DEBUG.txt");
  print FILEHANDLE @_;

}