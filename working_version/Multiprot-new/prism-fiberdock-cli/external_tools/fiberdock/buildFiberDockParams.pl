#!/usr/bin/perl -w

use strict;
use FindBin;

if ( $#ARGV < 12 ) {
  print "buildFiberDockParams.pl <receptorPDB> <ligandPDB> <(Un)bound> <(Un)bound> <complexType> <trans> <outFile> <RISCO/FISCO> <numRBOCycles> <radiiScale> <output refined PDBs> <lpMethod> <paramsFile> [optional:  <lambda> <recAnalyzedFile> <recNMs>] [optional: <ligAnalyzedFile> <ligNMs>]\n";
  print "\n";
  print "Help: \n";
  print "\n";
  print "(Un)bound        {B,U}          - B if the receptor is bound, U if the receptor is unbound.\n";
  print "complexType      {EI,AA,PR,     - EI for enzyme-inhibitor comlexes, AA for antibody-antigen, \n";
  print "                  Default}        PR for protein-RNA and Default for others.\n";
  print "trans            fileName       - Transformations of the input ligand.\n";
  print "outFile          fileName       - Output file for scoring.\n";
  print "RISCO/FISCO      {0,1}          - 0 - Restricted Side-Chain Optimization,\n";
  print "                                  1 - Full Side-Chain Optimization.\n";
  print "numRBOCycles     {unsigned int} - Number of Monte Carlo cycles in rigid-body minimization.\n";
  print "radiiScale       {0,1}          - Radii Scaling for soft rigid-body optimization and scoring.\n";
  print "output complexes {0,1}          - Output refined PDBs.\n";
  print "lpMethod         {cplex,glpk}   - choose lp-solver. \n";
  print "                                  You can use cplex ONLY if it is installed on your computer.\n";
  print "paramsFile    	  fileName       - Name of the parameters file.\n";
  print "lambda           [0-1]          - Scaling factor for the backbone refinement score.\n";
  print "                                  The lower it is the more backbone movement is allowed.\n";
  print "                                  Recomanded value: 0.05.\n";
  print "recAnalyzedFile  fileName       - The name of the receptor's PDB file on which the NMA was performed,\n";
  print "                                  For example - if the NMA was performed only on CA atoms\n";
  print "                                  then the file should contain only CA atoms.\n";
  print "recNMs           fileName       - The file with the calculated normal modes.\n";
  print "ligAnalyzedFile  fileName       - Same as recAnalyzedFile, but for the ligand\n";
  print "ligNMs           fileName       - Same as recNMs, but for the ligand\n";
  print "\n";
  print "Examples: \n";
  print "\n";
  print "1. Example for unbound enzyme-inhibitor. Here we run the coarse refinement (RISCO,radiiScale=0.8)\n";
  print "   buildFiberDockParams.pl 1ACB_r.pdb.CHB.pdb 1ACB_l.pdb.CHB.pdb U U EI 1ACB.zdock.trans 1ACB.firedock.out 0 50 0.8 0 cplex FiberDock_params.txt 0.05 1ACB_r.ca.pdb 1ACB_r.pdb 1ACB_l.ca.pdb 1ACB_l.pdb\n";
  print "\n";  
  print "2. Example for semi-unbound antibody-antigen. PatchDock swaps between the antibody and antigen. \n";
  print "   So, we should refer to them as PatchDock did. One of them is given bound. Here we run the full refinement on 25 top candidates\n"; 
  print "   (FISCO,radiiScale=0.85). In addition, we print the refined complexes.\n";
  print "   buildFiberDockParams.pl 2VIR_l_u.pdb.CHB.pdb 2VIR_r_u.pdb.CHB.pdb U B AA trans_top25.txt out 1 50 0.85 1 cplex FiberDock_params.txt 0.05 1ACB_r.ca.pdb 1ACB_r.pdb 1ACB_l.ca.pdb 1ACB_l.pdb\n";
  print "\n";
  exit;
}

my $home = "$FindBin::Bin";

my $recPdb = $ARGV[0];
my $ligPdb = $ARGV[1];
my $boundrec = $ARGV[2];
my $boundlig = $ARGV[3];
my $complexType = $ARGV[4];
my $trans = $ARGV[5];
my $energiesOut = $ARGV[6];
my $fullSCP = $ARGV[7];
my $rbmCycles = $ARGV[8];
my $radiiScale = $ARGV[9];
my $print = $ARGV[10];
my $RecFixedFile = "none";
my $LigFixedFile = "none";
my $RecFlexFile = "none";
my $LigFlexFile = "none";
my $lpMethod = $ARGV[11];
my $paramsFile = $ARGV[12];
my $RecAnalyzedFile;
my $ReceptorNMs;
my $LigAnalyzedFile;
my $LigandNMs;
my $lambda;

if ($#ARGV >12){
  $lambda = $ARGV[13];
  $RecAnalyzedFile = $ARGV[14];
  $ReceptorNMs = $ARGV[15];
}

if ($#ARGV >15){
  $LigAnalyzedFile = $ARGV[16];
  $LigandNMs = $ARGV[17];
}
open OUT, ">$paramsFile";

print OUT "### I/O\n";
print OUT "receptorPDBFileName $recPdb\n";
print OUT "ligandPDBFileName $ligPdb\n";
print OUT "## reference for rmsd calculations \n";
print OUT "# templateLigandPDBFileName $ligPdb\n";
print OUT "# transformations for refinement\n";
print OUT "transFileName $trans\n";
print OUT "## libraries files\n";
print OUT "rotamerLibFile  $home/lib/bbdep02.May.sortlib\n";
print OUT "protLib $home/lib/chem.lib\n";
print OUT "# pdbConventionFile $home/lib/Names.CHARMM.db\n";
print OUT "## output file\n";
print OUT "energiesOutFileName $energiesOut\n";
print OUT "\n";
print OUT "### Input Options\n";
print OUT "## these residues will be fixed unless appear also in receptorFlexibleResiduesFile\n";

if ($RecFixedFile eq "none") {
  print OUT "# receptorFixedResiduesFile rec_fixed_residues.txt\n";
} else {
  print OUT "  receptorFixedResiduesFile $RecFixedFile\n";
}
if ($LigFixedFile eq "none") {
  print OUT "# ligandFixedResiduesFile lig_fixed_residues.txt\n";
} else {
  print OUT "  ligandFixedResiduesFile $LigFixedFile\n";
}

print OUT "## defines for which chains to build surface residues rotamers. If empty than will build for all chains\n";
print OUT "## For bound receptor, you should uncomment it\n";
if ($boundrec eq "B") {
  print OUT "flexibleReceptorChains @\n";
} else {
  print OUT "# flexibleReceptorChains @\n";
}
print OUT "## For bound ligand, you should uncomment it \n";
if ($boundlig eq "B") {
  print OUT "flexibleLigandChains @\n";
} else {
  print OUT "# flexibleLigandChains @\n";
}
print OUT "## Residues specified here will be flexible\n";
print OUT "## The residue can be flexible if its chain was specified in flexibleReceptorChains or if flexibleReceptorChains is empty\n";

if ($RecFlexFile eq "none") {
  print OUT "# receptorFlexibleResiduesFile rec_flex_residues.txt\n";
} else {
  print OUT "  receptorFlexibleResiduesFile $RecFlexFile\n";
}
if ($LigFlexFile eq "none") {
  print OUT "# ligandFlexibleResiduesFile lig_flex_residues.txt\n";
} else {
  print OUT "  ligandFlexibleResiduesFile $LigFlexFile\n";
}

print OUT "\n";
print OUT "### Output Options\n";
print OUT "# to output refined complexes\n";
print OUT "printRefinedComplexes $print\n";
print OUT "## 1 - only energy caclulaltion is performed without refinement (works only for FiberDock pre)\n";
print OUT "onlyEnergyCalculation 0\n";

print OUT "\n";
if ($fullSCP eq "1") {
  print OUT "### side-chain optimization\n";
  print OUT "# 1 - only clashing residues are flexible, 0 - all residues are flexible\n";
  print OUT "receptorOnlyClashesMovable 0\n";
  print OUT "ligandOnlyClashesMovable 0\n";
  print OUT "# 0 - small rotamer set, 1 - extended rotamer set\n";
  print OUT "extraRotamers 1\n";
} else {
  print OUT "### side-chain optimization\n";
  print OUT "# 1 - only clashing residues are flexible, 0 - all residues are flexible\n";
  print OUT "receptorOnlyClashesMovable 1\n";
  print OUT "ligandOnlyClashesMovable 1\n";
  print OUT "# 0 - small rotamer set, 1 - extended rotamer set\n";
  print OUT "extraRotamers 0\n";
}
if ($lpMethod eq "glpk"){
  print OUT "lpMethod 2\n";
}
elsif ($lpMethod eq "cplex"){
  print OUT "lpMethod 1\n";
}

print OUT "exclusionsFile $home/lib/exclusions\n";

print OUT "\n";
print OUT "### rigid-body optimization (RBO)\n";
print OUT "# num of MC cycles for rigid-body minimization (if 0 - no RBO)\n";
print OUT "rigidBodyMinimizationCycles $rbmCycles\n";

print OUT "charmmParamFile $home/lib/par_all27_prot_na.prm\n";
print OUT "definitionsFile $home/lib/top_all27_prot_na.rtf\n";

print OUT "\n";

if ($complexType eq "EI") {
  print OUT "### weights for energy score (Enzyme-Inhibitor)\n";
  print OUT "attrVdWWeight   1.0\n";
  print OUT "repVdWWeight    0.95\n";
  print OUT "ACEWeight       1.6\n";
  print OUT "attrElWeight    0.07\n";
  print OUT "repElWeight     0.12\n";
  print OUT "l_attrElWeight  0.0\n";
  print OUT "l_repElWeight   0.3\n";
  print OUT "HBWeight        1.32\n";
  print OUT "pipiWeight      1.0\n";
  print OUT "catpiWeight     0.8\n";
  print OUT "aliphWeight     0.5\n";
  print OUT "insidenessWeight        1.55\n";
  print OUT "confProbWeight  0.0\n";

} elsif ($complexType eq "AA") {
  print OUT "### weights for energy score (Antibidy-Antigen)\n";
  print OUT "attrVdWWeight   1.5\n";
  print OUT "repVdWWeight    0.6\n";
  print OUT "ACEWeight       1.6\n";
  print OUT "attrElWeight    0.21\n";
  print OUT "repElWeight     0.21\n";
  print OUT "l_attrElWeight  0.46\n";
  print OUT "l_repElWeight   0.69\n";
  print OUT "HBWeight        1.2\n";
  print OUT "pipiWeight      1.0\n";
  print OUT "catpiWeight     0.7\n";
  print OUT "aliphWeight     2.5\n";
  print OUT "insidenessWeight        0.0\n";
  print OUT "confProbWeight  0.0\n";

} elsif ($complexType eq "PR") {
  print OUT "### weights for energy score (Protein-RNA)\n";
  print OUT "attrVdWWeight   1.5\n";
  print OUT "repVdWWeight    0.8\n";
  print OUT "ACEWeight       0.0\n";
  print OUT "attrElWeight    0.021\n";
  print OUT "repElWeight     0.015\n";
  print OUT "l_attrElWeight  0.069\n";
  print OUT "l_repElWeight   0.04\n";
  print OUT "HBWeight        3.0\n";
  print OUT "pipiWeight      1.0\n";
  print OUT "catpiWeight     0.0\n";
  print OUT "aliphWeight     0.0\n";
  print OUT "insidenessWeight        0.7\n";
  print OUT "confProbWeight  0.0\n";

} else {
  print OUT "### weights for energy score (Default)\n";
  print OUT "attrVdWWeight   1.5\n";
  print OUT "repVdWWeight    0.8\n";
  print OUT "ACEWeight       1.6\n";
  print OUT "attrElWeight    0.21\n";
  print OUT "repElWeight     0.21\n";
  print OUT "l_attrElWeight  0.0\n";
  print OUT "l_repElWeight   0.69\n";
  print OUT "HBWeight        1.2\n";
  print OUT "pipiWeight      1.0\n";
  print OUT "catpiWeight     0.7\n";
  print OUT "aliphWeight     2.5\n";
  print OUT "insidenessWeight        0.7\n";
  print OUT "confProbWeight  0.0\n";
}
print OUT "radiiScaling $radiiScale\n";

print OUT "\n\n### Backbone refinement parameters\n";
if ($#ARGV >12){
  print OUT "performBackboneRef 1\n";
  print OUT "Cycles 10\n";
  print OUT "RecAnalyzedFile $RecAnalyzedFile\n";
  print OUT "ReceptorNMs $ReceptorNMs\n";
  print OUT "Lambda $lambda\n";
}
if ($#ARGV >15){
  print OUT "LigAnalyzedFile $LigAnalyzedFile\n";
  print OUT "LigandNMs $LigandNMs\n";
}

#print OUT "#ReceptorAlignment <RecAlignFile>\n";
#print OUT "#LigandAlignment <LigAlignFile>\n";

print OUT "## LOG-FILE params\n";
print OUT "log-file fiberdock.log\n";
print OUT "log-level 0\n";

close (OUT);

