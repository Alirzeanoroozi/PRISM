#!/usr/bin/perl -w

#Ver 1.5


use strict;

my $home="/home/silly6/mol/demos/MultiProt/";

#detect the machine type
my $mtype=`uname -s`;
my $suffix;
if($mtype =~ /IRIX/){
    $suffix= ".IRIX";
}else{
    if($mtype =~ /Linux/){
        $suffix= ".Linux";
    }else{
        print "This type of machine: $mtype";
        print "is not supported.\n";
        print "The supported OS's are: Linux, IRIX and IRIX64\n";
        print "You can check the type of the system by: uname -s\n";
        exit(0);
    }
}


if( $#ARGV != 1) {
    print "Usage: trans_mult.pl <results.txt> <res_num>\n";
    exit(0);
}

open RES,$ARGV[0];

 
my @mols=();

my @selected_mols=();


my @match_list=();

my $filename;
my $referenceMol;
my $counter=0;
my $b_proc=0;
my $boolStartReadMatchList=0;
my @tmp;
my $tmp;
my $molid;

while(<RES>){
    chomp;
    if(/^Mol-/){
        ($tmp, $filename)=split(':',$_);
        $mols[$#mols+1]=$filename;
        next;
    }
    if( /^Solution Num :/){
	if(/^Solution Num : $ARGV[1]/){
	    $b_proc=1;
	}else{
	    $b_proc=0;
	} 
    }

    if($b_proc==1 and /^Molecule/){
	@tmp=split(':',$_);
        $filename=$mols[ ($tmp[1]+0)];

      	next;
    }

    if($b_proc==1 and /^Trans :/){
     	my $tline=$_;
     
        @tmp=split(':',$tline);
	my $params=$tmp[1]; 
	
	
#$params=$tmp[0]." ".$tmp[1]." ".$tmp[2]." ".$tmp[3]." ".$tmp[4]." ".$tmp[5];
        #system("$home/utils/getPDBheader.pl $filename > $tmp_file2");
        
	$tline=$filename.".trans.pdb";
        system("cat $filename | $home/utils/pdb_trans_all_atoms$suffix $params > $tline");
        print "Creating file $tline\n";
	next;
    }
}
