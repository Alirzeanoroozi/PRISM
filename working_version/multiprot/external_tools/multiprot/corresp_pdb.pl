#!/usr/bin/perl 

#Ver 1.9


use strict;

my $home="/export/home/aaytuna/research/multiprotv1.5/";

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
    print "Usage: corresp_pdb.pl <results.txt> <res_num>\n";
    exit(0);
}

open RES,$ARGV[0];

my $target_pdb="transTarget.pdb";
my $tmp_file="tmp_file_pdb.pdb";
my $tmp_file2="tmp_file_pdb2.pdb";
my $tmp_file3="tmp_file_pdb3.pdb";
#my $join_pdb="join.pdb";
#my $join_part_pdb="joinPart.pdb";
 
my $status=0; #for $target_pdb


my @mols=();
my @selected_mols=();
my @selected_mols_file_names=();
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
        ($tmp, $filename)=split(/\s*:\s*/,$_);
        $mols[$#mols+1]=$filename;
        next;
    }
    if( /^Solution Num :/){
      my $num;
      ($tmp, $num)=split(/\s*:\s*/,$_);
      
      if($num == $ARGV[1]){
	$b_proc=1;
      }else{
	$b_proc=0;
      }
    }

    if($b_proc==1 and /^Reference Molecule/){
	($tmp, $molid)=split(/\s*:\s*/,$_);
        
        $referenceMol=$mols[$molid];
                
        #system("cat $referenceMol > $join_pdb");
        #system("cat $referenceMol > $join_part_pdb");

       
        $selected_mols[$#selected_mols+1]=($#selected_mols+1)."$tmp_file2"; #$referenceMol;
        $match_list[$#match_list +1]={};#%empty_hash;

	#store molecule file name
	my @tmp_fn=split '/', $referenceMol;

	#print "$referenceMol:Split:",$tmp_fn[0],":", $tmp_fn[1], ":\n";

	$selected_mols_file_names[$#selected_mols_file_names+1]=$tmp_fn[$#tmp_fn];
		
	`cp -f $referenceMol $#selected_mols$tmp_file2`;
	next;
    }

    if($b_proc==1 and /^Molecule/){
	@tmp=split(/\s*:\s*/,$_);
        $filename=$mols[ ($tmp[1]+0)];

        $selected_mols[$#selected_mols+1]=($#selected_mols+1)."$tmp_file2";#$filename;
        $match_list[$#match_list +1]={};#%empty_hash;

	#store molecule file name
	my @tmp_fn=split '/', $filename;
	$selected_mols_file_names[$#selected_mols_file_names+1]=$tmp_fn[$#tmp_fn];
        
	next;
    }

    if($b_proc==1 and /^Trans :/){
     	my $tline=$_;
     
        @tmp=split(/\s*:\s*/,$tline);
        #@tmp=split(' ',$tmp[3]);
        my $params=$tmp[1]; # 3/1
        #chop($params);
        #$params=$tmp[0]." ".$tmp[1]." ".$tmp[2]." ".$tmp[3]." ".$tmp[4]." ".$tmp[5];
        #system("$home/utils/getPDBheader.pl $filename > $tmp_file2");
        
        system("cat $filename | $home/utils/pdb_trans_all_atoms$suffix $params > $#selected_mols$tmp_file2");
        
        #system("mv $tmp_file $join_pdb");
        #`rm -f $tmp_file2`;


        next;
    }
    if($b_proc==1 and /^Match List/){
       
        $boolStartReadMatchList=1;
        next;
    }
    if($b_proc==1 and /^End of Match List/){
      $boolStartReadMatchList=0;
      next;
    }
    if($b_proc==1 and $boolStartReadMatchList==1){
      @tmp=split(' ',$_);
        
        ++$counter;
        for(my $index=0;$index<= $#tmp; $index++){
	    my @match_entry=split('\.',$tmp[$index]);
	    my $key=$match_entry[0].".".$match_entry[2];
	    $match_list[$index]->{ $key } = 1;
        }
        next;
    }
    
}
 
print "There are $counter k-tuples.\n";

$tmp=$#selected_mols+1;
print "Creating correspondence for $tmp molecules.\n";



open OUT_PDB_joinFull, ">joinFull.pdb";
open OUT_PDB_RASMOL_SCRIPT, ">script_joinFull.rsm";
print OUT_PDB_RASMOL_SCRIPT "load joinFull.pdb\n define mcore ";

open OUT_PDB_join, ">join.pdb";
open OUT_PDB_joinPart, ">joinPart.pdb";


$counter=0;
my $key="";

my $cRes=-10000000;

for(my $index=0;$index<= $#selected_mols; $index++){
  open PDB_FILE, $selected_mols[$index];
  
  #open OUT_PDB, ">$tmp_file";
  
  print OUT_PDB_joinFull "MODEL     ";
  printf OUT_PDB_joinFull "%4d   $selected_mols_file_names[$index]\n",$index+1;
  
  print OUT_PDB_join "MODEL     ";
  printf OUT_PDB_join "%4d   $selected_mols_file_names[$index]\n",$index+1;
  
  print OUT_PDB_joinPart "MODEL     ";
  printf OUT_PDB_joinPart "%4d   $selected_mols_file_names[$index]\n",$index+1;
  
  while(<PDB_FILE>){
    if(/^ATOM/ or /^SIGATM/ or /^ANISOU/ or /^SIGUIJ/ or /^TER/ or /^HETATM/){
      
      my $res=0+substr($_,22,4);
      my $chainid=substr($_,21,1);
    
      if($chainid =~ / /){
	$key=".".$res;
      }else{
	$key=$chainid.".".$res;
      }
      
      if($match_list[$index]{ $key } == 1){
	if($index!=0){
	  print OUT_PDB_join $_;
	}
	print OUT_PDB_joinPart $_;
	
	if($cRes!=$res){
	  $cRes=$res;
	  
	  if($counter>11){
	    $counter=0;
	    print OUT_PDB_RASMOL_SCRIPT 
	      "\n define mcore mcore or ";
	  }
	  my $tmp_chainid=($chainid eq ' ')?'':$chainid;
	  print OUT_PDB_RASMOL_SCRIPT 
	    "(:$tmp_chainid:",$index+1," and $res) or ";
	  $counter++;
	}
      }

      if($index==0){
	print OUT_PDB_join $_;
      }
      print OUT_PDB_joinFull $_;
    }
  }

  print OUT_PDB_joinFull "ENDMDL \n";
  print OUT_PDB_join     "ENDMDL \n";
  print OUT_PDB_joinPart "ENDMDL \n";
  
  #close OUT_PDB;
  close PDB_FILE;
  
  #if($index==0){
  #  `cp -f $referenceMol $join_pdb`;
  #  
  #  `mv $tmp_file $join_part_pdb`; 
  #}else{
  #  
  #  `$home/utils/joinXpdb.pl $join_pdb $tmp_file > $tmp_file3`;
  #  `mv -f $tmp_file3 $join_pdb`;
  #  
  #  `$home/utils/joinXpdb.pl $join_part_pdb $tmp_file > $tmp_file3`;
  #  `mv -f $tmp_file3 $join_part_pdb`;
  #}               
  
  print "",($index+1),": $selected_mols_file_names[$index]\n";


  if( -r $selected_mols[$index]){
    `rm -f $selected_mols[$index]`;
  }
  if( -r $tmp_file){
    `rm -f $tmp_file`;
  }
}

print OUT_PDB_RASMOL_SCRIPT  "\n";

