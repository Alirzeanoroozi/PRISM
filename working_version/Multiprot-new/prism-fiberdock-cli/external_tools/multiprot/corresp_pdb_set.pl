#!/usr/bin/perl

#Ver 1.2


#$model="model.pdb";
#$target="target.pdb"

$home="/home/silly6/mol/demos/MultiProt/";

#detect the machine type
$mtype=`uname -s`;
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
    print "Usage: corresp_pdb_set.pl <results.txt> <res_num>\n";
    exit(0);
}

open RES,$ARGV[0];

$bool_do=0;   
$target_pdb="transTarget.pdb";
$tmp_file="tmp_file_pdb.pdb";
$tmp_file2="tmp_file_pdb2.pdb";
$join_pdb="join.pdb";
$join_part_pdb="joinPart.pdb";
$join_pdb_full="joinFull.pdb";


$status=0; #for $target_pdb


@mols=();

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

    if($b_proc==1 and /^Reference Molecule/){
        ($tmp, $molid)=split(':',$_);
        
        $referenceMol=$mols[$molid];
        
        system("cat $referenceMol > $join_pdb");
        system("cat $referenceMol > $join_part_pdb");
	system("cat $referenceMol > $join_pdb_full");
        next;
    }

    if($b_proc==1 and /^Molecule/){
	@tmp=split(':',$_);
        $filename=$mols[ ($tmp[1]+0)];

        next;
    }
    
    if($b_proc==1 and /^0-/){
        @tmp=split(':',$_);
        $params=$tmp[3];
        
        system("cat $filename | $home/utils/pdb_trans_frag$suffix $params > $tmp_file2");
	system("/home/silly6/mol/utils/joinXpdb/joinXpdb.pl $join_pdb_full $tmp_file2 > $tmp_file");
	system("mv $tmp_file $join_pdb_full");
                
        `rm -f $tmp_file2`;

	###############################################
	$params.=" ".$tmp[4];

	system("cat $filename | $home/utils/pdb_trans_frag$suffix $params > $tmp_file2");
        
        system("/home/silly6/mol/utils/joinXpdb/joinXpdb.pl $join_pdb $tmp_file2 > $tmp_file");
        
        system("mv $tmp_file $join_pdb");
                
        `rm -f $tmp_file2`;

        ###############################################
        system("cat $filename | $home/utils/pdb_trans_frag$suffix $params > $tmp_file2");
        
        system("/home/silly6/mol/utils/joinXpdb/joinXpdb.pl $join_part_pdb $tmp_file2 > $tmp_file");
        
        system("mv $tmp_file $join_part_pdb");
                
        `rm -f $tmp_file2`;

        next;
    }
}




