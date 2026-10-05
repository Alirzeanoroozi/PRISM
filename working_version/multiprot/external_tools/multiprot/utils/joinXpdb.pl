#!/usr/bin/perl

#ver 1.1

if($#ARGV<1){
    print("Usage:joinXpdb.pl pdbFile1 pdbFile2 ..\n");
    exit 0;
}

$chainID='A';

$currentChain='-';

foreach $file (@ARGV){
    if (! -r $file){
        print("File $file doesn't exist\n");
        exit(0);
    }
   

    open FILE, $file;
    while ($line=<FILE>){
        if($line =~ /^ATOM|HETATM/){
	    if (length ($chainID) >1){
		$chainID='A';
	    }
 
            if( !($currentChain =~ substr($line,21,1))){
                if(! ($currentChain =~ '-')){$chainID++;}
                $currentChain=" ".substr($line,21,1);
            }
            substr($line,21,1)=$chainID;
            
            print $line;
        }
    }
    $chainID++;
}
