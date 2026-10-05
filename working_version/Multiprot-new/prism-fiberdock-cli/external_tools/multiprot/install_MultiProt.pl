#!/usr/bin/perl

#installes MultiProt


chomp( $pwd=`pwd`);

#################################################
$text_lib=" \
#MultiProt start - don't remove this line \
set path=( $pwd \$path ) \
#MultiProt end - don't remove this line";




@text=`cat ~/.cshrc`;

$fname="~/.cshrc_before_changed_by_MultiProt";

`mv ~/.cshrc $fname`;

open RESORS, ">tmp.output";

$state=0;
foreach (@text){
  if(/#MultiProt start/){
    $state=1;
    next;
  }
  if(/#MultiProt end/){
    $state=0;
    next;
  }
  if($state==0){
    print RESORS $_;
  }
}
print RESORS $text_lib;

close RESORS;

`mv tmp.output ~/.cshrc`;


##############################################


open CORRESP_READ, "corresp_pdb.pl";
@text=<CORRESP_READ>;
close CORRESP_READ;

open CORRESP, ">corresp_pdb.pl";

foreach (@text){
  if(/my \$home=/){
    print CORRESP "my \$home=\"${pwd}/\";\n";
    next;
  }
  print CORRESP $_;
}

close CORRESP;


print "MultiProt was successfully installed.\n";

print "\nRun: source ~/.cshrc  \n";
print "or open a new terminal.\n";

print "\nYou can run MultiProt from any place:.\n";
print "multiprot.Linux pdb1 pdb2 ...\n";

print "\nFor further details see README.\n";
