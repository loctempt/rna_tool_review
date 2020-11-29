#! /usr/bin/perl
#using chimera to supperpose structures:
open(CMDI, ">superimpose.cmd") || die("unable to create file: supperpose");
print CMDI "open template.pdb\n";
print CMDI "open align.pdb\n";
print CMDI "mm #0 #1\n";
print CMDI "write #1 /home/leicq/recommendation_algo/data/superimpose/align_tmp.pdb\n";
close(CMDI);

 `/home/leicq/.local/UCSF-Chimera64-1.14/bin/chimera --nogui  superimpose.cmd`;



