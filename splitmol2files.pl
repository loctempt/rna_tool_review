#!/usr/bin/perl
##############################################################################
### SECTION:3: splits multi-MOLECULE input mol2 into individual mol2 files ###
##############################################################################

print("Splitting the multiple Ligand mol2 file into single mol2 files.\n");
#print("The single mol2 files will have the prefix: $prefix \n");
my($mol2_file) = "";
my($global_count) = 0;
my($local_count) = 0;

$mol2_file = $ARGV[0];
open (MOL2, $mol2_file)
    || die "\nError, cannot open ligand mol2 file: $mol2_file\n";

# assign basename for output files based on input MOL2 filename:

while (<MOL2>)
{
    chomp($_);
    if ($_ =~ /^(@<TRIPOS>MOLECULE)/)
    {
        $global_count++;
        open (OUT, ">$global_count.mol2")
            || die "\nError, cannot open output mol2 file: $global_count.mol2\n";
        print OUT "$_\n";
        $local_count++;
    }

    if ($_ =~ /^(@<TRIPOS>MOLECULE)/ && $local_count > 1)
    {
        $local_count = 1;
        next;
    }

    if ($_ !~ /^(@<TRIPOS>MOLECULE)/ && $local_count == 1 )
    {
        print OUT "$_\n" ;
        next;
    }
}
print("The number of single Ligand mol2 files generated is $global_count.\n");
close OUT;

