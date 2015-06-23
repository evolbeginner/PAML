#! /bin/env perl

use Getopt::Long;
use File::Basename;
use File::Copy;
use File::Temp qw/ tempfile tempdir /;
use Cwd;
use strict;
use 5.010;

my $yn00_ctl="$ENV{'HOME'}/software/phylo/paml4.7/yn00.ctl";
my $codon_alignment;

###############################################################
GetOptions(
    'i|in=s'   =>  \$codon_alignment,
) || die "illegal params!";

my $tempdir = tempdir( CLEANUP => 1 );
my $new_yn00_ctl="$tempdir/yn00.ctl";
my $MLmatrix;

###############################################################
my $basename_codon_alignment = basename($codon_alignment);
open(my $IN, '<', $yn00_ctl) || die "yn00_ctl $yn00_ctl could not be opened!";
open(my $OUT,'>', $new_yn00_ctl) || die "new_yn00_ctl $new_yn00_ctl could not be created!";
while(my $line=<$IN>){
    chomp($line);
    if ($line =~ /^(\s+seqfile\s+\=\s+)(.+)/){
        $line = $1.$basename_codon_alignment;
    }
    print $OUT $line."\n";
}
close $IN;
close $OUT;

my $cwd=getcwd;
copy($codon_alignment,"$tempdir/$basename_codon_alignment") || die "copy failed!";
chdir($tempdir);
`yn00`;
chdir($cwd);

###############################################
open(my $IN,'<',"$tempdir/yn") || die "cannot be opened!";
while(<$IN>){
    last if (/^seq\.\s+seq\./);
}
<$IN>;
while(<$IN>){
    if (
        m/^\s+(\d+)\s+  # seq #
        (\d+)\s+        # seq #
        (\d+(\.\d+))\s+ # S
        (\d+(\.\d+))\s+ # N
        (\d+(\.\d+))\s+ # t
        (\d+(\.\d+))\s+ # kappa
        (\d+(\.\d+))\s+ # omega
        \-??(\d+(\.\d+))\s+ # dN
        \+\-\s+
        \-??(\d+(\.\d+))\s+ # dN SE
        \-??(\d+(\.\d+))\s+ # dS
        \+\-\s+
        \-??(\d+(\.\d+))\s+ # dS SE
        /ox
      )
        {
            $MLmatrix = {
                'S'     => $3,
                'N'     => $5,
                't'     => $7,
                'kappa' => $9,
                'omega' => $11,
                'dN'    => $13,
                'dN_SE' => $15,
                'dS'    => $17,
                'dS_SE' => $19,
            };
        }
    last;
}

print join ("\t", map {$MLmatrix->{$_}} qw(dN dS omega))."\n";


