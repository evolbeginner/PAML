#! /usr/bin/perl

=cut
1.	bioperl_yn00.pl or run_yn00.pl is required!
=cut

use 5.010;
use Getopt::Long;
use Bio::SeqIO;
use Cwd;

local (%pep_seq,$CDS_seq,%pair,%para);
my ($temp_dir,$force_align,$force_pal2nal);
my ($count_pair);
my $cwd = getcwd;
my $bioperl_yn00="bioperl_yn00.pl";
my $run_yn00="run_yn00.pl";
my $separator='_';

#########################################################################################
($pep_file, $CDS_file, $duplicate_list, $outfile, $force_align, $force_pal2nal, $force_both, $seq_title_RegExp, $seq_title_RegExp_del, $separator) = &get_param();
($force_align, $force_pal2nal) = (1,1) if defined $force_both;

&check_param(qw(pep_file CDS_file duplicate_list outfile));

foreach (qw(pep_file CDS_file duplicate_list outfile)){
	print "$_ is\t${$_}\n";
}

open(my $OUTFILE, '>' , $outfile) || die "OUTFILE $outfile cannot be opened!\n";

*pep_seq = &read_fasta($pep_file,$seq_title_RegExp,$seq_title_RegExp_del);
*CDS_seq = &read_fasta($CDS_file,$seq_title_RegExp,$seq_title_RegExp_del);

(*pair,*para) = &generate_duplicate_rela($duplicate_list);

($temp_dir) = &create_temp_folder($duplicate_list);

for my $pair_name(sort keys %pair){
	&yn00($pair_name) if $yn00;
	&seq_similarity($pair_name) if $seq_similarity;
}

#system "rm -rf \"$temp_dir\"";

##########################################################################################
sub generate_duplicate_rela
{
my ($duplicate_list) = @_;
my (%pair,%para);
open(my $IN,'<',$duplicate_list) or die "cannot open duplicate_list file:$!";
while(<$IN>){
	chomp;
	@line=split;
	$pair{join $separator,sort @line[0,1]}=1;
	foreach(0..1){
		$para{$line[$_]}=$line[1-$_];
	}
}
close $IN;
return (\%pair,\%para);
}

#----------------------------------------------#
########	sub-line function	########
sub create_temp_folder{	
	my ($temp_dir,$basename,$dirname);
	my ($ori_name)=@_;
	$basename = `basename "$ori_name"`;
	$dirname  = `dirname "$ori_name"`;
	chomp($basename);
	chomp($dirname);
	$temp_dir="$dirname/$basename".'_temp';
	system "mkdir \"$temp_dir\"";
	#system "[ ! -e \"$temp_dir\" ] && mkdir \"$temp_dir\"" or last;
	return ("$temp_dir");
}

#------------------------------------------------------------------------------------#
sub get_param
{
GetOptions(
'yn00!'			=>	\$yn00,
'seq_similarity!'	=>	\$seq_similarity,
'pep=s' 		=>	\$pep_file,
'CDS=s' 		=>	\$CDS_file,
'dup_list=s'		=>	\$duplicate_list,
'force_align!'		=>	\$force_align,
'force_pal2nal!'	=>	\$force_pal2nal,
'force_both!'		=>	\$force_both,
'out=s'			=>	\$outfile,
'seq_title_RegExp=s'    =>      \$seq_title_RegExp,
'seq_title_RegExp_del=s'=>      \$seq_title_RegExp_del,
'sep|separator=s'           =>      \$separator,
) || die "illegal params!";
return($pep_file,$CDS_file,$duplicate_list,$outfile,$force_align,$force_pal2nal,$force_both,$seq_title_RegExp,$seq_title_RegExp_del,$separator);
}


sub read_fasta{
        my ($input_seq,$seq_title_RegExp,$seq_title_RegExp_del) = @_;
        my %seq;
        my $catchseq_seqio_obj;
        $catchseq_seqio_obj = Bio::SeqIO->new(-file=>$input_seq, -format=>'fasta');
        while(my $seq_obj=$catchseq_seqio_obj->next_seq()){
                my ($full_name);
                $full_name = $seq_obj->display_name;
                do {$full_name = $1 if $full_name =~ m/($seq_title_RegExp)/} if $seq_title_RegExp;
                $full_name =~ s/($seq_title_RegExp_del)// if $seq_title_RegExp_del;
                $seq{$full_name} = $seq_obj->seq;
        }
        return(\%seq);
}


sub check_param{
	if(not $yn00 and not $seq_similarity){
		die "yn00 or seq_similarity should be specified!\n";
	}
	if ($yn00){
		die "Arguments for yn00 have not been specified enough!\n" if ((not $CDS_file) or (not $pep_file) or (not $duplicate_list));
	}
	if ($seq_similarity){
		die "Arguments for seq_similarity have not been specified enough!\n" if (not $pep_file or not $dup_list);
	}
}


sub yn00{
	my ($pair_name) = @_;
	my $pal2nal_err=0;
	my @pair=split /$separator/,$pair_name;
	my ($pep_file,$CDS_file,$pep_aln_file,$codon_aln_file) = ("$temp_dir/$pair_name.pep","$temp_dir/$pair_name.cds","$temp_dir/$pair_name.pep.aln","$temp_dir/$pair_name.codon_aln");
	open($pep_OUT, '>', $pep_file);
	open($pep_OUT,'>>', $pep_file);
	open($CDS_OUT, '>', $CDS_file);
	open($CDS_OUT,'>>', $CDS_file);
	for (@pair){
		print $pep_OUT ">$_"."\n".$pep_seq{$_}."\n";
		print $CDS_OUT ">$_"."\n".$CDS_seq{$_}."\n";
	}
	close $CDS_OUT;
	close $pep_OUT;
	
	system "muscle -in \"$pep_file\" -out \"$pep_aln_file\" 1>/dev/null" if ($force_align or ! -e $pep_aln_file);
	system "pal2nal.pl \"$pep_aln_file\" \"$CDS_file\" -output fasta > \"$codon_aln_file\" ||  muscle_paml.pl \"$CDS_file\" > \"$codon_aln_file\"" if ($force_pal2nal or ! -e $codon_aln_file);
	print $OUTFILE $pair_name."\t";
	#my $haha=`$bioperl_yn00 $codon_aln_file`;
	my $haha=`$run_yn00 --in \"$codon_aln_file\"`;
	chomp($haha);
	print $OUTFILE "$haha\n";
	#print $count_pair."\n" if ++$count_pair % 50 == 0;
}


