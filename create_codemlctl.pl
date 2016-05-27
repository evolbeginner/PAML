#! /usr/bin/perl
use 5.010;
use warnings;

my $codemlctl_file="/home/sswang/software/phylo/paml4.7/codeml.ctl";
my @parameter=($model,$NSsites,$seqfile,$treefile,$outfile,$ncatG,$seqtype,$clock);
my @parameter_name = qw(model NSsites seqfile treefile outfile ncatG seqtype output_codemlctl_dir clock);
my $force;

&read_parameter(@ARGV);

&check_AND_define_parameter(@parameter_name);

&create_codemlctl($codemlctl_file);

&ShuChu_parameter(@parameter_name);


########################################
########################################
sub read_parameter
{
my @para=@_;
&usage and exit if not defined $para[0];
foreach (@para){
	given ($_){
	when (/^\-model\=(.+)/)	{$model=$1}
	when (/^\-seqfile\=(.+)/)	{$seqfile=$1}
	when (/^\-treefile\=(.+)/)	{$treefile=$1}
	when (/^\-outfile\=(.+)/)	{$outfile=$1}
	when (/^\-NSsites\=(.+)/)	{$NSsites=$1; $NSsites=&modify_NSsites($1)}
	when (/^\-ncatG\=(.+)/)	{$ncatG=$1}
	when (/^\-seqtype\=(.+)/)	{$seqtype=$1}
	when (/^\-clock\=(.+)/)	{$clock=$1;}
	#########################
	when (/^\-output_codemlctl_dir\=(.+)/)	{$output_codemlctl_dir=$1}
	when (/^\-\-?force/)	{$force='YES'}
	}
	}
}

sub check_AND_define_parameter
{
my @para=@_;
foreach (@para){
	if (not defined ${$_}){
		given ($_){
		when ('model')	{$model=0}
		when ('NSsites')	{$NSsites='0 3'}
		when ('seqfile')	{die "no seqfile given\n"}
		when ('treefile')	{die "no treefile given\n"}
		when ('outfile')	{$outfile='mlc'}
		when ('seqtype')	{$seqtype=1}
		when ('ncatG')	{$ncatG=3}
		when ('clock')	{$clock=0}
		when ('output_codemlctl_dir')	{$output_codemlctl_dir='.'; &check_output_codemlctl_dir('.')}
		}
		}
	else{
		given($_){
		when ('model')	{die "model parameter illegal" if $model !~ /^[012]$/}
		when ('NSsites')	{my $return_value=&check_NSsites($NSsites); die "NSsites parameter illegal\n" if $return_value eq 'die'}
		when ('seqfile')	{}
		when ('treefile')	{}
		when ('outfile')	{}
		when ('seqtype')	{die "seqtype parameter illegal\n" if $seqtype !~ /^[123]$/}
		when ('ncatG')		{die 'ncatG parameter illegal\n' if $ncatG !~ /^\d+$/}
		when ('output_codemlctl_dir')	{&check_output_codemlctl_dir($_)}
		}
		}
	}
}


sub check_NSsites
{
	my $return_value=0;
	my ($NSsites)=@_;
	if ($NSsites !~ /^\d+( \d+)*$/)	{return ('die')}
	else{ 
		my @a=split / /,$NSsites; 
		foreach(@a) {$return_value='die' if $_ !~ /^[0-9]$/};
		return ($return_value);
	}
}


sub create_codemlctl
{
my ($codemctl_file)=@_;
open(OUT,"> $output_codemlctl_dir/codeml.ctl");
open(IN,$codemctl_file);
while(my $line=<IN>)
{
	chomp($line);
	given($line){
	when(/^\s+seqfile/)	{$line=~s/(?<=seqfile \=)[^\*]+(?=\*)/ $seqfile /;}
	when(/^\s+treefile/)	{$line=~s/(?<=treefile \=)[^\*]+(?=\*)/ $treefile /}
	when(/^\s+seqtype/)	{$line=~s/(?<=seqtype \=)[^\*]+(?=\*)/ $seqtype /;}
	when(/^\s+model/)	{$line=~s/\=(.+)/\= $model/}
	when(/^\s+NSsites/)	{$line=~s/(?<=NSsites \=)[^\*]+(?=\*)/ $NSsites /}
	when(/^\s+ncatG/)	{$line=~s/(?<=ncatG \=)[^\*]+(?=\*)/ $ncatG /}
	when(/^\s+clock/)	{$line=~s/(?<=clock \=)[^\*]+(?=\*)/ $clock /}
	}
	print OUT $line."\n";
}
close OUT;
close IN;
}


sub modify_NSsites
{
my ($value)=@_;
my @p;
foreach (split /\,/,$value){
	push @p,$_;
	}
my $return_value=join (' ',@p);
return ($return_value);
}


sub check_output_codemlctl_dir
{
my ($output_codemlctl_dir)=@_;
################################
return ('OK') if defined $force;
################################
if (-e "$output_codemlctl_dir/codeml.ctl"){
	say "output_codemlctl_dir already exists";
	print "Would you like to cover the old one?\t";
	system "print_color.pl '[yes or no]?' 'bold green'";
	my $char=<STDIN>;
		if ($char =~ /^yes$/i)	{$return_value='OK'}
		else	{$return_value='die'; die "cannot cover the old codeml.ctl file\n";}
	}
else	{$return_value='OK'}
return ($return_value);
}


sub ShuChu_parameter
{
	@parameter_name = @_;
	my $parameter_file="codemlctl_parameter.txt";
	my $OUT2;
	open($OUT2,">$output_codemlctl_dir/$parameter_file");
	foreach (@parameter_name){
		printf $OUT2 "%20s\t%20s\n",$_,${$_};
		}
	close $OUT2;
}


sub usage{
    my $usage="usage:\t-model -NSsites -seqfile -treefile -outfile -ncatG -seqtype -output_codemlctl_dir";
    system "print_color.pl \"$usage\" 'bold blue'";
}


