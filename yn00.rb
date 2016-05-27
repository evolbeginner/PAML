#! /bin/env ruby


# The script is to calculate Ka, Ks and Ka/Ks for gene pairs.
# PAML, Perl, Ruby, Bioruby, MUSCLE and pal2nal.pl are needed. Please type 'ruby yn00.rb -h' to see the usage.
# The script is written by Sishuo Wang from University of British Columbia and is free for anyone.
# Please write e-mail to sishuowang@hotmail.com if you have any questions or want to report bugs. Your help will be highly appreciated.


#################################################################
require "getoptlong"
require "bio"


#################################################################
#Please specify the path to the following scripts or software
pal2nal = "pal2nal.pl"
muscle = "muscle"

dirname=File.dirname($0)
run_yn00 = File.join(dirname, "run_yn00.pl")

# You might also need yn00 installed and can be found in the default searching PATH.


#################################################################
def read_seq(seq_file)
  seq_objs = Hash.new
  Bio::FlatFile.open(seq_file, "r").each_entry do |f|
    seq_objs[f.definition] = f
  end
  return(seq_objs)
end


def read_pairs(pair_file, sep, is_sort=false)
  pairs = Hash.new
  File.open(pair_file, "r").each_line do |line|
    line.chomp!
    line_arr = line.split(sep)
    genes = is_sort ? line_arr.sort : line_arr
    raise "Not two genes are found in the infile line #{line}. Exiting ......" if genes.size != 2
    pairs[genes] = ""
  end
  return(pairs)
end


def output_fasta(genes, seq_objs, out_fasta)
  out_fh = File.open(out_fasta, 'w')
  genes.each do |gene|
    out_fh.puts ">" + gene
    out_fh.puts seq_objs[gene].seq
  end
  out_fh.close
end


def align(pairs, cds_seq_objs, pep_seq_objs, tmp_outdir, muscle, pal2nal)
  codon_alns = Array.new
  puts "Starting aligning sequences ......"
  pairs.each_pair do |genes, v|
    genes.each do |gene|
      raise "#{gene} is not found in cds file" if not cds_seq_objs.include?(gene)
      raise "#{gene} is not found in pep file" if not pep_seq_objs.include?(gene)
    end
    
    pair = genes.join("-")
    cds_fasta = File.join(tmp_outdir, pair+".cds.fasta")
    pep_fasta = File.join(tmp_outdir, pair+".pep.fasta")
    output_fasta(genes, cds_seq_objs, cds_fasta)
    output_fasta(genes, pep_seq_objs, pep_fasta)
    pep_aln = File.join(tmp_outdir, pair+".pep.aln")
    codon_aln = File.join(tmp_outdir, pair+".codon.aln")
    puts pair
    `#{muscle} -in #{pep_fasta} -out #{pep_aln} -quiet 2>/dev/null`
    `#{pal2nal} #{pep_aln} #{cds_fasta} -output paml > #{codon_aln}`
    codon_alns << codon_aln
  end
  return(codon_alns)
end


def run_yn00(run_yn00, codon_aln_files, outfile)
  out_fh = File.open(outfile, 'w') or raise "outfile #{outfile} cannot be created! Exiting ......"
  codon_aln_files.each do |codon_aln_file|
    basename = File.basename(codon_aln_file)
    basename =~ /(.+)\.codon\.aln$/
    pair = $1
    out_content = `perl #{run_yn00} --in #{codon_aln_file}`
    out_content.chomp!
    out_fh.puts [pair,out_content].join("\t")
  end
  out_fh.close
end


def help
  basename = File::basename($0)
  puts "Usage: ruby #{basename} <-o|--out> <-i|--in|-p|--pair infile> <--cds cds_file> <--pep pep_file> [options]"
  puts "Options:"
  puts "--sep\t\tseparator between gene pairs"
  puts "\t\tdefault\tTAB"
  puts "--tmp\t\ttmp_outdir. Alignments will written to this folder."
  puts "--sort\t\tsort gene names in the gene pair"
  puts "--force\t\tremove tmp_outdir if it exists."
  puts "\t\tdefault: disabled"
  puts "-h|--help\thelp information"
  puts
  puts "Please make sure that related scripts/software are available before the script is run."
  puts "Please specify the path the the scripts/software at the beginning of this script."
  puts "Please write e-mail to sishuowang@hotmail.com for any question/suggestions."
  puts
  exit
end


#################################################################
infile = nil
sep = "\t"
pep_file = nil
cds_file = nil
is_sort = false
tmp_outdir = "tmp"
is_force = false
outfile = nil


#################################################################
opts = GetoptLong.new(
  ['-i', '--in', "-p", "--pair", GetoptLong::REQUIRED_ARGUMENT],
  ['-o', '--out', GetoptLong::REQUIRED_ARGUMENT],
  ["--cds", GetoptLong::REQUIRED_ARGUMENT],
  ["--pep", GetoptLong::REQUIRED_ARGUMENT],
  ["--sep", GetoptLong::REQUIRED_ARGUMENT],
  ["--tmp", "--tmp_outdir", GetoptLong::REQUIRED_ARGUMENT],
  ["--sort", GetoptLong::NO_ARGUMENT],
  ["--force", GetoptLong::NO_ARGUMENT],
  ["-h", "--help", GetoptLong::NO_ARGUMENT],
)


help if ARGV.empty?

begin
  opts.each do |opt, value|
    case opt
      when '-i', '--in', '-p', '--pair'
        infile = value
      when '-o', '--out'
        outfile = value
      when '--cds'
        cds_file = value
      when '--pep'
        pep_file = value      
      when '--sep'
        sep = value
      when '--sort'
        is_sort = true
      when '--tmp', '--tmp_outdir'
        rmp_outdir = value
      when '--force'
        is_force = true
      when '-h', '--help'
        puts "Illegal arguments! Exiting ......"
        help
    end
  end
rescue Exception
  help
end


if outfile.nil?
  puts "outfile must be specified!"
  help
elsif cds_file.nil?
  puts "cds must be specified!"
  help
elsif pep_file.nil?
  puts "pep must be specified!"
  help
elsif infile.nil?
  puts "pair_file must be specified!"
  help
end


if File.exists?(tmp_outdir)
  if is_force
    `rm -rf #{tmp_outdir}`
  else
    raise "tmp_oudir #{tmp_outdir} has already existed! Please specify a new tmp_outdir or use --force. Exiting ......"
  end
end
`mkdir -p #{tmp_outdir}`


#################################################################
cds_seq_objs = read_seq(cds_file)

pep_seq_objs = read_seq(pep_file)

pairs = read_pairs(infile, sep, is_sort)

codon_alns = align(pairs, cds_seq_objs, pep_seq_objs, tmp_outdir, muscle, pal2nal)

run_yn00(run_yn00, codon_alns, outfile)


