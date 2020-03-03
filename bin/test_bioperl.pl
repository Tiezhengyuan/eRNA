#! /usr/bin/perl -w
use strict;
use warnings;
use List::Util;
use File::Find;


############
#import personal modules
use func_basic;
use func_common;
use func_bioperl;
use func_bioseq;
use func_data;


############
#human transcript id
#sub_bioperl::EntrezGene_query('NM_000014');
#sub_bioperl::EntrezGene_query('NM_007318');
#sub_bioperl::EntrezGene_query('NM_003988');
#human transcript id (mRNA)
#sub_bioperl::GenBank_query('NM_000014');
#virus accession number (DNA)
#sub_bioperl::GenBank_query('NC_000858');

#sub_bioperl::GenBank_query('NP_049561');

#swissprot
#virus
sub_bioperl::query_SwissProt('A0A126');
#human
#sub_bioperl::query_SwissProt('Q44238');
#sub_bioperl::query_SwissProt('NM_003980');

#test taxon_id
my $d='a';
#eval{ 	$d= sub_bioperl::taxon_id('52275', 'domain');	};
#print "=$d=\n";
#sub_bioperl::taxon_tree('145856');



my %par=('program'=>'blastn', 'expect'=>0.01, 'identity_percent'=>90, 'database'=>'/home/yuan/phip/blast/VirScan_v1', 
			'infile'=>'/home/yuan/phip/ref_seq/test.fa', 'outfile'=>'/home/yuan/phip/blast/test.bls',);
#sub_bioperl::local_blast(\%par);

###extract annotations from NCBI
#my @species=qw/human mouse rat cattle dog chicken maize zebrafish soybean rice/;
my @species=qw/rice/;
foreach my $specie(@species){
	#sub_bioperl::download_NCBI_genome('DNA_fa', $specie, '/home/yuan/data_preparation/');
	#sub_bioperl::download_NCBI_genome('RNA_fa', $specie, '/home/yuan/data_preparation/');
	#sub_bioperl::download_NCBI_genome('RNA_gbk', $specie, '/home/yuan/data_preparation/');
	#sub_bioperl::download_NCBI_genome('protein_fa', $specie, '/home/yuan/data_preparation/');
	#sub_bioperl::download_NCBI_genome('protein_gbk', $specie, '/home/yuan/data_preparation/');
	#sub_bioperl::download_NCBI_genome('annot_gff', $specie, '/home/yuan/data_preparation/');
}
#Ensemble
@species=qw/human mouse rat cattle dog chicken/;
#@species=qw/rat/;
foreach my $specie(@species){
	#sub_bioperl::download_Ensemble_genome('DNA_fa', $specie, '/home/yuan/data_preparation/');
	#sub_bioperl::download_Ensemble_genome('RNA_fa', $specie, '/home/yuan/data_preparation/');
	#sub_bioperl::download_Ensemble_genome('protein_fa', $specie, '/home/yuan/data_preparation/');
	#sub_bioperl::download_Ensemble_genome('annot_gff', $specie, '/home/yuan/data_preparation/');
}


######
my $pro_gbk_file='/home/yuan/data_preparation/NCBI/human/protein.gbk';
my $out_file='/home/yuan/data_preparation/NCBI/human/protein.txt';
#sub_bioperl::protein_gbk($pro_gbk_file, $out_file);

#my $gbk_file='/home/yuan/data_preparation/NCBI/human/rna.gbk';
#my $specie='human_NCBI';
#my $gbk_file='/home/yuan/data_preparation/NCBI/maize/rna.gbk';
#my $specie='maize_NCBI';
#my $result_dir=sub_common::file_operation($gbk_file, 'directory').'/RNA_annot';
#sub_bioperl::gbk_to_fa($gbk_file, $result_dir, $specie);
print "ok\n";