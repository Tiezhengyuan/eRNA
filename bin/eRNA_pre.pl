#! /usr/bin/perl

use strict;
use warnings;
use List::Util;
use List::MoreUtils;

#get subroutines
require "func_data.pm";
require "func_gui.pm";
require "func_basic.pm";
require "func_common.pm";
require "func_rna.pm";
require "func_mrna.pm";
require "func_ncrna.pm";

#get the directory of exe scripts involved in
my $dir_home=Cwd::abs_path(Cwd::getcwd().'/..' );
#store all parameters
my %variables=('dir_home'=>$dir_home.'/', 'dir_bin'=>$dir_home.'/bin/', 'dir_result'=>$dir_home.'/result/', 
				'dir_raw_data'=>$dir_home.'/raw_data/', 'specie'=>'human');

######################################

#$variables{dir_bowtie}='/home/yuan/eRNA/bowtie1/';
#$variables{software_aligner}='bowtie1';
#$variables{'index_name'}='hs_ref_GRCh38.p7_dna';

##################################### 
#1 directory of results
$variables{dir_result}=sub_common::stdin_dir($variables{'dir_result'}, 'Enter a directory for storing results');
$variables{file_var}=$variables{dir_result}.'variables.txt';
$variables{file_sample_info}=$variables{dir_result}.'sample_info.csv';

######################################
#2 directory of raw_data
while(1){
	#select parameters
	$variables{dir_raw_data}=sub_common::stdin_dir($variables{dir_raw_data}, 'Enter a directory for storing raw data in FASTQ');
	my @a=('single-end','paired-end');
	$variables{sequencing_end}=sub_common::stdin_select(\@a, 'Sequencing is single-end or double-end');
	
	#judge if exists raw data
	my $fq_files_pointer=sub_basic::SM_rawfiles_R1($variables{'dir_raw_data'});
	my @fq_files=@$fq_files_pointer;
	#print "@fq_files\n";
	if(@fq_files>0){#2
		#get hash of raw file name ~ sample name
		my $tag='NO';
		if(-f $variables{file_sample_info}){#3
			@a=('YES','NO');
			$tag=sub_common::stdin_select(\@a, 'Do you want RESERVE the sample_info.csv');
		}
		if($tag eq "NO"){#get fastq~sample_names
			my $rawfile_sample_pointer=sub_basic::SM_sample_from_raw($fq_files_pointer);
			$variables{rawfile_sample_pointer}=$rawfile_sample_pointer;
			printf("generate a new sample_info.csv: %s\n", $variables{file_sample_info});
			sub_common::hash_to_file($rawfile_sample_pointer, $variables{file_sample_info}, ',');
		}
		last;
	}#2
	else{
		print "\nError: No fastq files detected! Entry the rawdata directory again!\n";
	}
	#print "$variables{dir_raw_data}\n";
}


######################################
#3: select a pipeline, default is 'mrna';
my @pipelines=('mRNA', 'ncRNA','other');
$variables{'pipeline'}=sub_common::stdin_select(\@pipelines, 'Select a pipeline');
if ($variables{'pipeline'} eq 'ncRNA'){
	my $variables_pointer=sub_ncrna::initiate_ncRNA_variables(\%variables);
	%variables=%$variables_pointer;
}elsif ($variables{'pipeline'} eq 'mRNA'){
	my $variables_pointer=sub_mrna::initiate_mRNA_variables(\%variables);
	%variables=%$variables_pointer;
}

######################################
#4 specie
my @species=('human', 'mouse', 'rat', 'maize', 'yeast', 'none');
$variables{'specie'}=sub_common::stdin_select(\@species, 'Select specie');

######################################
print "Check bowtie index.\n";
my $file_names_pointer=sub_common::files_list($variables{'dir_bowtie'}, 'file_name');
#reference in fasta
my @fa_names=grep {$_=~/\.fa$|\.fasta$/} @$file_names_pointer;
sub_basic::check_bowtie_index($variables{dir_bowtie}, $variables{software_aligner});

##################################### 
#select files for mRNA-seq only
if ($variables{pipeline} eq 'mRNA'){
	$variables{fa_name}=sub_common::stdin_select(\@fa_names, 'Select reference sequences in FASTA');
	$variables{file_fa}=$variables{'dir_bowtie'}.$variables{'fa_name'};
	$variables{index_name}=sub_common::file_operation($variables{'fa_name'}, 'name_head');
	$variables{bowtie_index}=$variables{dir_bowtie}.$variables{index_name};
	$variables{trans_index}=$variables{dir_trans_index}.$variables{index_name};
	
	#select genome annotation in GTF/GFF3
	my @gtf_names=grep {$_=~/\.gtf$|\.gff3$|\.gff$/i} @$file_names_pointer;
	my @matched_gtf=grep {$_=~/$variables{'index_name'}/} @gtf_names;
	if(@matched_gtf==0){ #no gtf file detected.
		printf ("\nWarning: No gtf or gff file is detected in %s!\n", $variables{'dir_bowtie'});
	}else{
		$variables{'annot_name'}=sub_common::stdin_select(\@matched_gtf, 'Select genome annotation in GTF/GFF');
		$variables{'file_annot'}=$variables{'dir_bowtie'}.$variables{'annot_name'};
	}
	
	#create transcriptome index
	if(-f $variables{'file_annot'}){
		my $tophat=join('  ', $variables{exe_mapper}, '--transcriptome-index',  $variables{trans_index}, 
						'-G', $variables{file_annot}, $variables{bowtie_index}, );
		printf("Create transcriptome index: %s\n", $tophat);
		system("$tophat");
	}
	
	#tophat parameters
	my %lib_types=('First strand only'=>'fr-firststrand',  'Both strands'=>'fr-unstranded',  'Second strand only '=>'fr-secondstrand');
	my @names=sort keys %lib_types;
	my $type_name=sub_common::stdin_select(\@names, 'Select the type of the sequencing library');
	#export the tophat parameters formula
	$variables{tophat_options}=join('  ', $variables{exe_mapper}, '--library-type', $lib_types{$type_name}, 
									'--transcriptome-index', $variables{trans_index}, ); 
	$variables{cufflinks_options}=join('  ', $variables{exe_assembler}, '--library-type', $lib_types{$type_name}, '-G', $variables{file_annot}, );
}

######################################
#parameters only for ncRNA pipeline
if ($variables{'pipeline'} eq 'ncRNA'){
	#adapters trimming
	my @adapters_3=('TGGAATTCTCGGGTGCCAAGG', 'AGATCGGAAGAGCACACGTCT', 'NA');# /Illumina, NEB/
	my @adapters_5=('GTTCAGAGTTCTACAGTCCGACGATC', 'GTTCAGAGTTCTACAGTCCGACGATC', 'NA'); # /Illumina, NEB/
	$variables{'adapter_3'}=sub_common::stdin_select(\@adapters_3, 'Trim 3-end adapter sequences');
	$variables{'adapter_5'}=sub_common::stdin_select(\@adapters_5, 'Trim 5-end adapter sequences');
	
	#indexes for separate alignment
	my @index_names=sort map {$_=~s/\.fa$|\.fasta$//; $_} @fa_names;
	$variables{'index_seperate'}=sub_common::stdin_select(\@index_names, 'Select one or more bowtie index for seperate alignment');
	#indexes for interative alignment
	$variables{'index_iterative'}=sub_common::stdin_select(\@index_names, 'Select one or more bowtie index for iterative alignment');
	
	####set parameters of Bowtie
	#set mismatching
	my @mis=qw/1 2 3/;
	my $mismatch=sub_common::stdin_select(\@mis, 'Report alignments with at most mismatches');
	#suppress multiple alignments
	my @multiple=qw/1 2 3 4 5/;
	my $multi=sub_common::stdin_select(\@multiple, 'Suppress all valid alignments');
	#library type
	my %lib_types=('First strand only'=>'--norc',  'Both strands'=>' ',  'Second strand only '=>'--nofw');
	my @names=sort keys %lib_types;
	my $type_name=sub_common::stdin_select(\@names, 'Select the type of the sequencing library');
	#
	$variables{bowtie_options}=join('  ',  '-a --best --strata', '-v', $mismatch, $lib_types{$type_name}, '-m', $multi)
}

######################################
#multi-threading
my @threads=qw/1 2 4 8 16 24/;
$variables{'threads_num'}=sub_common::stdin_select(\@threads, 'Number of multiple threads:');

######################################
######################################

######################################

#export var file
sub_common::hash_to_file(\%variables, $variables{file_var}, '=');
print sub_data::print_hash(\%variables);
print "\n\nThe setup is done!\n\n";