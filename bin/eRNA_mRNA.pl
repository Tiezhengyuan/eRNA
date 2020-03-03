#! /usr/bin/perl -w
use strict;
use warnings;
use threads; 
use Cwd;
use List::Util;
use List::MoreUtils;

#get subroutines 
use func_bioperl;
use func_mrna;
use func_rna;
use func_basic;
use func_common;
use func_data;


################################################################################################################
#main program
################################################################################################################
print "Main running begins:\n\n\n";


#get E_RNA variables from the file E_RNA_pipeline.info
#it should be saved in the same directory with E_RNA.pl
my $var_file=$ARGV[0];
my $variables_pointer=sub_common::file_to_hash($var_file, '=');
my %variables=%$variables_pointer;

#check directories
$variables_pointer=sub_rna::initiate_directories(\%variables);
%variables=%$variables_pointer;
sub_data::print_hash(\%variables);#print %variables

#associate sample name with *.fastq files list
$variables_pointer=sub_basic::SM_sample_info(\%variables);
%variables=%$variables_pointer;
my $sample_info_pointer=$variables{sample_info_pointer};
my %sample_info=%$sample_info_pointer;


print "initiate and clear  the monitor.log and system_monitor.log\n\n";
my $total_beginning_time=join(',', localtime(time) ); #beginning time record
sub_basic::initiate_log_files(\%variables);


#get references information
#print "Get references information.\n\n\n";
my $ref_info_pointer=sub_bioperl::pre_references_info(\%variables);
$variables{$ref_info_pointer}=$ref_info_pointer;
	
#print "get genome annotation from $variables{file_annot}\n";
#my($annot_gene_pointer, $annot_transcript_pointer, $annot_exon_pointer)=sub_rna::pre_genome_annotation(\%variables);
#$variables{annot_gene_pointer}=$annot_gene_pointer;
#$variables{annot_transcript_pointer}=$annot_transcript_pointer;
#$variables{annot_exon_pointer}=$annot_exon_pointer;

############################################################
#main running
if ($variables{mRNA_genome_mapping} eq 'yes' or $variables{mRNA_transcripts_assembling} eq 'yes' or $variables{mRNA_sequencing_control} eq 'yes'){
	print "\n\n\n### Main running! ###\n\n\n";
	my @sample_names=split(',', $variables{sample_names});
	while(1){
		my $threads_num=sub_basic::read_log($var_file, 'threads_num');
		$variables{threads_num}=$threads_num unless length($threads_num)==0 or $threads_num=~/[^0-9]/;
		if(threads->list() < $variables{threads_num} and @sample_names>0 ){
			my $sample_name=shift @sample_names;
			#print "$sample_name\n";
			threads->create(\&sub_mrna::main_running, \%variables, $sample_name);
		}
		#recover all threads
		foreach my $sub_thread( threads->list() ){#2
			$sub_thread->join() if $sub_thread->is_joinable();
		}#2
		last if threads->list()==0 and @sample_names==0;
		sleep 1;
	}
}

#merge cufflinks
if ($variables{mRNA_merge_seq} eq 'yes'){
	sub_mrna::mRNA_S4_merge_seq(\%variables);
}
#export 
#sub_mrna::mRNA_S6_statistics(\%variables);


print "\n\n\n### Analysis is done! GOOD LUCK! ###\n\n\n";



