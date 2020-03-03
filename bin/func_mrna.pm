#! /usr/bin/perl -w
use strict;
use warnings;
use List::MoreUtils;
use List::Util;
use File::Find;
use constant false => 0;
use constant true  => 1;

##########################################################
#subroutines:
use func_common; #sub_common::
use func_basic; #sub_basic::
use func_rna; #sub_rna::
####################################################
#
package sub_mrna;
###########################################################



#############
sub initiate_mRNA_variables{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	

	#switches
	$variables{'mRNA_genome_mapping'}='yes';
	$variables{'mRNA_sequencing_control'}='yes';
	$variables{'mRNA_transcripts_assembling'}='yes';
	$variables{'mRNA_merge_seq'}='yes';
	$variables{'mRNA_cuffdiff'}='no';
	#alignment
	$variables{ 'software_aligner'}='bowtie2';
	$variables{'software_mapper'}='tophat2';
	$variables{'software_assembler'}='cufflinks';
	$variables{'UN_query'}='yes';

	
	#
	$variables{'quality_filter'}=13;
	$variables{'QC_compression'}=1000;
	$variables{'statistical_items'}='NA';
	$variables{'unaligned_read_counts_filter'}=3;
	$variables{'RC_background'}=1;


	#directory
	#$variables{dir_raw_data}='/home/yuan/data_1/rawdata/a,/home/yuan/data_1/rawdata/b,/home/yuan/data_1/rawdata/c';
	#$variables{dir_result_array}='/home/yuan/eRNA/result,/home/yuan/eRNA/result/a,/home/yuan/eRNA/result/b';
	#initiate variables in variables.txt
	$variables{dir_bowtie}=$variables{dir_home}.'bowtie2/';
	$variables{exe_aligner}=$variables{dir_bowtie}.'bowtie2';
	$variables{dir_mapper}=$variables{dir_home}.'tophat/';
	$variables{exe_mapper}=$variables{dir_mapper}.'tophat2';
	$variables{dir_assembler}=$variables{dir_home}.'cufflinks/';
	$variables{exe_assembler}=$variables{dir_home}.'cufflinks/cufflinks';
	$variables{dir_trans_index}=sub_common::format_directory($variables{dir_mapper}.'transcriptome_index/');

	
	#print $variables{file_var};
	sub_common::hash_to_file(\%variables, $variables{file_var}, '=');
	return(\%variables);
}
####################################################
#
sub main_running{ #1 main_running begin
	my ($variables_pointer, $sample_name)=@_;
	my %variables=%$variables_pointer;
	
	#initiate time and sample dir
	$variables_pointer=sub_rna::initiate_starting_time(\%variables, $sample_name);
	%variables=%$variables_pointer;
	$variables{sample_tophat_dir}=sub_common::format_directory($variables{sample_dir}.'tophat');
	$variables{sample_alignment_file}=$variables{sample_tophat_dir}.'accepted_hits.bam';
	
	#Paired-ends alignment of *.fastq files against genome sequences
	if ($variables{mRNA_genome_mapping} eq 'yes'){
		sub_mrna::mRNA_S1_genome_mapping(\%variables);
	}
	
	if ($variables{mRNA_sequencing_control} eq 'yes'){
		sub_mrna::mRNA_S2_sequencing_control(\%variables);
	}
	
	#transcripts assembling
	if ($variables{mRNA_transcripts_assembling} eq 'yes'){
		sub_mrna::mRNA_S3_transcripts_assembling(\%variables);
	}
	
	#record ending time of a sample analysis
	sub_basic::initiate_ending_time($variables{sample_log}, $variables{beginning_time});
	
	print "Analysis of main_running on $variables{sample_name} is done!\n\n";
}#1 main_running end



#################################################
#
sub mRNA_S1_genome_mapping{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my $sample_name=$variables{sample_name};
	
	#print "\n\n Genome mapping of $sample_name:\n";
	my $tophat=join(' ', $variables{tophat_options}, '-o', $variables{sample_tophat_dir}, $variables{bowtie_index});
						
	#
	if ($variables{'sequencing_end'} eq 'paired-end'){
		print "mapping R1_fastq and R2_fastq files of $sample_name to the reference genome using TopHat!\n\n";
		print "$tophat $sample_info{$sample_name}->{R1_files} $sample_info{$sample_name}->{R2_files}\n";
		system("$tophat $sample_info{$sample_name}->{R1_files} $sample_info{$sample_name}->{R2_files}");
	}
	else{
		print "mapping R1_fastq files of $sample_name to the reference genome using TopHat!\n\n";
		print "$tophat $sample_info{$sample_name}->{R1_files}\n";
		system("$tophat $sample_info{$sample_name}->{R1_files}");
	}
	#bam to sam file
	print "convert the bam file in $variables{sample_tophat_dir}\n";
	sub_rna::bam_to_sam($variables{sample_tophat_dir}, $variables{sample_out_file});
	
	print "Alignment of $sample_name is done!\n";
}

###############################
sub mRNA_S2_sequencing_control{
	my($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my $sample_name=$variables{sample_name};
	my $sample_dir=$variables{sample_dir};
	my $annot_exon_pointer=$variables{annot_exon_pointer};
	my %annot_exon=%$annot_exon_pointer;
	
	print "Analysis of sequencing quality of $variables{sample_name}.\n";
	my $sam_file=$variables{sample_out_file}.'.sam';
	my $target_file=$variables{sample_out_file}.'_genes.txt';
	sub_rna::analyze_sam_file($sam_file, $target_file, \&sub_rna::target_references, \%annot_exon);
	
	print "\tSequencing quality control : saturation analysis\n";
	sub_rna::QC_saturation($target_file, $variables{sample_out_file}.'.SD');
	
	print "\tGet read counts\n";
	my $num_out_pointer=sub_rna::count_reads($target_file,  $variables{sample_out_file}.'_read_counts.txt');
	my %num_out=%$num_out_pointer;
	#raw reads
	my @R1_files=split(',', $sample_info{$sample_name}->{R1_files});
	my @R1_lines=map {sub_common::file_info($_, 'lines_number')} @R1_files;
	$num_out{raw_num}=(List::Util::sum @R1_lines)/4;
	sub_basic::refresh_log_hash($variables{sample_log}, \%num_out, '=');

}

#################################################
#
sub mRNA_S3_transcripts_assembling{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_name=$variables{sample_name};
	my $cufflinks_output_dir=$variables{sample_dir}.'cufflinks/';
	my $bam_file=$variables{sample_alignment_file};

	print "\n\n Transcripts assembling of $sample_name:\n";
	my $cufflinks=join('  ', $variables{cufflinks_options}, '-o', $cufflinks_output_dir, $bam_file);
	#assemble expressed genes and transcripts using cufflinks
	print "###$cufflinks\n";
	system("$cufflinks");
	
	#compare assembled transcripts to a reference annotation
	my $cuffcompare_exe=$variables{dir_assembler}.'cuffcompare';
	my $cuff_gtf_file=$cufflinks_output_dir.'transcripts.gtf';
	my $output_prefix=$cufflinks_output_dir.'cuffcomp';
	my $cuffcompare_options=join(" ", $cuff_gtf_file, '-o', $output_prefix);
	print "###$cuffcompare_exe $cuffcompare_options\n";
	system("$cuffcompare_exe $cuffcompare_options");
	print "Transcripts assembling of $sample_name is done!\n";
	#
}

##########################################
#merge  several cufflinks assembliers
sub mRNA_S4_merge_seq{
	my($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	
	########cuffmerge
	print "create assemblies.txt for each sample.\n";
	my $out_dir=sub_common::format_directory($variables{dir_result}.'cuffmerge');
	my $assemblies_file=$out_dir.'assemblies.txt';
	open my ($ASS), ">", $assemblies_file or die;
	foreach my $sample_name(sort(keys %sample_info)){
		my $sample_dir=$sample_info{$sample_name}->{sample_dir};
		print $ASS $sample_dir.'cufflinks/transcripts.gtf', "\n";
	}
	close($ASS);
	
	#run cuffmerge to create a signle merged transcriptome annotation
	my $cuffmerge=join(" ", $variables{dir_assembler}.'cuffmerge',  '-g', $variables{file_annot}, '-s', $variables{file_fa}, 
						'-p', $variables{threads_num}, '-o', $out_dir, );
	print "###$cuffmerge $assemblies_file \n\n\n";
	system(" $cuffmerge $assemblies_file");
	print "###cuffmerge is done \n\n";
	
	#run cuffquant
	#my $cuffquant=join(" ", $variables{dir_assembler}.'cuffquant',  '-p', $variables{threads_num}, '-o', $variables{dir_result}.'cuffquant', $variables{file_annot},);
	#print "###$cuffquant $sam_file \n\n\n";
	#system(" $cuffquant $sam_file");
	#print "####Cuffquant is done\n\n";
	
	print "####export delete and insertion\n\n";
	sub_mrna::tophat_cufflinks_export(\%variables, 'tophat/deletions.bed', 'deletions_counts.csv', '0,1', '4');
	sub_mrna::tophat_cufflinks_export(\%variables, 'tophat/insertions.bed', 'insertions_counts.csv', '0,1,2,3', '4');
	
	print "####export FPKM\n\n";
	sub_mrna::tophat_cufflinks_export(\%variables, 'cufflinks/genes.fpkm_tracking', 'genes_fpkm.csv', '0', '9');
	sub_mrna::tophat_cufflinks_export(\%variables, 'cufflinks/isoforms.fpkm_tracking', 'isoforms_fpkm.csv', '0', '9');
	
	print "####export read count\n\n";
	sub_rna::export_read_counts(\%variables);
	#
}

###########################################
sub  mRNA_S5_cuffdiff{
	my($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	
	print "\n\n Differential transcripts identification:\n";
	if ($variables{'diff_A'} eq 'NA' or $variables{'diff_B'} eq 'NA'){
		my @sample_candidates=split(',', $variables{diff_candidates});
		my @bam_files=map {$sample_info{$_}->{sample_dir}.'tophat/accepted_hits.bam'} @sample_candidates;
		my $bam_files_str=join(" ", @bam_files);
		$variables{cuffdiff_options}=join(" ", $variables{cuffdiff_options}, '-L', $variables{diff_candidates}, $bam_files_str);
	}
	else{
		my @A_sample_names=split(',', $variables{'diff_A'});
		my @A_bam_files;
		foreach(@A_sample_names){
			my $sample_dir=$sample_info{$_}->{sample_dir};
			push(@A_bam_files, $sample_dir.'tophat/accepted_hits.bam');
		}
		my $A_bam_files_str=join(',', @A_bam_files);
		#
		my @B_sample_names=split(',', $variables{'diff_B'});
		my @B_bam_files;
		foreach(@B_sample_names){
			my $sample_dir=$sample_info{$_}->{sample_dir};
			push(@B_bam_files, $sample_dir.'tophat/accepted_hits.bam');
		}
		my $B_bam_files_str=join(',', @B_bam_files);
		$variables{cuffdiff_options}=join(" ", $variables{cuffdiff_options}, '-L A,B', $A_bam_files_str, $B_bam_files_str);
	}

	#run cuffdiff 
	my $cuffdiff_file=$variables{dir_assembler}.'cuffdiff';
	print "#########\n $cuffdiff_file $variables{cuffdiff_options} \n\n\n";
	system(" $cuffdiff_file $variables{cuffdiff_options}");
	print "Differential transcripts identification is done!\n";
	#
}
#######################################
##############################################################################
#export: statistics.txt 
sub mRNA_S6_statistics{#1
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my @sample_names=split(",", $variables{sample_names});
	
	#read sample log
	my %counting_hash;
	foreach my $sample_name(@sample_names){#2
		my $sample_dir=$sample_info{$sample_name}->{sample_dir};
		my $sample_log=$sample_dir.$sample_name.'.log';
		open my($IN), "<", $sample_log or die "Error: canot open $sample_log!\n";
		while(<$IN>){#3
			chomp($_);
			my($key, $value)=split("=", $_);
			$counting_hash{$key}->{$sample_name}=$value;
		}#3
		close($IN);
		$counting_hash{raw_files}->{$sample_name}=split(',', $sample_info{$sample_name}->{'raw_files'});
	}#2
	
	#generate statistics result 
	sub_common::hash2_to_file(\%counting_hash, $variables{file_statistics});
	#
}


##########################################
#
sub tophat_cufflinks_export {
	my ($variables_pointer, $infile_tail, $outfile_tail, $row_index, $col_index)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my @sample_names=split(",", $variables{sample_names});
	my @row_arr=split(',', $row_index);
	my @col_arr=split(',', $col_index);
	#
	my %hash;
	foreach my $sample_name(@sample_names){#2
		my $sample_dir=$sample_info{$sample_name}->{sample_dir};
		printf("%s:%s\n", $sample_name, $sample_dir.$infile_tail);
		open my($IN), "<", $sample_dir.$infile_tail or die;
		my $first_line=<$IN>;#remove the first line
		while(<$IN>){#3
			chomp($_);
			my @items=split("\t", $_);
			my $row_name=join(':', map {$items[$_]} @row_arr);
			my $value=join(':', map {$items[$_]} @col_arr);
			$hash{$row_name}->{$sample_name}=$value;
		}#3
		close($IN);
	}#2
	
	#export
	my $outfile=$variables{dir_stat}.$outfile_tail;
	sub_common::hash2_to_file(\%hash, $outfile, ',');
	#
}

########################################
1;  # make sure the file returns true or require will not succeed!#
