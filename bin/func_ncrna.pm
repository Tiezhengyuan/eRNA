#! /usr/bin/perl -w
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq::Quality;
use List::MoreUtils;
use List::Util;
use File::Find;


##########################################################
#get the directory of perl scripts involved in Pscore
our $perl_dir=Cwd::getcwd();
#subroutines:
require "func_common.pm"; #sub_common::
require "func_basic.pm"; #sub_basic::
require "func_rna.pm"; #sub_rna::
#require $perl_dir."/functions_common.pm"; #sub_common::
#require $perl_dir."/functions_basic.pm"; #sub_basic::
#require $perl_dir."/functions_rna.pm"; #sub_rna::
####################################################
#
package sub_ncrna;
##############################################


#############
sub initiate_ncRNA_variables{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	
	#switches
	$variables{'ncRNA_read_rawdata'}='yes';
	$variables{'ncRNA_adapter_removal'}='yes';
	$variables{'ncRNA_seperate_alignment'}='yes';
	$variables{'ncRNA_iterative_alignment'}='yes';
	$variables{'ncRNA_counting'}='yes';
	
	#alignment
	$variables{'UN_query'}='yes';
	$variables{'query_len'}=16;
	#adapter trimming
	$variables{'match_len'}=8;
	$variables{'mismatch_allowed'}='yes';
	$variables{'mismatch_len'}=12;
	
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
	$variables{dir_bowtie}=$variables{dir_home}.'bowtie1/';
	$variables{software_aligner}='bowtie1';
	$variables{exe_aligner}=$variables{dir_bowtie}.'bowtie1';
	$variables{dir_mirdeep}=$variables{dir_home}.'mirdeep/';
	$variables{dir_mirspring}=$variables{dir_home}.'mirspring/';

	#
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
	
	#read raw data determined by sequencing analyzer
	if ($variables{ncRNA_read_rawdata} eq 'yes'){
		print "\n\n Read raw data files of $sample_name:\n";
		if ($variables{sequencing_end}==1){
				sub_ncrna::ncRNA_S1_read_rawdata_1(\%variables);
		}
		else{
				sub_ncrna::ncRNA_S1_read_rawdata_2(\%variables);
		}
		print "Raw data ($sample_name) reading is done!\n";
	}
		
	#adapter removal
	if ($variables{ncRNA_adapter_removal} eq 'yes'){
		print  "\n\n Remove adapter sequences of reads in $sample_name:\n";
		sub_ncrna::ncRNA_S2_adapter_removal(\%variables);
		print  "Adapter remove in $sample_name is done\n";
	}
		
	#alignment one by one using bowtie
	if ($variables{ncRNA_seperate_alignment} eq 'yes'){
		print  "\n\n Sequence alignments of $sample_name begin:\n";
		sub_ncrna::ncRNA_S3_seperate_alignment(\%variables);
		print  "Sequence alignments of $sample_name is done.\n";
	}


	#iterative alignment using bowtie
	if ($variables{ncRNA_iterative_alignment} eq 'yes'){
		print  "\n\n Iterative alignments of $sample_name begin:\n";
		sub_ncrna::ncRNA_S4_iterative_alignment(\%variables) ;
		print  "Iterative alignments of $sample_name is done.\n";
	}
		
	#record ending time of a sample analysis
	sub_basic::initiate_ending_time($variables{sample_log}, $variables{beginning_time});
	
	print "Analysis of $sample_name is done!\n\n";
}#1 main_running end



#####################################################
#read all fastq files and remove adapter sequences
#only R1
sub ncRNA_S1_read_rawdata_1{ #1
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my $sample_name=$variables{sample_name};
	
	my $raw_reads_num=0;
	my @R1_files=split(',', $sample_info{$sample_name}->{R1_files});
	#export pair ends of read sequences of R1 and R2
	open my($OUT), ">", $variables{sample_out_file}.".RC" or die; 
	foreach my $R1_file(@R1_files){#2
		#read R1 raw data
		print "\t read $R1_file!\n";
		my $in_obj = Bio::SeqIO->new(-file => $R1_file, -format => 'fastq');
		while (my $seq_obj = $in_obj->next_seq() ) {#3
				my @display_id=split(/:/, $seq_obj->display_id());
				my $coordinate=$display_id[-2]."_".$display_id[-1];
				my $R1_seq=$seq_obj->seq();
				$raw_reads_num++;
				print $OUT join("\t", $coordinate, $raw_reads_num, $R1_seq, 'NA', 'NA'), "\n";
		}#3
	}#2
	close($OUT);
	sub_basic::refresh_log($variables{sample_log}, "raw_reads_num", $raw_reads_num);
	
	#Sequencing quality of R1_file
	sub_rna::QC_Qvalues(\@R1_files, $variables{sample_out_file}.'_R1.QC');

	#return (\%variables);
}#1

###
#read all fastq files and remove adapter sequences
#R1 and R2
sub ncRNA_S1_read_rawdata_2{ #1
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my $sample_name=$variables{sample_name};
	
	my $raw_reads_num=0;
	my @R1_files=split(',', $sample_info{$sample_name}->{R1_files});
	#export pair ends of read sequences of R1 and R2
	open my($OUT), ">", $variables{sample_out_file}.".RC" or die; 
	foreach my $R1_file(@R1_files){#2
		my %pairs; # grabs the FASTQ parser, specifies the Illumina variant
		#read R1 raw data
		print "\t read $R1_file!\n";
		my $in_obj = Bio::SeqIO->new(-file => $R1_file, -format => 'fastq');
		while (my $seq_obj = $in_obj->next_seq() ) {#3
				my @display_id=split(/:/, $seq_obj->display_id());
				my $coordinate=$display_id[-2]."_".$display_id[-1];
				$pairs{$coordinate}->{R1_seq}=$seq_obj->seq();
				$raw_reads_num++;
				$pairs{$coordinate}->{R1_No}=$raw_reads_num;
				$pairs{$coordinate}->{R2_seq}="NA";
				$pairs{$coordinate}->{R2_No}="NA";
		}#3
		
		#read R2 raw data
		my $R2_file=$R1_file;
		$R2_file=~s/_R1/_R2/;
		print "\t read $R2_file!\n";
		$in_obj = Bio::SeqIO->new(-file => $R2_file, -format => 'fastq');
		while (my $seq_obj = $in_obj->next_seq() ) {#4
				my @display_id=split(/:/, $seq_obj->display_id());
				my $coordinate=$display_id[-2]."_".$display_id[-1];
				$pairs{$coordinate}->{R2_seq}=$seq_obj->seq();
				$pairs{$coordinate}->{R2_No}=$raw_reads_num;
				$raw_reads_num++;
		}#4
		
		#export
		foreach (keys %pairs){#3
				print $OUT join("\t", $_, $pairs{$_}->{R1_No}, $pairs{$_}->{R1_seq}, $pairs{$_}->{R2_No}, $pairs{$_}->{R2_seq}), "\n";
		}#3
	}#2
	close($OUT);
	sub_basic::refresh_log($variables{sample_log}, "raw_reads_num", $raw_reads_num);
	
	#Sequencing quality of R1_file
	sub_rna::QC_Qvalues(\@R1_files, $variables{sample_out_file}.'_R1.QC');
	my @R2_files=split(',', $sample_info{$sample_name}->{'R2_files'});
	sub_rna::QC_Qvalues(\@R2_files, $variables{sample_out_file}.'_R2.QC');

	#return (\%variables);
}#1

##########################################################
#adapter removal at the 3' end  
sub ncRNA_S2_adapter_removal{#1
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my $sample_name=$variables{sample_name};
	
	#generate the possible adapters for seq truncation
	my $adapters_pointer=sub_rna::adapters(\%variables);
	my %adapters=%$adapters_pointer;

	my (%insert, %insert_len);
	open my($RC), "<", $variables{sample_out_file}.".RC" or die;
	open my($INSERT), ">", $variables{sample_out_file}.".insert_RC" or die;
	while (<$RC>){#2
		chomp($_);
		my ($coordinate, $R1_seq_No, $R1_read_seq, $R2_seq_No, $R2_read_seq)=split("\t", $_);
		#truncate R1_seq
		my $R1_seq_trunc=sub_rna::adapter_3_truncation(\%variables, $R1_read_seq, \%adapters);
		my $R1_trunc_len=0;
		$R1_trunc_len=length($R1_seq_trunc) unless $R1_seq_trunc eq 'NA';
		$insert_len{$R1_trunc_len}++;
		if(exists $insert{$R1_seq_trunc}){
			$insert{$R1_seq_trunc}->{seq_No_info} .= ",".$R1_seq_No;
			$insert{$R1_seq_trunc}->{read_counts}++;
		}
		else{
			$insert{$R1_seq_trunc}->{seq_No_info}=$R1_seq_No;
			$insert{$R1_seq_trunc}->{read_counts}=1;
		}
		#truncate R2_seq
		if ($R2_read_seq eq 'NA'){
			print $INSERT "$coordinate\t", "$R1_seq_No\t", "$R1_seq_trunc\t", "NA\t", "NA\n";
		}
		else{
			my $R2_seq_trunc=sub_rna::adapter_5_truncation(\%variables, $R2_read_seq, \%adapters);
			my $R2_trunc_len=0;
			$R2_trunc_len=length($R2_seq_trunc) unless $R2_seq_trunc eq 'NA';
			$insert_len{$R2_trunc_len}++;
			if(exists $insert{$R2_seq_trunc}){
				$insert{$R2_seq_trunc}->{seq_No_info} .= ",".$R2_seq_No;
				$insert{$R2_seq_trunc}->{read_counts}++;
			}
			else{
				$insert{$R2_seq_trunc}->{seq_No_info}=$R2_seq_No;
				$insert{$R2_seq_trunc}->{read_counts}=1;
			}
			print $INSERT join("\t", $coordinate, $R1_seq_No, $R1_seq_trunc, $R2_seq_No, $R2_seq_trunc), "\n";
		}

	}#2
	close($RC);
	close($INSERT);
 
	#export the result of adapter removal process  
	#initiate variables
	my $insert_reads_num=0;
	my $query_reads_num=0; 
	my $nr_insert_reads_num=0;
	my $nr_query_reads_num=0; 
	open my($query), ">", $variables{sample_out_file}."_query.fa" or die;
	open my($TXT), ">", $variables{sample_out_file}.".query_RC" or die;
	delete $insert{'NA'};
	foreach my $seq(keys %insert){#2
		if (length($seq)>=$variables{query_len}){#4
			my @seq_No_info=split(",", $insert{$seq}->{seq_No_info});
			my $typical_seq_No=List::Util::min @seq_No_info;
			my $display_id=$sample_name.":".$typical_seq_No."_".$insert{$seq}->{read_counts};
			print $query ">$display_id\n", "$seq\n";
			print $TXT "$typical_seq_No\t", "$insert{$seq}->{read_counts}\t", "$seq\t", "$insert{$seq}->{seq_No_info}\n";
			$nr_query_reads_num ++;
			$query_reads_num += $insert{$seq}->{read_counts};
		}#4
		$insert_reads_num +=$insert{$seq}->{read_counts};
		$nr_insert_reads_num++;
	}#2
	close($query);
	close($TXT);
	sub_basic::refresh_log($variables{sample_log}, "insert_reads_num", $insert_reads_num);
	sub_basic::refresh_log($variables{sample_log}, "nr_insert_reads_num", $nr_insert_reads_num);
	sub_basic::refresh_log($variables{sample_log}, "query_reads_num", $query_reads_num);
	sub_basic::refresh_log($variables{sample_log}, "nr_query_reads_num", $nr_query_reads_num);
	
	#export insert length distribution
	sub_common::hash_to_file(\%insert_len, $variables{sample_out_file}.'.IL');
	#
	return(\%variables);
}#1


###################################################
#sequence alignment or genome mapping using bowtie 
#the software named bowtie-build and bowtie should be copied into the bowtie directory
sub ncRNA_S3_seperate_alignment{#1
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my $sample_name=$variables{sample_name};
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	
	#sequence alignment using bowtie
	my @index_names=split(",", $variables{index_seperate}) ;
	foreach my $index_name(@index_names){#2  #index name circyling
		#get all read files with the name tail *_combined.fasta 
		my $bowtie_index=$variables{dir_bowtie1}.$index_name;
		my $fasta_file=$variables{sample_out_file}."_query.fa";
		my $alignment_sam_output=$variables{sample_out_file}."_".$index_name.'.sam';
		print "Seperate alignment: $variables{bowtie_options} $bowtie_index -f $fasta_file -S $alignment_sam_output\n\n";
		system("$variables{bowtie_options} $bowtie_index -f $fasta_file -S $alignment_sam_output ");
		my $alignment_out_pointer=sub_rna::read_bowtie1_sam($fasta_file, $alignment_sam_output);
		my %alignment_out=%$alignment_out_pointer;
		
		#export miRNA quanlification (read counts of known miRNA versus mappable query read counts) 
		my $mappable_counts_pointer=$alignment_out{mappable_counts_pointer};
		my %mappable_counts=%$mappable_counts_pointer;
		open my($COU), ">", $variables{sample_out_file}."_".$index_name.".alignment" or die; 
		foreach my $ref_name( sort (keys %mappable_counts) ){#3
			print $COU join("\t", $sample_name, $index_name, $ref_name, 
				$mappable_counts{$ref_name}->{top_RC}, $mappable_counts{$ref_name}->{middle_RC},
				$mappable_counts{$ref_name}->{bottom_RC}, $mappable_counts{$ref_name}->{seq_info}, ), "\n";
			#print "seperate alignment: $sample_name: $ref_name.\n";
		}#3
		close($COU);
		
		#update sample log file
		my $num_out_pointer=$alignment_out{num_out_pointer};
		my %num_out=%$num_out_pointer;
		foreach my $num_name(keys %num_out){
			sub_basic::refresh_log($variables{sample_log}, $index_name.':'.$num_name, $num_out{$num_name}, '=');
		}
		print  "mappable redundant number: $num_out{nr_query_reads_num}\t";
		print  "$num_out{nr_mappable_reads_num}\t", "$num_out{nr_unmappable_reads_num}\n";
		
		#export unalignment and multiple alignment sequences
		my $un_alignment_pointer=$alignment_out{un_alignment_pointer};
		my %un_alignment=%$un_alignment_pointer;
		open my($UN), ">", $variables{sample_out_file}."_".$index_name."_unalignment.fa" or die; 
		foreach my $query_name(sort {$un_alignment{$b}->{RC}<=>$un_alignment{$a}->{RC}} (keys %un_alignment)){
			print $UN ">$query_name\n", "$un_alignment{$query_name}->{seq}\n";
		}
		close($UN);
		my $multiple_alignment_pointer=$alignment_out{multiple_alignment_pointer};
		my %multiple_alignment=%$multiple_alignment_pointer;
		open my($MUL), ">", $variables{sample_out_file}."_".$index_name."_multialignment.fa" or die; 
		foreach my $query_name(sort {$multiple_alignment{$b}->{RC}<=>$multiple_alignment{$a}->{RC}} (keys %multiple_alignment)){
			print $MUL ">$query_name\n", "$multiple_alignment{$query_name}->{seq}\n";
		}
		close($MUL);
		
		#export sequencing depth(SD)
		my @ref_names;
		my $ref_num=0;
		my $SD_counts_pointer=$alignment_out{SD_counts_pointer};
		my %SD_counts=%$SD_counts_pointer;
		$SD_counts{'1'}=0 unless exists $SD_counts{'1'};
		open my($SD), ">", $variables{sample_out_file}."_".$index_name.".SD" or die; 
		foreach my $seq_No(sort {$a<=>$b} keys(%SD_counts) ){
			my @sub_ref_names=split(";;", $SD_counts{$seq_No});
			for my $ref_name(@sub_ref_names){
				push(@ref_names, $ref_name) unless List::Util::first {$ref_name eq $_} @ref_names;
			}
			$ref_num=@ref_names;
			print $SD "$seq_No\t", "$ref_num\n";
		}
		my $seq_No_max=List::Util::max keys %SD_counts;
		my $raw_reads_num=sub_basic::read_log($variables{sample_log}, 'raw_reads_num');
		print $SD "$raw_reads_num\t", "$ref_num\n" if $raw_reads_num > $seq_No_max;
		close($SD);
	}#2 index name circyling end
  
	#return(\%variables);
}#1


##############################################################
#iterativesequence alignment using bowtie
sub ncRNA_S4_iterative_alignment{#1
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my $sample_name=$variables{sample_name};
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	
	#initiate the sample log
	my %hash;
	open my($IN), "<", $variables{sample_log} or die;	#read old data
	while (<$IN>){
		chomp($_);
		my($name, $value)=split("=", $_);
		$hash{$name}=$value;
	}
	close($IN);
	open my($OUT), ">", $variables{sample_log} or die;  #export refreshed data
	foreach my $key( sort {$a cmp $b} (keys %hash ) ){
		print $OUT "$key=$hash{$key}\n" unless $key=~/^iterative/;
	}
	close($OUT);
	
	#initiate temporary query fasta file and temporary sam output file
  	my $fasta_file=$variables{sample_dir}.$sample_name."_query.fa";
  	my $fasta_tmp=$fasta_file.'_tmp';
	system("cp $fasta_file $fasta_tmp");
	my $sam_tmp=$variables{sample_dir}.$sample_name.'_iterative.sam_tmp';
	
	#iterative alignment
	my (%total_num_out, %SD_counts); #total number counting
	my @index_names=split(",", $variables{index_iterative});
	open my($SAM), ">", $variables{sample_dir}.$sample_name.'_iterative_combined.sam' or die;
	open my($COU), ">", $variables{sample_dir}.$sample_name.'_iterative.alignment' or die;
	for (my $i=0; $i<@index_names; $i++){#2  index_name circyling
		#bowtie alignment
		my $index_name=$index_names[$i];
		my $bowtie_index=$variables{dir_bowtie1}.$index_name;
		print  "Iterative alignment of $sample_name against $bowtie_index!\n";
		system("$variables{bowtie_options} $bowtie_index -f $fasta_tmp -S $sam_tmp ");
		my $alignment_out_pointer=sub_rna::read_bowtie1_sam($fasta_tmp, $sam_tmp);
		my %alignment_out=%$alignment_out_pointer;
		
		#export miRNA quanlification (read counts of known miRNA versus mappable query read counts) 
		my $mappable_counts_pointer=$alignment_out{mappable_counts_pointer};
		my %mappable_counts=%$mappable_counts_pointer;
		foreach my $ref_name( sort (keys %mappable_counts) ){#3
			print $COU join("\t", $sample_name, 'iterative_'.$i.'_'.$index_name, $ref_name, 
				$mappable_counts{$ref_name}->{top_RC}, $mappable_counts{$ref_name}->{middle_RC},
				$mappable_counts{$ref_name}->{bottom_RC}, $mappable_counts{$ref_name}->{seq_info}, ), "\n";
			#print "Iterative alignment: $sample_name: $ref_name.\n";
		}#3
		#export unique alignment sequences into sam file combined
		my $alignment_pointer=$alignment_out{alignment_pointer};
		my %alignment=%$alignment_pointer;
		foreach my $line(keys %alignment){
			print $SAM "$alignment{$line}->{sam_line}\n";
		}
		
		#export unalignment alignment sequences
		my $un_alignment_pointer=$alignment_out{un_alignment_pointer};
		my %un_alignment=%$un_alignment_pointer;
		open my($UN), ">", $variables{sample_dir}.$sample_name.'_iterative_'.$i.'_'.$index_name."_unalignment.fa" or die; 
		open my($TMP), ">", $fasta_tmp or die; 
		foreach my $query_name(sort {$un_alignment{$b}->{RC}<=>$un_alignment{$a}->{RC}} (keys %un_alignment) ){
			print $UN ">$query_name\n", "$un_alignment{$query_name}->{seq}\n";
			print $TMP ">$query_name\n", "$un_alignment{$query_name}->{seq}\n";
		}
		close($UN);
		close($TMP);
		
		#export multiple alignment sequences
		my $multiple_alignment_pointer=$alignment_out{multiple_alignment_pointer};
		my %multiple_alignment=%$multiple_alignment_pointer;
		open my($MUL), ">", $variables{sample_dir}.$sample_name.'_iterative_'.$i.'_'.$index_name."_multialignment.fa" or die; 
		foreach my $query_name(sort {$multiple_alignment{$b}->{RC}<=>$multiple_alignment{$a}->{RC}} (keys %multiple_alignment) ){
			print $MUL ">$query_name\n", "$multiple_alignment{$query_name}->{seq}\n";
		}
		close($MUL);
		
		#update %SD_counts
		my $SD_counts_pointer=$alignment_out{SD_counts_pointer};
		my %tmp_SD=%$SD_counts_pointer;
		while (my($tmp_seq_No, $tmp_ref_name)=(each %tmp_SD)){
			if (exists $SD_counts{$tmp_seq_No}){		$SD_counts{$tmp_seq_No} .= ';;'.$tmp_ref_name;	}
			else{	$SD_counts{$tmp_seq_No}=$tmp_ref_name;	}
		}
		#update sample log file
		my $num_out_pointer=$alignment_out{num_out_pointer};
		my %num_out=%$num_out_pointer;
		while( my($num_name, $value)=(each %num_out) ){
			sub_basic::refresh_log($variables{sample_log}, 'iterative_'.$i.'_'.$index_name.':'.$num_name, $value);
			$total_num_out{'iterative:'.$num_name} +=$value;
		}
		print  "Iterative mappable redundant number against $index_name: $num_out{nr_query_reads_num}\t";
		print  "$num_out{nr_mappable_reads_num}\t", "$num_out{nr_unmappable_reads_num}\n";
		
	}#2 index name circyling end
	close($COU);
	close($SAM);
	
	#export sequencing depth(SD)
	my @ref_names;
	my $ref_num=0;
	$SD_counts{'1'}=0 unless exists $SD_counts{'1'};
	open my($SD), ">", $variables{sample_dir}.$sample_name."_iterative.SD" or die; 
	foreach my $seq_No(sort {$a<=>$b} keys(%SD_counts) ){
		my @sub_ref_names=split(";;", $SD_counts{$seq_No});
		for my $ref_name(@sub_ref_names){
			push(@ref_names, $ref_name) unless List::Util::first {$ref_name eq $_} @ref_names;
		}
		$ref_num=@ref_names;
		print $SD "$seq_No\t", "$ref_num\n";
	}
	my $seq_No_max=List::Util::max keys %SD_counts;
	my $raw_reads_num=sub_basic::read_log($variables{sample_log}, 'raw_reads_num');
	print $SD "$raw_reads_num\t", "$ref_num\n" if $raw_reads_num > $seq_No_max;
	close($SD);
		
	#total alignment analysis
	while( my($num_name, $value)=(each %total_num_out) ){
		sub_basic::refresh_log($variables{sample_log}, $num_name, $value);
	}
	
	#clear temporary files
	system("rm $sam_tmp $fasta_tmp");
	print  "\n Iterative sequences alignment of $sample_name using bowtie is done.\n\n\n";

}#1

##############################################################################
#
sub ncRNA_S5_counting{#1
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my @sample_names=split(",", $variables{sample_names});
	
	#get sample log information
	my (%statistics, @statistical_names);
	foreach my $sample_name(@sample_names){
		my $sample_log=$variables{dir_result}.$sample_name.'/'.$sample_name.'.log';
		open my($IN), "<",  $sample_log or die;
		while(<$IN>){
			chomp($_);
			my($statistical_name, $value)=split("=", $_);
			$statistics{$sample_name}->{$statistical_name}=$value;
			push(@statistical_names, $statistical_name) unless List::Util::first {$_ eq $statistical_name} @statistical_names;
		}
		close($IN);	
	}
	
	print "\n\nIntegrate and export results into statistics.txt!\n\n\n";
	open my($STI), ">", $variables{dir_statistics}.'statistics.txt' or die;
	print $STI join("\t", 'Name_of_counting', @sample_names), "\n";
	foreach my $statistical_name(@statistical_names){
		my @values;
		foreach my $sample_name(@sample_names){
			$statistics{$sample_name}->{$statistical_name}=0 unless exists $statistics{$sample_name}->{$statistical_name};
			push(@values, $statistics{$sample_name}->{$statistical_name});
		}
		print $STI join("\t", $statistical_name, @values), "\n";
	}
	
	print "generate data frame of expression level.\n";
	my (%RC_df, %seq_df);
	my $files_pointer=sub_common::files_list($variables{dir_result}, 'incrusive_file');
	my @files=@$files_pointer;
	my @alignment_files=grep(/\.alignment$/, @files);
	foreach my $alignment_file(@alignment_files){#2
		#print "$alignment_file\n";
		open my($IN), "<", $alignment_file or die;
		while(my $line=<$IN>){#4
			chomp($line);
			my ($sample_name, $index_name, $ref, $top_RC, $middle_RC, $bottom_RC, $seq_info)=split("\t", $line);
			$index_name='iterative' if $index_name=~/iterative/;
			$RC_df{$index_name}->{$ref}->{$sample_name}=$middle_RC;
		}#4
		close($IN);
	}#2
	
	#export %RC_df
	foreach my $index_name(keys %RC_df){
		print "$index_name\n";
		my $pointer=$RC_df{$index_name};
		my %sample_df=%$pointer;
		
		#export %data_frame
		#my %hash=();
		open my($RC), ">", $variables{dir_statistics}.'df_totalRC_'.$index_name.'.txt' or die;
			print $RC join("\t", 'Ref_names', @sample_names), "\n";
		open my($MRC), ">", $variables{dir_statistics}.'df_mappableRC_'.$index_name.'.txt' or die;
			print $MRC join("\t", 'Ref_names', @sample_names), "\n";
		open my($IRC), ">", $variables{dir_statistics}.'df_insertRC_'.$index_name.'.txt' or die;
			print $IRC join("\t", 'Ref_names', @sample_names), "\n";
		foreach my $ref_name(sort(keys %sample_df)){#3
			my (@RC, @MRC, @IRC);
			foreach my $sample_name(@sample_names){#4
				#print "##$ref_name:$sample_name\n";
				$sample_df{$ref_name}->{$sample_name}=0 unless exists $sample_df{$ref_name}->{$sample_name};
				my $RC=int($sample_df{$ref_name}->{$sample_name}+0.5);
				push(@RC, $RC );
				my $mappableRC=int( $RC * $variables{normalization_scaling} / $statistics{$sample_name}->{$index_name.':mappable_reads_num'} +0.5 );
				push(@MRC, $mappableRC );
				my $insertRC=int( $RC * $variables{normalization_scaling} / $statistics{$sample_name}->{'insert_reads_num'} +0.5 );
				push(@IRC, $insertRC );
			}#4
			#print "$ref_name\n" unless @RC==192 and @MRC==192 and @IRC==192;
			print $RC  join("\t", $ref_name, @RC), "\n" if ((List::Util::sum @RC)/@RC) >=$variables{RC_background};
			print $MRC  join("\t", $ref_name, @MRC), "\n" if ((List::Util::sum @MRC)/@MRC) >=$variables{RC_background};
			print $IRC  join("\t", $ref_name, @IRC), "\n" if ((List::Util::sum @IRC)/@IRC) >=$variables{RC_background};
		}#3
		close($RC);
		close($MRC);
		close($IRC);
	}#2


}#1


#####################
#download reference from mirbase
sub fetch_mirbase{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	
	my @file_names=('hairpin.fa.gz','mature.fa.gz', 'organisms.txt.gz');
	#download
	my $mirbase_url="ftp://mirbase.org/pub/mirbase/CURRENT/";
	foreach my $file_name(@file_names){
		sub_basic::download_url_file($mirbase_url, $file_name, $variables{dir_bowtie1});
	}
	#
	#return();
}


#######################################################################
#generate specie miRNA fasta file based on miRBase fasta files
sub specie_miRNA{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my$precursor_fasta_file=$variables{ncRNA_mirbase_precursor_file};
	my $mature_fasta_file=$variables{ncRNA_mirbase_mature_file};
	my $specie=$variables{ncRNA_species};
	my $dir = sub_common::format_directory($variables{dir_bowtie1});
	my $specie_precursor_fasta_file=$dir.$specie.'_precursor_miRNA.fa';
	my $other_precursor_fasta_file=$dir.'other_precursor_miRNA.fa';
	my $specie_mature_fasta_file=$dir.$specie.'_mature_miRNA.fa';
	my $other_mature_fasta_file=$dir.'other_mature_miRNA.fa';
	
	#extract precursor miRNA
	my (%specie_precursor, %other_precursor);
	my $in_obj=Bio::SeqIO->new(-file=> $precursor_fasta_file, -format=>"fasta");
	while (my $seq_obj=$in_obj->next_seq){
		my $seq=$seq_obj->seq();
		$seq=~tr/a-z/A-Z/;
		$seq=~s/U/T/g;
		my @array=split(" ", $seq_obj->display_id());
		my $display_id=$array[0];
		if ($display_id=~/$specie/){
			$specie_precursor{$seq}=(exists $specie_precursor{$seq}) ? $specie_precursor{$seq}.','.$display_id : $display_id;
		}
		else{
			$other_precursor{$seq}=(exists $other_precursor{$seq}) ? $other_precursor{$seq}.','.$display_id : $display_id;
		}
	}
	#export precursor miRNA
	open my($out1), ">", $specie_precursor_fasta_file or die;
	while( my($seq, $miRNA_name)=each(%specie_precursor) ){
		print $out1 ">precursor_miRNA:$miRNA_name\n", "$seq\n";
	}
	close($out1);
	open my($out2), ">", $other_precursor_fasta_file or die;
	while( my($seq, $miRNA_name)=each(%other_precursor) ){
		print $out2 ">precursor_miRNA:$miRNA_name\n", "$seq\n";
	}
	close($out2);
	
	#extract mature miRNA
	my (%specie_mature, %other_mature);
	$in_obj=Bio::SeqIO->new(-file=> $mature_fasta_file, -format=>"fasta");
	while (my $seq_obj=$in_obj->next_seq){
		my $seq=$seq_obj->seq();
		$seq=~tr/a-z/A-Z/;
		$seq=~s/U/T/g;
		my @array=split(" ", $seq_obj->display_id());
		my $display_id=$array[0];
		if ($display_id=~/$specie/){
			$specie_mature{$seq}=(exists $specie_mature{$seq}) ? $specie_mature{$seq}.','.$display_id : $display_id;
		}
		else{
			$other_mature{$seq}=(exists $other_mature{$seq}) ? $other_mature{$seq}.','.$display_id : $display_id;
		}
	}
	#export precursor miRNA
	open my($out3), ">", $specie_mature_fasta_file or die;
	while( my($seq, $miRNA_name)=each(%specie_mature) ){
		print $out3 ">mature_miRNA:$miRNA_name\n", "$seq\n";
	}
	close($out3);
	open my($out4), ">", $other_mature_fasta_file or die;
	while( my($seq, $miRNA_name)=each(%other_mature) ){
		print $out4 ">mature_miRNA:$miRNA_name\n", "$seq\n";
	}
	close($out4);
	
	#my $out=join(" ", $specie_matured_fasta_file, $other_matured_fasta_file, $specie_precursor_fasta_file);
	#return($out);
}




########################################
1;  # make sure the file returns true or require will not succeed!#
