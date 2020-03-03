#! /usr/bin/perl -w
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq::Quality;
use List::MoreUtils;
use List::Util;
use File::Find;

##########################################################
#subroutines:
use func_basic;
use func_common;
use func_bioseq; # sub_bioseq
#require $perl_dir."/functions_basic.pm"; #sub_basic::
####################################################
#


#the file constains all subroutines required for running E_RNA_pipeline

package sub_rna;
###################################################
#used in GUI
sub initiate_directories{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	$variables{dir_home}=sub_common::format_directory($variables{dir_home});
	$variables{dir_result_array}=$variables{dir_result} unless exists $variables{dir_result_array};
	
	#default stat file
	$variables{dir_stat}=sub_common::format_directory($variables{dir_result}.'statistics');
	$variables{file_stat}=$variables{dir_stat}.'statistics.txt';
	
	#default log file
	$variables{dir_log}=sub_common::format_directory($variables{dir_result}.'sample_log');
	$variables{file_total_log}=$variables{dir_log}.'Total.log';
	$variables{file_log}=$variables{dir_result}.'output.log';
	$variables{file_time_monitor_log}=$variables{dir_result}.'time_monitor.log'; #sub_basic::initial_log_files()
	$variables{file_system_monitor_log}=$variables{dir_result}.'system_monitor.log';
	$variables{file_storage_log}=$variables{dir_result}.'storage.log';
	
	#default
	$variables{file_references_txt}=$variables{dir_result}.'references.txt'; #sub_basic::pre_ref_info()
	#$variables{file_targets_df}=$variables{dir_result}.'expression_profilings.txt';
	#rule the directories used for saving results determined by eRNA when large data analysis
	

 	return(\%variables);
}
######################
#initiate start/end time and sample directory
sub initiate_starting_time{
	my ($variables_pointer, $sample_name)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	
	#refresh sample_dir
	$variables{sample_name}=$sample_name; #sample name
	$variables{sample_dir}=$sample_info{$sample_name}->{sample_dir};
	$variables{sample_out_file}=$variables{sample_dir}.$sample_name;
	$variables{sample_log}=$variables{sample_out_file}.'.log';
	#Note: The result directory is temporarily changed to $variables{sample_dir} for $variables{sample_name}

	#refresh sample_log
	$variables{beginning_time}=join(',', localtime(time) );
	#unit is second
	$variables{supposed_time}=int($sample_info{$variables{sample_name} }->{'supposed_time'}/$sample_info{$variables{sample_name} }->{'parallel_parts'}); 
	sub_basic::refresh_log($variables{sample_log}, "sample_name", $variables{sample_name});
	sub_basic::refresh_log($variables{sample_log}, "beginning_time" , $variables{beginning_time});
	sub_basic::refresh_log($variables{sample_log}, "supposed_time" , $variables{supposed_time});
	sub_basic::refresh_log($variables{sample_log}, "status" , 'on');
	
	#
	return(\%variables);
}

##############################################
sub bam_to_sam{
	my ($dir, $out_file)=@_;
	
	#merge both $sam_file
	my $accepted_bam_file=$dir.'accepted_hits.bam';
	my $unmapped_bam_file=$dir.'unmapped.bam';
	my $tmp_bam_file=$dir.'tmp.bam';
	print "samtools merge $tmp_bam_file $accepted_bam_file $unmapped_bam_file\n";
	system("samtools merge -f $tmp_bam_file $accepted_bam_file $unmapped_bam_file");
	
	#sort bam file, and convert bam file to sam file
	my $sam_file=$out_file.'.sam';
	print "samtools sort $tmp_bam_file -o $sam_file\n";
	system("samtools sort $tmp_bam_file -o $sam_file");
	#
	#system("rm $tmp_bam_file");
	
	#return($sam_file)
}
################
sub sam_to_bam{
	my ($sample_dir)=@_;
	
	#print "Convert *.sam files into *bsam files\n";
	my $sam_file=$sample_dir.'accepted_hits.sam';
	my $file=$sample_dir.'accepted_hits.tmp';
	my $bam_file=$sample_dir.'accepted_hits';
	system("samtools view -bS $sam_file > $file");
	system("samtools sort $file $bam_file");
	system("rm $file");
	#
	$sam_file=$sample_dir.'unmapped.sam';
	$file=$sample_dir.'unmapped.tmp';
	$bam_file=$sample_dir.'unmapped';
	system("samtools view -bS $sam_file > $file");
	system("samtools sort $file $bam_file");
	system("rm $file");
	
}

######################################
sub read_sam_line{
	my ($sam_line)=@_;
	chomp($sam_line);
	my @items=split("\t", $sam_line);
	
	my %single;
	$single{query}=$items[0];
	#query name
	my @query_info=split(":", $single{query});
	$single{coordinate}=($query_info[-2]) ? $query_info[-2].":".$query_info[-1] : $query_info[-1];
	$single{flag}=sprintf("%08b", $items[1]); #convert to binary number
	#if(substr($single{flag}, 5, 1)==1){#no reported alignment (+4)
	#	$single{ref}='*';
	#	$single{start}=0; 
	#}
	$single{ref}=$items[2];
	$single{start}=$items[3]; 
	$single{MapQ}=$items[4]; 
	$single{CIGAR}=$items[5]; 
	$single{mate_ref}=$items[6]; 
	$single{mate_offset}=$items[7];
	#
	if(substr($single{flag}, 3, 1)==1){	#reversed reference strand	
		$single{seq}=sub_bioseq::reverse_complement($items[9]);
		$single{alignment}='revcom';
	}
	else{
		$single{seq}=$items[9];
		$single{alignment}='itself';
	}
	$single{len}=length($items[9]);
	$single{Qvalues}=$items[10]; 
	
	# split CIGAR
	$single{match_str}=sub_rna::CIGAR_splicing($single{start}, $single{len}, $single{CIGAR});
	
	#option fields
	for (my $i=11; $i<@items; $i++){
		my($a, $b, $c)=split(':', $items[$i]);
		$single{$a}=$c;
	}
	#unique: $single{NH}==1
	#multiple: $single{NH}>1
	#no alignment: $single{NH}==0
	$single{NH}=0 unless $single{NH};
	
	return(\%single);
}

###############################
#
sub CIGAR_splicing{
	my($start, $len, $CIGAR)=@_;
	my $match_str='NA';
	
	if ($CIGAR eq '*'){	#no alignment
		$match_str=$start.'-'.($start+$len-1);
	}
	else{#2
		my @cigar_pos=split(/[A-Z]/, $CIGAR);
		my @cigar_op=split(/\d+/, $CIGAR);
		shift @cigar_op;
		#
		my $end=$start;
		for (my $i=0; $i<@cigar_op; $i++){#3
			if ($cigar_op[$i] eq 'M' or $cigar_op[$i] eq 'D'){#M is alignment match. D is deletion
				$end += $cigar_pos[$i];
			}
		
			#junction exon
			if ($cigar_op[$i] eq 'N' or $i==(@cigar_op-1)){ #N is gap
				$match_str = ($match_str eq 'NA' ) ? $start.'-'.($end-1) : $match_str.','.$start.'-'.($end-1);
				$start=$end+$cigar_pos[$i];
				$end=$start;
			}
		}#3
	}#2
	return($match_str);
}
#################################
sub read_gtf_line{
	my ($gtf_line)=@_;
	chomp($gtf_line);
	
	my %gtf;
	my @items=split("\t", $gtf_line);
	$gtf{chr}=$items[0];
	$gtf{type}=$items[2];
	$gtf{start}=$items[3];
	$gtf{start_len}=int($gtf{start}/1e4);
	$gtf{end}=$items[4];
	$gtf{end_len}=int($gtf{end}/1e4);
	#
	my $annot=$items[8];
	$annot=~s/;$//;
	my @annot=split(/;\s/, $annot);
	foreach my $ann(@annot){
		my($type, $id)=split(/\s/, $ann);
		$id=~s/"//g;
		$gtf{$type}=$id;
		#print "#$gtf{$type}#\n";
	}

	return(\%gtf);
}
#####################
sub overlapped_genes{
	my ($genes_pointer)=@_;
	my %genes=%$genes_pointer;
	
	#assign regions 
	my $pointer=sub_bioseq::region_fragments(\%genes);
	my %region_genes=%$pointer;
	
	my %overlapped;
	foreach my $region(keys %region_genes){#2
		my $pointer=$region_genes{$region};
		my %sub_genes=%$pointer;
		my %poss_genes=%sub_genes;
		foreach my $id(keys %sub_genes){#3
			my $start=$genes{$id}->{start};
			my $end=$genes{$id}->{end};
			foreach my $poss_id(keys %poss_genes){#4
				my $poss_start=$poss_genes{$poss_id}->{start};
				my $poss_end=$poss_genes{$poss_id}->{end};
				my $judging=sub_common::overlapping($poss_start, $poss_end, $start, $end);
				if($judging eq 'cross' or $judging eq 'in' or $judging eq 'cover'){#5
					$overlapped{$id}= (exists $overlapped{$id}) ? $overlapped{$id}.','.$poss_id : $poss_id;
					$overlapped{$poss_id}= (exists $overlapped{$poss_id}) ? $overlapped{$poss_id}.','.$id : $id;
					#print join("\t", $id, $start, $end, $judging, $poss_id, $poss_start, $poss_end ), "\n";
				}#5
			}#4
			delete $poss_genes{$id} ;
		}#3
	}#2

	return(\%overlapped);
}
#####################################
sub pre_gene_annotation{
	my($gtf_file, $out_file)=@_;
	
	my %annot_info;
	open my($IN), "<", $gtf_file or die;
	while(<$IN>){#2
		my $gtf_pointer=sub_rna::read_gtf_line($_);
		my %gtf=%$gtf_pointer;
		if ($gtf{type} eq 'gene'){
			$annot_info{$gtf{chr}}->{$gtf{gene_id}}->{gene_name}=$gtf{gene_name};
			$annot_info{$gtf{chr}}->{$gtf{gene_id}}->{gene_biotype}=$gtf{gene_biotype};
			$annot_info{$gtf{chr}}->{$gtf{gene_id}}->{start}=$gtf{start};
			$annot_info{$gtf{chr}}->{$gtf{gene_id}}->{end}=$gtf{end};
		}
		elsif($gtf{type} eq 'transcript'){
			$annot_info{$gtf{chr}}->{$gtf{gene_id}}->{transcript_id} .= ','.$gtf{transcript_id};
			$annot_info{$gtf{chr}}->{$gtf{gene_id}}->{transcript_id}=~s/^\,//;
		}
		elsif($gtf{type} eq 'exon'){
			$annot_info{$gtf{chr}}->{$gtf{gene_id}}->{exons} .= ','.$gtf{start}.'-'.$gtf{end};
			$annot_info{$gtf{chr}}->{$gtf{gene_id}}->{exons}=~s/^\,//;
		}
	}#2
	close($IN);
	
	#export
	open my($OUT), ">", $out_file or die;
	print $OUT join("\t", 'chr', 'gene_id', 'gene_name', 'start', 'end', 'biotype', 'transcript_id', 'exons', 'overlapped', ), "\n";
	foreach my $chr(keys %annot_info){#2
		my $pointer=$annot_info{$chr};
		my %hash=%$pointer;
		my $over_genes_pointer=sub_rna::overlapped_genes(\%hash);
		my %over_genes=%$over_genes_pointer;
		foreach my $id(keys %hash){#3
			my $gene_name=$hash{$id}->{gene_name} ? $hash{$id}->{gene_name} : 'NA';
			my $biotype=$hash{$id}->{gene_biotype};
			my $start=$hash{$id}->{start};
			my $end=$hash{$id}->{end};
			my $transcript_id=$hash{$id}->{transcript_id};
			my $exons=$hash{$id}->{exons};
			my $overlapped=(exists $over_genes{$id}) ? $over_genes{$id} : 'NA';
			$annot_info{$chr}->{$id}->{overlapped}=$overlapped;
			print $OUT join("\t", $chr, $id, $gene_name, $start, $end, $biotype, $transcript_id, $exons, $overlapped, ), "\n";
		}#3
	}#2
	close($OUT);
	
	return(\%annot_info);
}
########################################
sub pre_transcript_annotation{
	my($gtf_file, $out_file)=@_;
	
	my %annot_info;
	open my($IN), "<", $gtf_file or die;
	while(<$IN>){#2
		my $gtf_pointer=sub_rna::read_gtf_line($_);
		my %gtf=%$gtf_pointer;
		if($gtf{type} eq 'transcript'){
			$annot_info{$gtf{chr}}->{$gtf{transcript_id}}->{gene_id} = $gtf{gene_id};
			$annot_info{$gtf{chr}}->{$gtf{transcript_id}}->{gene_name}=$gtf{gene_name};
			$annot_info{$gtf{chr}}->{$gtf{transcript_id}}->{transcript_name}=$gtf{transcript_name};
			$annot_info{$gtf{chr}}->{$gtf{transcript_id}}->{transcript_biotype}=$gtf{transcript_biotype};
			$annot_info{$gtf{chr}}->{$gtf{transcript_id}}->{start}=$gtf{start};
			$annot_info{$gtf{chr}}->{$gtf{transcript_id}}->{end}=$gtf{end};
		}
		elsif($gtf{type} eq 'exon'){
			$annot_info{$gtf{chr}}->{$gtf{transcript_id}}->{exons} .= ','.$gtf{start}.'-'.$gtf{end};
			$annot_info{$gtf{chr}}->{$gtf{transcript_id}}->{exons}=~s/^\,//;
		}
	}#2
	close($IN);
	
	#export
	open my($OUT), ">", $out_file or die;
	print $OUT join("\t", 'chr', 'transcript_id', 'transcript_name', 'start', 'end', 'biotype', 'exons', 'overlapped', 'gene_id', 'gene_name', ), "\n";
	foreach my $chr(keys %annot_info){#2
		my $pointer=$annot_info{$chr};
		my %hash=%$pointer;
		my $over_transcripts_pointer=sub_rna::overlapped_genes(\%hash);
		my %over_transcripts=%$over_transcripts_pointer;
		foreach my $id(keys %hash){#3
			my $gene_id=$hash{$id}->{gene_id};
			my $gene_name=$hash{$id}->{gene_name} ? $hash{$id}->{gene_name} : 'NA';
			my $start=$hash{$id}->{start};
			my $end=$hash{$id}->{end};
			my $exons=$hash{$id}->{exons};
			my $transcript_name=$hash{$id}->{transcript_name} ? $hash{$id}->{transcript_name} : 'NA';
			my $transcript_biotype=$hash{$id}->{transcript_biotype};
			my $overlapped=(exists $over_transcripts{$id}) ? $over_transcripts{$id} : 'NA';
			$annot_info{$chr}->{$id}->{overlapped}=$overlapped;
			print $OUT join("\t", $chr, $id, $transcript_name, $start, $end, $transcript_biotype, $exons, $overlapped, $gene_id, $gene_name, ), "\n";
		}#3
	}#2
	close($OUT);
	
	return(\%annot_info);
}
#####################################
sub pre_exon_annotation{
	my($gtf_file, $out_file)=@_;
	
	my %annot_info;
	open my($IN), "<", $gtf_file or die;
	while(<$IN>){#2
		my $gtf_pointer=sub_rna::read_gtf_line($_);
		my %gtf=%$gtf_pointer;
		if($gtf{type} eq 'exon'){#3
			$annot_info{$gtf{chr}}->{$gtf{exon_id}}->{gene_id} = $gtf{gene_id};
			$annot_info{$gtf{chr}}->{$gtf{exon_id}}->{gene_name}=$gtf{gene_name};
			#$annot_info{$gtf{chr}}->{$gtf{exon_id}}->{gene_biotype}=$gtf{gene_biotype};
			$annot_info{$gtf{chr}}->{$gtf{exon_id}}->{transcript_id}=$gtf{transcript_id};
			$annot_info{$gtf{chr}}->{$gtf{exon_id}}->{transcript_name}=$gtf{transcript_name};
			#$annot_info{$gtf{chr}}->{$gtf{exon_id}}->{transcript_biotype}=$gtf{transcript_biotype};
			$annot_info{$gtf{chr}}->{$gtf{exon_id}}->{start}=$gtf{start};
			$annot_info{$gtf{chr}}->{$gtf{exon_id}}->{end}=$gtf{end};
		}#3
	}#2
	close($IN);
	
	#export
	open my($OUT), ">", $out_file or die;
	print $OUT join("\t", 'chr',  'exon_id', 'start', 'end', 'transcript_id', 'transcript_name', 'gene_id', 'gene_name', 'overlapped'), "\n";
	foreach my $chr(keys %annot_info){#2
		my $pointer=$annot_info{$chr};
		my %hash=%$pointer;
		my $over_id_pointer=sub_rna::overlapped_genes(\%hash);
		my %over_id=%$over_id_pointer;
		foreach my $id(keys %hash){#3 #exon_id
			my $gene_id=$hash{$id}->{gene_id};
			my $gene_name=$hash{$id}->{gene_name} ? $hash{$id}->{gene_name} : 'NA';
			my $start=$hash{$id}->{start};
			my $end=$hash{$id}->{end};
			my $transcript_id=$hash{$id}->{transcript_id};
			my $transcript_name=$hash{$id}->{transcript_name} ? $hash{$id}->{transcript_name}  : 'NA';
			my $overlapped=(exists $over_id{$id}) ? $over_id{$id} : 'NA';
			$annot_info{$chr}->{$id}->{overlapped}=$overlapped;
			print $OUT join("\t", $chr, $id, $start, $end, $transcript_id, $transcript_name, $gene_id, $gene_name, $overlapped), "\n";
		}#3
	}#2
	close($OUT);
	
	return(\%annot_info);
}
#########################################
#analyze gtf file 
sub pre_genome_annotation{
	my($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	
	print "\tgene annotation\n";
	my %gene_info;
	my $annot_file=$variables{dir_result}.'gene_annotation.txt';
	if (-f $annot_file){
		my $gene_info_pointer=sub_common::file_to_multihash($annot_file, "\t");
		%gene_info=%$gene_info_pointer;
	}
	else{
		my $gene_info_pointer=sub_rna::pre_gene_annotation($variables{file_annot}, $annot_file);
		%gene_info=%$gene_info_pointer;
	}
	#assign genes regions
	my $pointer=sub_bioseq::chr_region_fragments(\%gene_info);
	%gene_info=%$pointer;
	
	print "\ttranscript annotation\n";
	my %transcript_info;
	$annot_file=$variables{dir_result}.'transcript_annotation.txt';
	if (-f $annot_file){
		my $transcript_info_pointer=sub_common::file_to_multihash($annot_file, "\t");
		%transcript_info=%$transcript_info_pointer;
	}
	else{
		my $transcript_info_pointer=sub_rna::pre_transcript_annotation($variables{file_annot}, $annot_file);
		%transcript_info=%$transcript_info_pointer;
	}
	#assign genes regions
	$pointer=sub_bioseq::chr_region_fragments(\%transcript_info);
	%transcript_info=%$pointer;
	
	print "\texon annotation\n";
	my %exon_info;
	$annot_file=$variables{dir_result}.'exon_annotation.txt';
	if (-f $annot_file){
		my $exon_info_pointer=sub_common::file_to_multihash($annot_file, "\t");
		%exon_info=%$exon_info_pointer;
	}
	else{
		my $exon_info_pointer=sub_rna::pre_exon_annotation($variables{file_annot}, $annot_file);
		%exon_info=%$exon_info_pointer;
	}
	$pointer=sub_bioseq::chr_region_fragments(\%exon_info);
	%exon_info=%$pointer;
	
	return(\%gene_info, \%transcript_info, \%exon_info);
	#return(\%exon_info);
}
########################################
#This function is common functions for reading sam file
#pass the function as parameters
sub analyze_sam_file{
	my($sam_file, $out_file, $analytic_func, $annot_pointer)=@_;
	my %annot_info=%$annot_pointer;
	
	my @sam_arr;
	my $head;
	my $num=1;
	open my($IN), "<", $sam_file or die "can't open $sam_file\n";
	open my($OUT), ">", $out_file or die "can't write to $out_file\n";
	#remove annottion line in sam file
	while(1){
		my $line=<$IN>;
		unless ($line=~/^\@/){
			my @items=split("\t", $line);
			my $single_pointer=sub_rna::read_sam_line($line);
			my %single=%$single_pointer;
			push(@sam_arr, \%single);
			$head=$items[0];
			last;
		}
	}
	#read sam line
	while (<$IN>){#2
		my @items=split("\t", $_);
		my $single_pointer=sub_rna::read_sam_line($_);
		my %single=%$single_pointer;
		unless ($head eq $items[0]){#3
			my $output=$analytic_func->(\@sam_arr, \%annot_info);
			print $OUT join("\t", $num, $head, $output), "\n";
			$head=$items[0];
			@sam_arr=();
			$num++;
			#print "$head\n";
		}#3
		push(@sam_arr, \%single);
	}#2
	close($IN);
	close($OUT);
	
}

#####################################################
#read Q values 
sub QC_Qvalues{ #1
	my ($raw_files_pointer, $out_file)=@_;
	my @raw_files=@$raw_files_pointer;
	#$variables{QC_compression}
	
	my @quals;
	my $num=0;
	open my($QC), ">", $out_file or die;
	foreach my $raw_file(@raw_files){#2
		print "\tCheck sequencing quality of $raw_file.\n";
		my $in_obj = Bio::SeqIO->new(-file => $raw_file, -format => 'fastq');
		while (my $seq = $in_obj->next_seq() ) {#3
				my $quality_pointer=$seq->Bio::Seq::Quality::qual;
				my @quality=@$quality_pointer;
				@quals= ($num==0) ? @quality : List::MoreUtils::pairwise {$a+$b} @quals, @quality;
				$num++;
				if ($num== 1000 or eof ){#4
						@quals=map{int($_/$num+0.5)} @quals;
						print $QC join("\t", @quals), "\n";
						$num=0;
						undef @quals;
				}#4
		}#3
	}#2
	close($QC);

}#1
###############################################
#
sub target_references{
	my($sam_arr_pointer, $annot_info_pointer)=@_;
	my @sam_arr=@$sam_arr_pointer;
	my %annot_info=%$annot_info_pointer;
	
	my @refs;
	foreach my $single_pointer(@sam_arr){#2
		my %single=%$single_pointer;
		my$chr=$single{ref};
		my @matches=split(',', $single{match_str});
		my $ids='NA';
		foreach my $se(@matches){#3
			my($frag_start, $frag_end)=split('-', $se);
			my $region=int($frag_start/1e5);
			#junction exon
			if (exists $annot_info{$chr}->{$region}){#4
				my $region_pointer=$annot_info{$chr}->{$region};
				my %region_annot_info=%$region_pointer;
				foreach my $id(keys %region_annot_info){#5
					my $start=$region_annot_info{$id}->{start};
					my $end=$region_annot_info{$id}->{end};
					my $gene_id=$region_annot_info{$id}->{gene_id};
					my $judging=sub_common::overlapping($start, $end, $frag_start, $frag_end);
					$ids=sub_basic::most_intersections($ids, $gene_id) if $judging eq 'in' or $judging eq 'same';
				}#5
			}#4
		}#3
		push(@refs, $ids) unless $ids eq 'NA';
	}#2
	
	my $references=sub_basic::most_unions(\@refs);
	#print "$references\n";
	return($references);
}
############################################
#analyze sam file and give saturation for estimation of  sequencing depth 
sub QC_saturation{
	my($target_file, $out_file)=@_;
	
	my %SD;
	my $ref1=0;
	my $ref10=0;
	my $ref20=0;
	open my($IN), "<", $target_file or die "can't open $target_file";
	open my($OUT), ">", $out_file or die "can't open $target_file";
	print $OUT join("\t", 'raw_num', 'RC1', 'RC10', 'RC20'), "\n";
	while(<$IN>){#2
		chomp($_);
		my($num, $query_info, $id)=split("\t", $_);
		if ($id cmp 'NA' and $id !~/\,/){#2
			$SD{$id}++;
			if ($SD{$id}==1){
				$ref1++;
				print $OUT join("\t", $num, $ref1, $ref10, $ref20), "\n";
			}
			elsif ($SD{$id}==10){
				$ref10++;
				print $OUT join("\t", $num, $ref1, $ref10, $ref20), "\n";
			}
			elsif ($SD{$id}==20){
				$ref20++;
				print $OUT join("\t", $num, $ref1, $ref10, $ref20), "\n";
			}
		}#2
	}#2
	close($IN);
	close($OUT);

}
############################################
#count number of reads
sub count_reads{
	my($target_file, $out_file)=@_;
	
	my (%RCmin, %RCmid, %RCmax, %multiple, %num_out);
	open my($IN), "<", $target_file or die "can't open $target_file";
	while(<$IN>){#2
		chomp($_);
		my($num, $query_info, $id_str)=split("\t", $_);
		if ($id_str eq 'NA'){#3
			$num_out{unmapped_raw_num}++;
		}#3
		else{#3
			my @ids=split(',', $id_str);
			if (@ids==1){#4
				$RCmin{$ids[0]}++;
				$RCmid{$ids[0]}++;
				$RCmax{$ids[0]}++;
				$num_out{unique_mapped_raw_num}++;
			}#4
			else{#4
				foreach my $id(@ids){
					$RCmax{$id}++;
				}
				$multiple{$id_str}++;
				$num_out{multiple_mapped_raw_num}++;
			}#4
		}#3
	}#2
	close($IN);
	
	#calcualte RCmid based on the ratios of RCmin
	foreach my $ids_str(keys %multiple){#2
		my @multiple_ids=split(',', $ids_str);
		my %RC_ratio;
		my $unique_total=0;
		foreach my $id(@multiple_ids){#3
			$RC_ratio{$id}=($RCmin{$id}) ? $RCmin{$id} : 1;
			$unique_total += $RC_ratio{$id};
		}#3
		#assign RC based on their ratio
		my $multiple_RC=$multiple{$ids_str};
		foreach my $id(keys %RC_ratio){#3
			my $assigned_RC=$multiple_RC*($RC_ratio{$id}/$unique_total);
			$RCmid{$id}=($RCmid{$id}) ? ($RCmid{$id}+$assigned_RC) : $assigned_RC;
		}#3
	}#2
	
	#export
	open my($OUT), ">",  $out_file or die "can't open $target_file";
	print $OUT join("\t", 'gene', 'RCmin', 'RCmid', 'RCmax'), "\n";
	foreach my $id(keys %RCmax){
		my $max=$RCmax{$id};
		my $mid=($RCmid{$id}) ? int($RCmid{$id} +0.5) : 0;
		my $min=($RCmin{$id}) ? $RCmin{$id} : 0;
		print $OUT join("\t", $id, $min, $mid, $max), "\n";
	}
	close($OUT);
	
	return(\%num_out);
}

############################################
#count number of reads
sub export_read_counts{
	my($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my @sample_names=split(',', $variables{sample_names});
	
	#get read counts of samples
	my(%RCmin, %RCmid, %RCmax);
	foreach my $sample_name(@sample_names){#2
		my $sample_dir=$sample_info{$sample_name}->{sample_dir};
		my $file=$sample_dir.$sample_name.'_read_counts.txt';
		open my($IN), "<", $file or die "can't open $file";
		my $first_line=<$IN>;
		while(<$IN>){
			chomp($_);
			my($id, $min, $mid, $max)=split("\t", $_);
			$RCmin{$id}->{$sample_name}=$min;
			$RCmid{$id}->{$sample_name}=$mid;
			$RCmax{$id}->{$sample_name}=$max;
		}
		close($IN);
	}#2
	
	#export
	sub_common::hash2_to_file(\%RCmin, $variables{dir_stat}.'minimum_read_counts.csv', ',');
	sub_common::hash2_to_file(\%RCmid, $variables{dir_stat}.'middle_read_counts.csv', ',');
	sub_common::hash2_to_file(\%RCmax, $variables{dir_stat}.'maximum_read_counts.csv', ',');
	#
}
##################################################################################
#only mismatch on the 5 end and im the middle of the read sequence
#3' end mismatched is not allowed due to the short sequences
#output:%adapters

sub adapters{#1
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;

	my (%adapters, @adapter_3_5end_back, @adapter_3_5end_trim, @adapter_5_3end_back, @adapter_5_3end_trim);
	#generate exactly matched adapters
	#3' end adapter without mismatch
	foreach (1..5){ #matched length is 11 nt
		push(@adapter_3_5end_back, substr($variables{adapter_3},$_,$variables{match_len}+3) );
	}
	foreach (4..7){
		push(@adapter_3_5end_trim, substr($variables{adapter_3},0,$_) );
	}
	#5' end adapter without mismatch
	foreach (1..5){
		my $str=substr($variables{adapter_5}, -$variables{match_len}-3-$_);
		push(@adapter_5_3end_back, substr($str, 0, $variables{match_len}+3) );
	}
	foreach (4..7){
		push(@adapter_5_3end_trim, substr($variables{adapter_5}, -$_) );
	}


	#generate mis-matched adapter 3
	my @mismatch_adapter_3_5end;
	# one mismatched among 9nt
	my $seq=substr($variables{adapter_3}, 0, $variables{match_len}+1);
	my $seq_len=length($seq);
	foreach (0..$seq_len-2){#2
		my $str=substr($seq, 0, $_).".".substr($seq, $_+1, $seq_len-$_);
		push( @mismatch_adapter_3_5end, $str);
		#print "$variables{adapter_3}:$str\n";
	}#2 
  
	#one gap allowed in the first 11nt
	$seq=substr($variables{adapter_3}, 0, $variables{match_len}+3);
	$seq_len=length($seq);
	foreach (1..$seq_len-2){#2
		my $str=substr($seq,0,$_).substr($seq, $_+1, $seq_len-$_);
		push(@mismatch_adapter_3_5end, $str);
		#print "$variables{adapter_3}:$str\n";
	}#2 
   
	#one insert allowed in the first 11nt
	foreach (1..$seq_len-2){#2
		my $str=substr($seq,0,$_).".".substr($seq, $_, $seq_len-$_);
		push(@mismatch_adapter_3_5end, $str);
		#print "$variables{adapter_3}:$str\n";
	}#2 
     
	#one mismatch of 10nt or 11nt when 5' end is backward sided 1nt or 2nt of adapter 3
	my @mismatch_adapter_3_5end_back;
	foreach my $matched_adapter(@adapter_3_5end_back){#2
		my $len=length($matched_adapter);
		foreach (0..$len-2){#3
			my $str=substr($matched_adapter,0,$_).".".substr($matched_adapter, $_+1, $len-$_);
			push(@mismatch_adapter_3_5end_back, $str);
			#print "$variables{adapter_3}:$str\n";
		}#3 
	}#2

	#generate mis-matched adapter 5
	my @mismatch_adapter_5_3end;
	# one mismatched among 9nt
	$seq=substr($variables{adapter_5}, -$variables{match_len}-1);
	$seq_len=length($seq);
	foreach (0..$seq_len-2){#2
		my $str=substr($seq, 0, $_).".".substr($seq, $_+1, $seq_len-$_);
		push( @mismatch_adapter_5_3end, $str );
		#print "$variables{adapter_5}:$str\n";
	}#2 

	#one gap allowed in the first 11nt
	$seq=substr($variables{adapter_5}, -$variables{match_len}-3);
	$seq_len=length($seq);
	foreach (1..$seq_len-2){#2
		my $str=substr($seq,0,$_).substr($seq, $_+1, $seq_len-$_) ;
		push(@mismatch_adapter_5_3end, $str);
		#print "$variables{adapter_5}:$str\n";
	}#2 
   
	#one insert allowed in the first 11nt
	foreach (1..$seq_len-2){#2
		my $str=substr($seq,0,$_).".".substr($seq, $_, $seq_len-$_);
		push(@mismatch_adapter_5_3end, $str);
		#print "$variables{adapter_5}:$str\n";
	}#2 

	#one mismatch of 10nt or 11nt when 3' end is backward sided 1nt or 2nt of adapter 5
	my @mismatch_adapter_5_3end_back;
	foreach my $matched_adapter(@adapter_5_3end_back){#2
		my $len=length($matched_adapter);
		foreach (0..$len-2){#3
			my $str=substr($matched_adapter,0,$_).".".substr($matched_adapter, $_+1, $len-$_);
			push(@mismatch_adapter_5_3end_back, $str);
			#print "$variables{adapter_5}:$str\n";
		}#3 
	}#2

	#Note: mismatch of adapter 5 with 5' end of read and of adapter 3 with 3' end of read is not allowed 
	#used for sequence truncation
	$adapters{adapter_3_5end}=substr($variables{adapter_3}, 0, $variables{match_len});
	$adapters{adapter_5_3end}=substr($variables{adapter_5}, -$variables{match_len});
	$adapters{adapter_3_5end_back}=\@adapter_3_5end_back;
	$adapters{adapter_3_5end_trim}=\@adapter_3_5end_trim; 
	$adapters{adapter_5_3end_back}=\@adapter_5_3end_back; 
	$adapters{adapter_5_3end_trim}=\@adapter_5_3end_trim;
	$adapters{mismatch_adapter_3_5end}=\@mismatch_adapter_3_5end;
	$adapters{mismatch_adapter_5_3end}=\@mismatch_adapter_5_3end;
	$adapters{mismatch_adapter_3_5end_back}=\@mismatch_adapter_3_5end_back;
	$adapters{mismatch_adapter_5_3end_back}=\@mismatch_adapter_5_3end_back;        
	return(\%adapters);
}#1


##################################################################
#truncate sequence of adapter from read
sub adapter_3_truncation{#1
  my $variables_pointer=$_[0];
  my %variables=%$variables_pointer;
  my $seq=$_[1];
  $seq=~s/^N|N$//g;
  my $adapters_pointer=$_[2];
  my %adapters=%$adapters_pointer;
  my $adapter_3_5end_back_pointer=$adapters{adapter_3_5end_back};
  my @adapter_3_5end_back=@$adapter_3_5end_back_pointer;
  my $adapter_3_5end_trim_pointer=$adapters{adapter_3_5end_trim};
  my @adapter_3_5end_trim=@$adapter_3_5end_trim_pointer;
  my $mismatch_adapter_3_5end_pointer=$adapters{mismatch_adapter_3_5end};
  my @mismatch_adapter_3_5end=@$mismatch_adapter_3_5end_pointer;
  my $mismatch_adapter_3_5end_back_pointer=$adapters{mismatch_adapter_3_5end_back};
  my @mismatch_adapter_3_5end_back=@$mismatch_adapter_3_5end_back_pointer;

  my $seq_trunc='NA';
  my $found=0;
  if ($seq=~/^$adapters{adapter_3_5end}/){#2
    $found=1;
    #print "exact_match:$variables{adapter_3}:$adapters{adapter_3_5end}\n$seq\n($seq_trunc)\n\n\n";
  }#2
  elsif ($seq=~s/$adapters{adapter_3_5end}/#/){#2
    my @array=split("#", $seq); 
    $seq_trunc=$array[0];
    $found=1;
    #print "exact_match:$variables{adapter_3}:$adapters{adapter_3_5end}\n$_[1]\n$seq\n($seq_trunc)\n\n\n";
  }#2
  else{#2     #3' end of read matching
    foreach my $end5(@adapter_3_5end_trim){#3
      if ($seq eq $end5){#4
        $found=1;
    #print "3' adapter trim:$variables{adapter_3}:$adapters{adapter_3_5end}\n$_[1]\n$seq\n($seq_trunc)\n\n\n";
        last;
      }#4
      elsif ($seq=~s/$end5$/#/){#4
        my @array=split("#", $seq); 
        $seq_trunc=$array[0];
        $found=1;
    #print "3' adapter trim:$variables{adapter_3}:$adapters{adapter_3_5end}\n$_[1]\n$seq\n($seq_trunc)\n\n\n";
        last;
      }#4
    } #3
  }#2
 
  #read sequence is all sequences of adapter 3'
  if ($found==0){#2
    foreach my $back(@adapter_3_5end_back){#3
      if ($seq=~/^$back/){#4
        $found=1;
    #print "all 3'adapter:$variables{adapter_3}:$adapters{adapter_3_5end}\n$_[1]\n$seq\n($seq_trunc)\n\n\n";
        last;
      }#4   
    }#3
  }#2


  #one mismatch allowed
  if ($variables{mismatch_allowed} eq "yes" and $found==0){#2        mismatch circyling
    foreach my $mismatch(@mismatch_adapter_3_5end){#3
      if ($seq=~/^$mismatch/ ){#4
        $found=1;
    #print "miss match:$variables{adapter_3}:$adapters{adapter_3_5end}\n$seq\n($seq_trunc)\n\n\n";
      }#4   
      elsif ($seq=~s/$mismatch/#/ ){#4
        my @array=split("#", $seq); 
        $seq_trunc=$array[0];
        $found=1;
    #print "miss match:$variables{adapter_3}:$adapters{adapter_3_5end}\n$_[1]\n$seq\n($seq_trunc)\n\n\n";
      }#4   
    } #3
 
    #read sequence is all sequences of adapter 3'
    if ($found==0){#3
      foreach my $back(@mismatch_adapter_3_5end_back){#4
        if ($seq=~/^$back/){#5
          $found=1;
          #print "mis_match:$variables{adapter_3}:$adapters{adapter_3_5end}\n$_[1]\n$seq\n($seq_trunc)\n\n\n";
          last;
        }#5   
      }#4
    }#3
  }#2                 mismatch circyling

  #print "match:\t$variables{adapter}\n$seq_original\n$seq_trunc\n\n";
  return($seq_trunc);
}#1
 


#######################################################################################
sub adapter_5_truncation{#1
  my $variables_pointer=$_[0];
  my %variables=%$variables_pointer;
  my $seq=$_[1];
  $seq=~s/^N|N$//g;
  my $adapters_pointer=$_[2];
  my %adapters=%$adapters_pointer;
  my $adapter_5_3end_back_pointer=$adapters{adapter_5_3end_back};
  my @adapter_5_3end_back=@$adapter_5_3end_back_pointer;
  my $adapter_5_3end_trim_pointer=$adapters{adapter_5_3end_trim};
  my @adapter_5_3end_trim=@$adapter_5_3end_trim_pointer;
  my $mismatch_adapter_5_3end_pointer=$adapters{mismatch_adapter_5_3end};
  my @mismatch_adapter_5_3end=@$mismatch_adapter_5_3end_pointer;
  my $mismatch_adapter_5_3end_back_pointer=$adapters{mismatch_adapter_5_3end_back};
  my @mismatch_adapter_5_3end_back=@$mismatch_adapter_5_3end_back_pointer;

    #print "Before++++:$variables{adapter_5}:$adapters{adapter_5_3end}\n$seq\n\n\n";

  my $seq_trunc="NA";
  my $found=0;
  if($seq=~/$adapters{adapter_5_3end}$/ or $seq eq "NA"){#2
    $found=1;
    #print "all_adapter:\t$variables{adapter_5}:$adapters{adapter_5_3end}\n$seq\n($seq_trunc)\n\n" ;
  }#2
  elsif($seq=~s/$adapters{adapter_5_3end}/#/){#2
    my @array=split("#", $seq); 
    $seq_trunc=$array[-1];
    $found=1;
    #print "exact_match:\t$variables{adapter_5}:$adapters{adapter_5_3end}\n$seq\n$array[-1]\n($seq_trunc)\n\n" ;
  }#2
  else{#2     #5' end of read matching
    foreach my $end3(@adapter_5_3end_trim){#3
      if ($seq eq $end3){#4
        $found=1;
        #print "5 end match:\t$variables{adapter_5}:$end3\n$_[-1]\n$seq\n($seq_trunc)\n\n";
        last;
      }#4
      elsif ($seq=~s/^$end3/#/){#4
        my @array=split("#", $seq); 
        $seq_trunc=$array[-1];
        $found=1;
        #print "5 end match:\t$variables{adapter_5}:$end3\n$_[-1]\n$seq\n($seq_trunc)\n\n";
        last;
      }#4
    } #3  
  }#2
 
  #read sequence is all sequences of adapter 5'
  if ($found==0){#2
    foreach my $back(@adapter_5_3end_back){#3
      if ($seq=~/$back$/){#4
        $found=1;
        #print "all adapter 5':\t$variables{adapter_5}:$back\n$seq\n($seq_trunc)\n\n\n";
        last;
      }#4   
    }#3
  }#2


  #one mismatch allowed
  if ($variables{mismatch_allowed} eq "yes" and $found==0){#2        mismatch circyling
    foreach my $mismatch(@mismatch_adapter_5_3end){#3
      if ($seq=~/$mismatch$/){#4
        $found=1;
        #print "one mis of all adapter:\t$variables{adapter_5}:$mismatch\n$seq\n($seq_trunc)\n\n" ;
      } #4
      elsif ($seq=~s/$mismatch/#/ ){#4
        my @array=split("#", $seq); 
        $seq_trunc=$array[-1];
        $found=1;
        #print "one mis:\t$variables{adapter_5}:$mismatch\n$_[-1]\n$seq\n($seq_trunc)\n\n" ;
      } #4
    }#3

    #read sequence is all sequences of adapter 5'
    if ($found==0){#3
      foreach my $back(@mismatch_adapter_5_3end_back){#4
        if ($seq=~/$back$/){#5
          $found=1;
          #print "one mis:\t$variables{adapter_5}:$back\n$seq\n($seq_trunc)\n\n";
          last;
        }#5   
      }#4
    }#3
  }#2                 mismatch circyling

  #print "after:\t$seq\n$seq_trunc\n\n";
  return($seq_trunc);
}#1

####################################################
#read SAM file determined by bowtie1 and bowtie2
#1S12_GGCTAC:4316023_1   0       precursor_miRNA:hsa-mir-451a    17      255     18M     *       0       0       AGACCGTTACCATTACTG      IIIIIIIIIIIIIIIIII      XA:i:1  MD:Z:1A16       NM:i:1
sub read_bowtie1_sam{
	my ($fasta_file, $alignment_output)=@_;

	my (%alignment_out, %alignment, %un_alignment, %multiple_alignment, %mappable_hits, %mappable_counts, %SD_counts);
	my %num_out=(query_reads_num=>0, nr_query_reads_num=>0, 
				mappable_reads_num=>0, nr_mappable_reads_num=>0, 
				unmappable_reads_num=>0, nr_unmappable_reads_num=>0, 
				multiple_mappable_reads_num=>0, nr_multiple_mappable_reads_num=>0, 
				mappable_ref_num=>0, );
	my $n=1;
	open my($IN), "<", $alignment_output or die;
	while (my $line=<$IN>){#2 file circyling
		unless ($line=~/^@/){#3
			chomp($line);
			my @items=split("\t", $line);
			my $query_name=$items[0];
			my ($sample_name, $query_str)=split(':', $query_name);
			my ($seq_No, $query_read_counts)=split('_', $query_str);
			my $flags=$items[1];
			my $ref_name=$items[2];
			my $ref_offset=$items[3];
			my $read_seq=$items[9];
			my %optional_fields;
			foreach (@items){
				if ($_=~/:i:|:Z:/){
					my ($a, $b, $c)=split(":", $_);
					$optional_fields{$a}=$c;
				} 
			}
			
			#update %mappable_hits
			if( exists $mappable_hits{$query_name} ){
				$mappable_hits{$query_name} ++;
			}
			else{
				$mappable_hits{$query_name}=1;
				$num_out{query_reads_num} += $query_read_counts;
				$num_out{nr_query_reads_num} ++;
			}

			#alignment, un_alignment, and multiple_alignment
			if( exists $optional_fields{XM} and $optional_fields{XM}==0 ){	#4 #no alignment
				$num_out{unmappable_reads_num} += $query_read_counts;
				$num_out{nr_unmappable_reads_num} ++;
				$un_alignment{$query_name}->{RC}=$query_read_counts;
				$un_alignment{$query_name}->{seq}=$read_seq;
			}#4
			elsif( exists $optional_fields{XM} and $optional_fields{XM}>0 ){#4 #multiple alignment
				$multiple_alignment{$query_name}->{RC}=$query_read_counts;
				$multiple_alignment{$query_name}->{seq}=$read_seq;
				$num_out{multiple_mappable_reads_num} += $query_read_counts;
				$num_out{nr_multiple_mappable_reads_num} ++;
			}#4
			else{#4 unique alignment
				if ($flags==0)		{	$alignment{$n}->{mapping_strand}= 'forward';	}
				elsif ($flags==16)	{	$alignment{$n}->{mapping_strand}= 'reverse';		}
				#update %alignment
				$alignment{$n}->{sam_line}=$line;
				$alignment{$n}->{query_name}=$query_name;
				$alignment{$n}->{seq_No}=$seq_No;
				$alignment{$n}->{read_counts}=$query_read_counts;
				$alignment{$n}->{ref_name}=$ref_name;
				$alignment{$n}->{ref_offset}=$ref_offset;
				$alignment{$n}->{read_seq}=$read_seq;
				#update %mappable_counts
				if (exists $mappable_counts{$ref_name}->{seq_info}){
					$mappable_counts{$ref_name}->{seq_info} .= ';'.$query_str.','.$read_seq;
					$mappable_counts{$ref_name}->{query_name_str} .= ','.$query_name;
				}
				else{	
					$mappable_counts{$ref_name}->{seq_info} =$query_str.','.$read_seq;
					$mappable_counts{$ref_name}->{query_name_str} =$query_name;
				}
				#update %SD_counts
				if(exists $SD_counts{$seq_No}) {
					$SD_counts{$seq_No} .= ';;'.$ref_name;
				}
				else{
					$SD_counts{$seq_No} = $ref_name;
				}
				#update %num_out
				if($mappable_hits{$query_name}==1){
					$num_out{mappable_reads_num} += $query_read_counts;
					$num_out{nr_mappable_reads_num} ++;
				}
			}#4
		}#3
		$n++;
	}#2 file circyling
	close($IN);
	
	#generate RC of unique alignment
	foreach my $ref_name(keys %mappable_counts){
		my @query_names=split(',', $mappable_counts{$ref_name}->{query_name_str});
		$mappable_counts{$ref_name}->{top_RC}=0;
		$mappable_counts{$ref_name}->{middle_RC}=0;
		$mappable_counts{$ref_name}->{bottom_RC}=0;
		foreach my $query_name(@query_names){
			my $query_hits=$mappable_hits{$query_name};
			my ($sample_name, $query_str)=split(':', $query_name);
			my ($seq_No, $query_read_counts)=split('_', $query_str);
			$mappable_counts{$ref_name}->{top_RC} += $query_read_counts;
			$mappable_counts{$ref_name}->{middle_RC} += int(($query_read_counts*100/$query_hits)+0.5)/100;
			$mappable_counts{$ref_name}->{bottom_RC} += $query_read_counts if $query_hits==1;
		}
	}
	$num_out{mappable_ref_num}=keys %mappable_counts;
	
	$alignment_out{alignment_pointer}=\%alignment;
	$alignment_out{un_alignment_pointer}=\%un_alignment;
	$alignment_out{multiple_alignment_pointer}=\%multiple_alignment;
	$alignment_out{mappable_counts_pointer}=\%mappable_counts; #only the counting of unique mappable reads
	$alignment_out{SD_counts_pointer}=\%SD_counts;
	$alignment_out{num_out_pointer}=\%num_out;
	return(\%alignment_out);
}

#####################################
#export SD after read sam file
sub export_sequncing_depth{
	my($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $SD_counts_pointer=$variables{SD_counts_pointer};
	my %SD_counts=%$SD_counts_pointer;
	my $index_name='xx';
	
	#export sequencing depth(SD)
	my @ref_names;
	my $ref_num=0;
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
		my $raw_reads_num=read_log($variables{sample_log}, 'raw_reads_num');
	print $SD "$raw_reads_num\t", "$ref_num\n" if $raw_reads_num > $seq_No_max;
	close($SD);
}
###########################################################
#return the index of sample names used for comparison based on parameter
sub sample_index{
  my $variables_pointer=$_[0];
  my %variables=%$variables_pointer;
  my $sample_info_pointer=$variables{sample_info_pointer};
  my @sample_info=@$sample_info_pointer;
  my @index_parameter=split(";", $_[1]);  

  my %sample_all;
  my @index=split("\t", shift @sample_info);
  for (my $i=0; $i<@index; $i++){#2
    foreach (@sample_info){#3
      my @array=split("\t", $_);
      $sample_all{$array[0]}->{$index[$i]}=$array[$i];
    }#3
  }#2

  my %sample_indexed=%sample_all;
  foreach (@index_parameter){#2
    my($index_name, $para)=split("=", $_);
    foreach my $sample_name(keys %sample_all){#3
      if ($para=~/_/){
        my($para_1, $para_2)=split("_", $para);
        delete $sample_indexed{$sample_name} unless $sample_all{$sample_name}->{$index_name}>=$para_1 and $sample_all{$sample_name}->{$index_name}<=$para_2;
      }
      else{
        delete $sample_indexed{$sample_name} unless $sample_all{$sample_name}->{$index_name} eq $para;  
      } 
    }#3
  }#2

  my @samples;
  foreach (keys %sample_indexed){
    push (@samples, $_);
  }
  return(\@samples);
}
########################################
sub copy_bowtie_index{
	my ($variables_pointer, $version)=@_;
	my %variables=%$variables_pointer;
	
	my %index_info;
	$index_info{bowtie1}->{tails}=join(',', '.1.ebwt', '.2.ebwt','.3.ebwt','.4.ebwt','.rev.1.ebwt','.rev.2.ebwt');
	$index_info{bowtie1}->{dir}=$variables{dir_bowtie1};
	$index_info{bowtie1}->{script}=$variables{dir_bowtie1}.'bowtie-build';
	$index_info{bowtie2}->{tails}=join(',', '.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2');
	$index_info{bowtie2}->{dir}=$variables{dir_bowtie2};
	$index_info{bowtie2}->{script}=$variables{dir_bowtie2}.'bowtie2-build';
	
	#copy bowtie index
	my @tails=split(',', $index_info{$version}->{tails});
	foreach my $tail(@tails){
		my $index_file=$variables{genome_index}.$tail;
		system("cp $index_file $variables{dir_bowtie_index}");
	}
	
	#copy fasta file
	system("cp $variables{genome_fasta_file}  $variables{dir_bowtie_index}" );
}
##################
sub refresh_cforest_line{
	my($file, $chr1, $new)=@_;
	
	print length($chr1), "\n";
	my @lines;
	open my($IN), "<", $file or die;
	while(<$IN>){
		chomp($_);
		if ($_=~/$chr1/){
			push(@lines, $new);
			print $_,"\n", $new, "\n";
		}
		else{
			push(@lines, $_);
		}
	}
	close($IN);

	open my($OUT), ">", $file or die;
	foreach(@lines){
		print $OUT "$_\n";
	}
	close($OUT);
}


########################################
1;  # make sure the file returns true or require will not succeed!#






