#! /usr/bin/perl -w
use strict;
use warnings;
use Bio::SeqIO;
use List::MoreUtils;
use List::Util;
use File::Find;

require "func_basic.pm"; #sub_basic::
require "func_common.pm"; #sub_common::

####################################################
#
#the file constains all subroutines required 
package sub_bioseq;

##############################


########################
####################
sub translate_DNA{
		my($dna)=@_;
		$dna=uc $dna;
		my %codons=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'.',
						'TAG'=>'.','TGC'=>'C','TGT'=>'C','TGA'=>'.','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P',
						'CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I',
						'ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K',
						'AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A',
						'GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');
		my @aa;
		for (my $i=0; $i<length($dna); $i+=3){#2
				my $g=substr($dna, $i, 3);
				#print $i;
				if( $codons{$g} ) {
					push(@aa, $codons{$g});
				}
				else{
					push(@aa, '-');
				}
		}#2
		my $aa_seq=join('', @aa);
		return($aa_seq);
}
###########################
sub format_DNA{
	my($DNA)=@_;
	
	$DNA=~tr/a-z/A-Z/;
	$DNA=~s/[^A-Z]//g;
	$DNA=~s/U/T/g;
	return($DNA)
}
##########################################
sub reverse_complement {
	my $dna_seq = $_[0];

	# reverse the DNA sequence
	my $revcom_seq = reverse($dna_seq);
	# complement the reversed DNA sequence
	$revcom_seq =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;

	return ($revcom_seq);
}
##################################################################
#truncate sequence of adapter from read
sub truncate_end3_seq{#1
	my ($adapter, $seq)=@_;
	my $adapter_len=length($adapter);
	my $seq_len=length($seq);
	
	my $seq_trunc=$seq;
	unless ($seq_len<=10 or $seq_trunc=~s/$adapter.*//){#2  all adapter sequence
		my $end_cut=($seq_len>$adapter_len) ? $adapter_len: $seq_len;
		for (my$i=$end_cut; $i>4; $i--){
			my $adapter_index=substr($adapter, 0, $i);
			last if $seq_trunc=~s/$adapter_index$//;
		}
	}#2

	return($seq_trunc);
}#1
################
#split hybrid sequence by enzyme site
sub split_hybrid_seq{
	my($sequence, $known_seq, $string, $type)=@_;
	
	my $trunc_seq=$sequence;
	if($known_seq=~/$sequence/){#2
		$trunc_seq='NA';
	}#2
	elsif($type eq 'head'){#2
		my $pos = 0;
		for (my $pos = index($sequence, $string, $pos); $pos >= 0; $pos = index($sequence, $string, $pos) ) {
				my $head_seq=substr($sequence, 0, $pos+length($string));
				my $tail_seq=substr($sequence, $pos );
				#print "$pos:$head_seq\n$tail_seq\n\n";
				if($known_seq=~/$tail_seq/ and length($tail_seq)>10 ){
					$trunc_seq=$head_seq;
					last;
				}
				$pos += length $string;
		}
	}#2
	elsif($type eq 'tail'){#2
		my $pos = length($sequence)-1;
		for ( my $pos = rindex($sequence, $string, $pos); $pos >= 0; $pos = rindex($sequence, $string, $pos) ) {
			my $head_seq=substr($sequence, 0, $pos+length($string));
			my $tail_seq=substr($sequence, $pos );
			 if($known_seq=~/$head_seq/ and length($head_seq)>10){
				$trunc_seq=$tail_seq;
				#print "$known_seq\n$pos:$head_seq\n$tail_seq\n\n";
				last;
			}
			$pos -= length $string;
		}
	}#2
	
	return($trunc_seq);
}
#####################
#split sequence by enzyme sites
sub split_seq{
	my ($seq, $enzyme_site)=@_;
	
	my $head_seq=$seq;
	my $tail_seq='NA';
	my @a=split(/$enzyme_site/, $seq);
	if (@a>=2){
		$tail_seq=pop @a;
		$head_seq=join($enzyme_site, @a);
		$head_seq .= $enzyme_site;
	}
	
	return($head_seq);
}
##########################
#get genome sequences from fasta file
sub genome_sequences{
	my($genome_fasta_file)=@_;
	
	my %genome;
	my $in_obj = Bio::SeqIO->new(-file => $genome_fasta_file, -format => 'fasta');
	while (my $seq_obj = $in_obj->next_seq() ) {
		my $chr=$seq_obj->display_id;
		my $seq=$seq_obj->seq;
		my $seq_len=$seq_obj->length;
		$genome{$chr}=$seq;
		print "\t$chr......\n";
	}
	
	return(\%genome);
}

##########################
#get sequences from fasta file
sub fasta_to_hash{
	my($fasta_file)=@_;
	
	my %hash;
	my $in_obj = Bio::SeqIO->new(-file => $fasta_file, -format => 'fasta');
	while (my $seq_obj = $in_obj->next_seq() ) {
		my $display_id=$seq_obj->display_id;
		my $seq=$seq_obj->seq;
		my $seq_len=$seq_obj->length;
		$hash{$display_id}=$seq;
	}
	return(\%hash);
}

###################################
#
sub numeric_seq{
	my $seq=$_[0];
	
	$seq=uc($seq);
	$seq=~s/A/1/g;
	$seq=~s/T/2/g;
	$seq=~s/G/3/g;
	$seq=~s/C/4/g;
	$seq=~s/[A-Z]/0/g;
	my @seq_arr=split("", $seq);
	my @power=(0..(@seq_arr-1));
	@power=reverse @power;
	@power=map {5**$_} @power;
	my @numeric_seq_arr=List::MoreUtils::pairwise {our($a, $b); $a*$b} @seq_arr, @power;
	my $numeric_seq=List::Util::sum(@numeric_seq_arr);
	
	return($numeric_seq);
}
#####################
sub read_fasta{
	my($fasta_file)=@_;
	
	my $in_obj = Bio::SeqIO->new(-file => $fasta_file, -format => 'fasta');
	while (my $seq_obj = $in_obj->next_seq() ) {#3
		my $displayid=$seq_obj->display_id();
		my $seq_len=$seq_obj->length();
		print "$displayid=$seq_len\n";
	}#3
	
}
#################################################################################
#used for conversion between %hash and fasta file
sub hash_to_fasta{
	my $hash_pointer=$_[0];
	my %hash=%$hash_pointer;
	my $fasta_file=$_[1];
  
	open my($OUT), ">", $fasta_file or die;
	foreach my $display_id(sort (keys %hash)){
		my $read_seq=$hash{$display_id}->{read_seq};
		if (exists $hash{$display_id}->{ref_name}) { print $OUT ">", "$display_id|$hash{$display_id}->{ref_name}\n", "$read_seq\n"; }
		else {  print $OUT ">", "$display_id\n", "$read_seq\n";  }
	}
	close($OUT);
}



###################################
#use for regioning fragment within chromosome 
#key1 $chr
#key2 id
sub chr_region_fragments{
	my($pointer, $ratio)=@_;
	my %hash=%$pointer;
	$ratio=1e5 unless $ratio;
	
	my %region_hash;
	foreach my $key1(keys %hash){#2
		my $pointer1=$hash{$key1};
		my %hash1=%$pointer1;
		foreach my $key2(keys %hash1){#3
			my $pointer2=$hash1{$key2};
			my %hash2=%$pointer2;
			my $start_region=int($hash1{$key2}->{start}/$ratio);
			my $end_region=int($hash1{$key2}->{end}/$ratio)+1;
			foreach my $region($start_region..$end_region){#4
				foreach my $key3(keys %hash2){#5
					$region_hash{$key1}->{$region}->{$key2}->{$key3}=$hash2{$key3};
					#print "$key1:$region:$key2:$key3:$hash2{$key3}\n";
				}#5
			}#4
		}#3
		#print "\n";
	}#2
	return(\%region_hash);
}
###################################
#use for regioning fragment within chromosome 
#key=$id
sub region_fragments{
	my($pointer, $ratio)=@_;
	my %hash=%$pointer;
	$ratio=1e5 unless $ratio;
	
	my %region_hash;
	foreach my $key(keys %hash){#2
		my$pointer=$hash{$key};
		my %hash1=%$pointer;
		my $start_region=int($hash{$key}->{start}/$ratio);
		my $end_region=int($hash{$key}->{end}/$ratio);
		foreach my $region($start_region..$end_region){#3
			foreach my $key1(keys %hash1){#4
				$region_hash{$region}->{$key}->{$key1}=$hash1{$key1};
				#print "
			}#4
		}#3
			
	}#2
	return(\%region_hash);
}
##############################
#return shared sequencing with two equal length sequences
sub shared_seq2{
	my($a, $b)=@_;
	my @a = split //, $a;
	my @b = split //, $b;
	my $mis=0;
	my $new_b = '';
	#printf ("%s\n%s\n", $a, $b);
	for(my $i = 0; $i < scalar(@a); $i++) {#2
		if ($a[$i] eq $b[$i]){
			$new_b .=  $b[$i];
		}
		else{
			$new_b .= '-';
		$mis++;
		}
	}#2
	#print "$mis: $new_b\n";
	return($new_b);
}
##############################
#return shared sequences among multiple equal length sequences
sub shared_seq{
	my $pointer=shift @_;
	my @seq_arr=@$pointer;
	my $share=shift @seq_arr;
	#
	foreach my $new(@seq_arr){
		$share=shared_seq2($new, $share);
	}
	return($share);
}
##############################################33
#get shared sequencing, and return number of mismatches
sub seq_mismatch{
	my $pointer=shift @_;
	#get shared sequence
	my $shared_seq=shared_seq($pointer);
	#
	my @matched=($shared_seq=~/-/g);
	my $mis=@matched;
	return($mis);
}

################################
#split FASTQ file into mulitple parts
sub split_fastq_file{
	my ($in_dir, $out_dir, $lines_num)=@_;
	$out_dir=sub_common::format_directory($out_dir);
	$lines_num //= 16000000; # the file size is 1GB
	
	#achieve all fastq files
	my $fastq_files_pointer=sub_common::files_list($in_dir, 'file');
	my @fastq_files=grep {$_=~/\.fq$|\.fastq$/} @$fastq_files_pointer;
	
	#split fastq file
	foreach my $fq_file(@fastq_files){
		print "Split $fq_file\n";
		my $fq_head=sub_common::file_operation($fq_file, 'name_head');
		$fq_head .= '_' unless $fq_head=~/\_$/;
		my $out_fq_head=$out_dir.$fq_head;
		#print "split -d -a 3 -b 1G $fq_file $fq_head\n";
		system("split -d -a 3 -l $lines_num $fq_file $out_fq_head");
	}

	#rename fastq files
	my $out_files_pointer=sub_common::files_list($out_dir, 'file');
	foreach my $file_head(@$out_files_pointer){
		unless ($file_head=~/\.fq$|\.fastq$/){
			my $out_file=$file_head.'.fq';
			system("mv $file_head $out_file");
		}
	}
	
	#
}

###########################33
1;