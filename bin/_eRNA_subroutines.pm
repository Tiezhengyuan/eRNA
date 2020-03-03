#! /usr/bin/perl -w
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq::Quality;
use List::MoreUtils;
use List::Util;
use File::Find;
#use Carp;
#use Venn::Chart;

#the file constains all subroutines required for running E_RNA_pipeline

package E_RNA;

################################################################
#get date and time
sub get_time{
	my $begin_pointer=$_[0];
	my @begin=@$begin_pointer;
	my $end_pointer=$_[1];
	my @end=@$end_pointer;
	
	my ($sec1,$min1,$hour1,$monthday1,$month1,$year1,$weekday1,$yearday1,$isdaylight1)=@begin;
	my ($sec2,$min2,$hour2,$monthday2,$month2,$year2,$weekday2,$yearday2,$isdaylight2)=@end;
	my $duration=($yearday2-$yearday1)*24*60*60 + ($hour2-$hour1)*60*60 + ($min2-$min1)*60 + ($sec2-$sec1); #second
	my $duration_m = int($duration/60);  #minutes
	my $duration_s = $duration % 60;  #seconds
	$year1 += 1900;
	$year2 += 1900;
	$month1 ++;
	$month2 ++;
	my $begin_time=$year1."/".$month1."/".$monthday1.", ".$hour1.":".$min1;
	my $end_time=$year2."/".$month2."/".$monthday2.", ".$hour2.":".$min2;
	my $my_time=$begin_time."-----".$end_time.". Duration: ".$duration_m."min".$duration_s."sec";
	return($my_time);
}

###################################################
##intersection of two arrays
sub two_intersection{#1
  	my ($a_pointer, $b_pointer)=@_;
	my @a=@$a_pointer;
	my @b=@$b_pointer;

	my %out=(a_num=>0, b_num=>0, a_only_num=>0, b_only_num=>0, ab_only_num=>0,  total_num=>0); 
	$out{a_num}=@a;
	$out{b_num}=@b;
	my @ab=(@a, @b);
	@ab=grep($_, @ab);
	@ab=List::MoreUtils::uniq @ab;
	$out{total_num}=@ab;
	foreach my $item(@ab) {#3
		if (List::Util::first {$item eq $_} @a and List::MoreUtils::none {$item eq $_} @b){ $out{a_only_num}++; }
		if (List::Util::first {$item eq $_} @a and List::Util::first {$item eq $_} @b){ $out{ab_only_num}++; }
		else{  $out{b_only_num}++; }
	}#3
	
	return(\%out);
}#1
#############
#intersection of three arrays
sub three_intersection{#1
	my ($a_pointer, $b_pointer, $c_pointer)=@_;
	my @a=@$a_pointer;
	my @b=@$b_pointer;
	my @c=@$c_pointer;
	my %out=(a_num=>0, b_num=>0, c_num=>0, a_only_num=>0, b_only_num=>0, c_only_num=>0, 
           ab_only_num=>0, ac_only_num=>0, bc_only_num=>0, abc_only_num=>0, total_num=>0); 

	$out{a_num}=@a;
	$out{b_num}=@b;
	$out{c_num}=@c;
	my @abc=(@a, @b, @c);
	@abc=grep($_, @abc);
	@abc=List::MoreUtils::uniq @abc;
	$out{total_num}=@abc;
	foreach my $item(@abc) {#3
		if (List::Util::first {$item eq $_} @a) {#4
			if (List::Util::first {$item eq $_} @b and List::Util::first {$item eq $_} @c){
				$out{abc_only_num}++;
			}
			elsif (List::Util::first {$item eq $_} @b and List::MoreUtils::none {$item eq $_} @c){
				$out{ab_only_num}++;
			}
			elsif (List::Util::first {$item eq $_} @c and List::MoreUtils::none {$item eq $_} @b){
				$out{ac_only_num}++;
			}
			else{  $out{a_only_num}++; }
		}#4
		elsif (List::Util::first {$item eq $_} @b) {#4
			if ( List::Util::first {$item eq $_} @c){	$out{bc_only_num}++;	}
			else{  $out{b_only_num}++; }
		}
		else{  
			$out{c_only_num}++; 
			#print "$item\n";
		}
	}#3

	return(\%out);
}#1


################################################################
#group files based on the threads number
sub thread_files{#1
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my $files_pointer=$_[1];
	my @files=@$files_pointer;
	
	my (%files_hash,$Thread_No);
	foreach my $num(0..@files-1){#2
	
		my $thread_no = $num % $variables{threads_num} + 1;
		if($files_hash{$thread_no}){#3
			$files_hash{$thread_no} .=",".$files[$num];
		}
		else{
			$files_hash{$thread_no}=$files[$num];
		}#3
		$num++;
	}#2
    
	while(my($key,$value)=each(%files_hash)){
		#print "($key:$value)\n";
	}
  
	return(\%files_hash);  
}#1

###########################################
#get files information
# $type: 'first_line', 'lines_number', 'bytes', 'KB', 'MB'
sub file_info{
	my ($file, $type)=@_;
	
	my $return_value=0;
	if($type eq 'first_line'){
		open my ($IN), "<", $file or die;
		$return_value=<$IN>;
		close($IN);
	}
	elsif($type eq 'lines_number'){
		open my($COM), "wc -l $file |" or die;
		while(<$COM>){
			chomp($_);
			my @a=split(" ", $_);
			$return_value=$a[0];
		}
		close($COM);
	}
	elsif($type eq 'bytes'){
		my @a=stat($file);
		$return_value=$a[7];
	}
	elsif($type eq 'KB'){
		my @a=stat($file);
		$return_value=int(($a[7]/1024)+0.5);
	}
	elsif($type eq 'MB'){
		my @a=stat($file);
		$return_value=int(($a[7]/1024/1024)+0.5);
	}
	return($return_value);
}

###########################################
#judge if a string exists in  one file
sub line_judging{
	my ($file, $str)=@_;
	
	my $judging='F';
	if (-f $file){
		open my ($IN), "<", $file or die;
		while(<$IN>){
			if ($_=~/$str/){
				$judging='T';
				last;
			}
		}
		close($IN);
	}
	return($judging);
}

#########################################################
#hash into two columns txt file
sub hash_to_txt{
  my $hash_pointer=$_[0];
  my %hash=%$hash_pointer;
  my $txt_file=$_[1];
  
  open my ($OUT), ">", $txt_file or die; 
  foreach my $key (sort(keys %hash)){
    print $OUT "$key=$hash{$key}\n";
  }
  close($OUT);
}


#########################################################
# two columns txt file into hash
sub txt_to_hash{
  my $txt_file=$_[0];
  my %hash;
  open my ($IN), "<", $txt_file or die; 
  while(my $line=<$IN>){
    chomp($line);
    my($key, $value)=split("=", $line);
    $hash{$key}=$value;
  }
  close($IN);
  return(\%hash);
}
####################################################
#
sub scan_hash{
  my $hash_pointer=$_[0];
  my %hash=%$hash_pointer;
  my $redundant_num=0;
  while (my($key, $value)=each(%hash)){
    my @items=split("_", $key);
    $redundant_num +=$items[-1];
    #print "$key\t$value\n";
  }
  return($redundant_num);
}

#######################################
#inspect bowtie index or build it
sub check_bowtie_index{
	my($bowtie_dir, $version, $fasta_file, $index_name)=@_;
	my @bowtie1_tail=('.1.ebwt', '.2.ebwt','.3.ebwt','.4.ebwt','.rev.1.ebwt','.rev.2.ebwt');
	my @bowtie2_tail=('.1.bt2','.2.bt2','.3.bt2','.4.bt2','.rev.1.bt2','.rev.2.bt2',);
	
	my $judging=1;
	if ($version=~/bowtie1/gi){
		foreach(@bowtie1_tail){
			my $index_file=$index_name.$_;
			$judging=0 unless -f $index_file;
		}
		system("$bowtie_dir/bowtie-build $fasta_file $index_name") if $judging==0;
	}
	elsif ($version=~/bowtie2/gi){
		foreach(@bowtie2_tail){
			my $index_file=$index_name.$_;
			$judging=0 unless -f $index_file;
		}
		system("$bowtie_dir/bowtie2-build $fasta_file $index_name") if $judging==0;
	}
	return($judging);
}
###################################################################
#function
sub pre_get_sample_file{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	
	#associate sample name with *.fastq files list
	my @fastq_files;
	my @rawdata_dirs=split(',', $variables{dir_rawdata});
	foreach my $rawdata_dir(@rawdata_dirs){
		my $files_pointer=E_RNA::files_list($rawdata_dir, 'incrusive_files', "fastq,fq");
		my @files=@$files_pointer;
		foreach(@files){
			push(@fastq_files, $_);
		}
	}
	#
	my %sample_info;
	my @R1_fastq_files=grep(/_R1/, @fastq_files);
	foreach my $R1_file(@R1_fastq_files){
		my @a1=split("/", $R1_file);
		my $R1_file_name=$a1[-1];
		my $R1_file_name_head=$a1[-1];
		$R1_file_name_head=~s/\.fastq|\.fq//;
		my @a2=split("_", $R1_file_name_head);
		my $sample_name=$a2[0];
		$sample_name .= '_'.$a2[1] if $a2[1];
		$sample_info{$sample_name}=$sample_name;
	}
	#export
	open my($OUT), ">", $variables{file_sample_info}, or die;
	print $OUT join(";", 'sample_name', 'fastq_name'), "\n";
	foreach my $name(sort (keys %sample_info)){
		print $OUT join(";", $name, $sample_info{$name}), "\n";
	}
	close($OUT);
	
}

########################
#ge the hash %sample_info
sub pre_sample_info{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	
  	#read sample_info file
	my (%sample_info, @col_names);
	my $row=0; 
	open my ($IN), "<", $variables{file_sample_info} or die;
	while(<$IN>){#2
		chomp($_);
		my @array=split(";", $_);
		$row++;
		if($row==1){#3
			@col_names=@array;
		}#3
		else{#3
			for(my $col=1; $col<@array; $col++){#4
				my $sample_name=$array[0];
				my $col_name=$col_names[$col];
				$col_name='fastq_names' if $col==1; #the second column
				$sample_info{$sample_name}->{$col_name}=$array[$col];
			}#4
		}#3
	}#2
	close($IN);
	#get @sample_names ordered by file size
	my @sample_names=keys %sample_info;
	@sample_names=sort @sample_names;
	
	print "read raw data from $variables{dir_rawdata}\n";
	my %total_fastq;
	my @rawdata_dirs=split(',', $variables{dir_rawdata});
	foreach my $rawdata_dir(@rawdata_dirs){#2
		my $files_pointer=E_RNA::files_list($rawdata_dir, 'incrusive_files', "fastq,fq");
		my @files=@$files_pointer;
		foreach my $file(@files){
			my @a=split('/', $file);
			my $file_name=$a[-1];
			$total_fastq{$file}=$file_name;
		}
	}#2
	
	#
	foreach my $sample_name(keys %sample_info){#2
		#assign fastq file to $sample
		my @fastq_files;
		my @fastq_names=split(',', $sample_info{$sample_name}->{'fastq_names'});
		foreach my $name(@fastq_names){#3
			foreach my $file(keys %total_fastq){
				my $file_name=$total_fastq{$file};
				#$name must be the head of fastq files seperated by '_'
				if ($file_name=~/^$name/){
					push(@fastq_files, $file);
					delete $total_fastq{$file};
				} 
			}
		}#3
		#add sample information into %sample_info
		my @fastq_R1_files=grep(/_R1/, @fastq_files);
		my @fastq_R2_files=grep(/_R2/, @fastq_files);
		$sample_info{$sample_name}->{'fastq_files'}=(@fastq_files>0) ? join(',', @fastq_files) : 'NA';
		$sample_info{$sample_name}->{'fastq_R1_files'}=(@fastq_R1_files>0) ? join(',', @fastq_R1_files) : 'NA';
		$sample_info{$sample_name}->{'fastq_R2_files'}=(@fastq_R2_files>0) ? join(',', @fastq_R2_files) : 'NA';
		#return the size of all fastq files
		$sample_info{$sample_name}->{'files_size'}=0;
		foreach(@fastq_files){
			my @args=stat($_);
			$sample_info{$sample_name}->{'files_size'} +=$args[7]; #byte
		}
		print "$sample_name:\n", "\t$sample_info{$sample_name}->{'fastq_files'}\n";
	}#2
	
	#assign the directory per $sample_name;
	my @dir_array;
	open my($OUT), ">", $variables{file_storage_log} or die;
	foreach my $sample_name(@sample_names){ #2
		@dir_array=split(',', $variables{dir_result_array}) if @dir_array==0;
		my $out_dir=shift @dir_array;
		$out_dir .= '/' unless $out_dir=~/\/$/;
		my $sample_dir=$out_dir.$sample_name;
		$sample_info{$sample_name}->{sample_dir}=$sample_dir;
		print $OUT "$sample_name=$sample_dir\n";
	}#2
	close($OUT);
	
	return(\%sample_info);
}
########################################################################
#function: the specific files list in the specific direcotry


#list_type: files, incrusive_files, incrusive_file_names, file_names, and sample_names
sub files_list{#1
	my ($dir, $list_type, $file_tail)=@_;
	$dir .= '/' unless $dir=~/\/$/;
	my @file_tails=split(',', $file_tail) if $file_tail;
	
	#read information under directory
	opendir(DIR, $dir);
	my @all=readdir(DIR);
	closedir(DIR);
	my (@file_names, @dir_names);
	foreach (@all){
		unless ($_ =~ /^\.|~/){
			if (-d $dir.$_) {	push(@dir_names, $_) ;		}
			else {	push(@file_names, $_) ;		}
		}
	}
	#print @file_names, "\n";
	
	my @files;
	if ($list_type eq 'incrusive_directories'){
		File::Find::find(sub {
			my $file=$File::Find::name;
			push(@files, $file) if -d $file;
		}, $dir);
		@files=grep(/$file_tail/, @files) if $file_tail;
	}
	elsif ($list_type eq 'directories'){
		@files=@dir_names;
		@files=grep(/$file_tail/, @files) if $file_tail;
	}
	elsif ($list_type eq 'incrusive_files'){
		File::Find::find(sub {
			my $file=$File::Find::name;
			push(@files, $file) if (-f $file) and List::Util::first {$file=~/\.$_$/} @file_tails;
		}, $dir);
	}
	elsif ($list_type eq 'files'){
		if ($_[2]){
			my @sub_file_names=grep(/\.$file_tail$/, @file_names) ;
			@files=map {$dir.$_} @sub_file_names;
		}
		else{	@files=map {$dir.$_} @file_names;		}
	}
	elsif ($list_type eq 'incrusive_file_names'){
		File::Find::find(sub {
			my $file=$File::Find::name;
			my @array=split("/", $file);
			my $name=$array[-1];
			push(@files, $name) if (-f $file) and $name=~/\.$file_tail$/;
		}, $dir);
		@files=List::MoreUtils::uniq(@files);
		@files=sort @files;
	}
	elsif ($list_type eq 'file_names'){
		@files=$file_tail ? grep(/\.$file_tail$/, @file_names) : @file_names;
		@files=List::MoreUtils::uniq(@files);
		@files=sort @files;
	}
	elsif($list_type eq "sample_names"){ #get all of sample name in the dir based on fastq files
		my %sample_size;
		File::Find::find(sub {
				my $file=$File::Find::name;
				if (-f $file and $file=~/\.fastq$|\.fq$/){#4
					my @array1=split("/", $file);
					my $file_name=$array1[-1];
					my @array2=split("_", $file_name);
					my $sample_name=$array2[0];
					$sample_name .= "_".$array2[1] if $array2[1];
					$sample_name .= "_".$array2[2] if $array2[2];
					my @args=stat($file);
					$sample_size{$sample_name} += $args[7];
				}#4
			}, $dir);
		foreach my $sample_name(sort { $sample_size{$b}<=>$sample_size{$a} } (keys %sample_size) ){
			push(@files, $sample_name);
		}
	}

	#print "@files\n";
	return(\@files);
}#1

##############################################
sub calculate_dir_room{
	my ($dir)=@_;
	
	#calculate the left room of @dir_array 
	my %dir_left;
	my @dir_array=split(',', $dir);
	foreach my $dir(@dir_array){
		#get the unused room(GB) of $dir;
		my @room_info=readpipe ("df -hl $dir") or die 'wrong df command!\n';
		my @items=split(" ", $room_info[1]);
		my $left_room=$items[3];
		if($left_room=~s/G$//){	$left_room=$left_room;				}
		elsif($left_room=~s/T$//){	$left_room=int($left_room*1024);		}
		else{	$left_room=1;		}
		$dir_left{$dir}=$left_room;
	}
	
	return(\%dir_left);
}
############################################################
sub pre_references_info{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my %ref_info;
	
	#read old log file
	if( -f $variables{file_ref_log}){ 
		open my($IN), "<", $variables{file_ref_log}or die;
		while(<$IN>){
			chomp($_);
			my($index_name, $seq_num, $total_base, $ave_len)=split("\t", $_);
			$ref_info{$index_name}->{seq_num}=$seq_num;
			$ref_info{$index_name}->{total_base}=$total_base;
			$ref_info{$index_name}->{ave_len}=$ave_len;
		}
		close($IN);
	}
	
	#generate the new one
	my $files_pointer=E_RNA::files_list($variables{dir_ref_seq}, 'file_names', 'fa');
	my @fasta_names=@$files_pointer;

	open my($OUT), ">", $variables{file_ref_log} or die;
	print $OUT join("\t", 'index_name', 'seq_num', 'total_base', 'ave_len'), "\n";
	foreach my $fasta_name(@fasta_names){#2
		print "\t$fasta_name\n";
		my $index_name=$fasta_name;
		$index_name=~s/\.fasta$|\.fa$//;
		unless (exists $ref_info{$index_name}){
			my $seq_num=0;
			my $total_base=0;
			my $in_obj = Bio::SeqIO->new(-file => $variables{dir_ref_seq}.'/'.$fasta_name, -format => 'fasta');
			while (my $seq_obj = $in_obj->next_seq() ) {#3
				my $displayid=$seq_obj->display_id();
				my $seq_len=$seq_obj->length();
				$seq_num++;
				$total_base += $seq_len;
			}#3
			$ref_info{$index_name}->{seq_num}=$seq_num;
			$ref_info{$index_name}->{total_base}=$total_base;
			$ref_info{$index_name}->{ave_len}=int($total_base/$seq_num+0.5);
		}
		print $OUT join("\t", $index_name, $ref_info{$index_name}->{seq_num}, 
						$ref_info{$index_name}->{total_base}, $ref_info{$index_name}->{ave_len}), "\n";
	}#2
	close($OUT);

	$variables{ref_info_pointer}=\%ref_info;
	return(\%variables);
}

############################################################
sub pre_check_bowtie_index{
	my ($variables_pointer, $version)=@_;
	my %variables=%$variables_pointer;
	
	#check index files
	my $found=1;
	if($version eq 'bowtie2'){#2
		my @tails=('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2');
		foreach (@tails){
			$found=0 unless -f $variables{genome_index}.$_;
		}
		#build index 
		if($found==0){
			my $bowtie_build_script=$variables{dir_bowtie2}.'/bowtie2-build';
			system("$bowtie_build_script $variables{genome_fasta_file} $variables{genome_index}");
		}
	}#2
	elsif($version eq 'bowtie1'){#2
		my @index_names;
		my @a=split(",", $variables{index_seperate}.','.$variables{index_iterative});
		foreach my $name(@a){
			unless ($name eq 'no'){
				push(@index_names, $name) unless List::Util::first {$_ eq $name} @index_names;
			}
		}
		
		foreach my $index_name(@index_names){#3
			my $bowtie_index=$variables{dir_bowtie1}.'/'.$index_name;
			my $fasta_file=$variables{dir_ref_seq}.'/'.$index_name.'.fa';
			#check index
			my @tails=('.1.ebwt', '.2.ebwt', '.3.ebwt', '.4.ebwt', '.rev.1.ebwt', '.rev.2.ebwt');
			foreach (@tails){
				$found=0 unless -f $bowtie_index.$_;
			}
			#build index 
			if($found==0){
				my $bowtie_build_script=$variables{dir_bowtie1}.'/bowtie-build';
				system("$bowtie_build_script $fasta_file $bowtie_index");
			}
		}#3
	}#2
	#
}
#################################################
#
sub mRNA_tophat_mapping{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my $sample_name=$variables{sample_name};
	my $sample_dir=$variables{sample_dir};
	my $sample_log=$variables{sample_log};
	my @fastq_R1_files=split(',', $sample_info{$sample_name}->{fastq_R1_files});
	my @fastq_R2_files=split(',', $sample_info{$sample_name}->{fastq_R2_files});
	my @lines_arr=map {E_RNA::file_info($_, 'lines_number')} @fastq_R1_files;
	my $raw_pairs=(List::Util::sum @lines_arr)/4;
	E_RNA::refresh_log($sample_log, "raw_pairs", $raw_pairs);
	
	my $tophat_file=$variables{dir_mapper}.'/tophat2';
	my $tophat_ouput_dir=$sample_dir.'/tophat';
	my $tophat_annotation_dir=$tophat_ouput_dir.'/genome_annotation';
	mkdir($tophat_ouput_dir, 0755) unless -d $tophat_ouput_dir;
	mkdir ($tophat_annotation_dir, 0755) unless -d $tophat_annotation_dir;
	
	#copy fa into annotation dir
	my @a=split('/', $variables{genome_fasta_file});
	my $tophat_fasta_file=$tophat_annotation_dir.'/'.$a[-1];
	system("cp $variables{genome_fasta_file} $tophat_fasta_file");
	
	#copy gtf/gff into annotation dir
	my $tophat_annotation_file=$tophat_annotation_dir.'/'.$variables{genome_index_name}.'.gff';
	system("cp $variables{genome_gtf_file} $tophat_annotation_file");
	
	#copy bowite index into annotation dir
	my $found=1;
	my @tails=('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2');
	foreach (@tails){
		my $source_file=$variables{genome_index}.$_;
		my $dest_file=$tophat_annotation_dir.'/'.$variables{genome_index_name}.$_;
		system("cp $source_file $dest_file");
	}

	print "\n\n Genome mapping of $sample_name:\n";
	my $tophat_options=join(' ', $variables{tophat_options}, 
						'--transcriptome-index', $tophat_annotation_dir.'/'.$variables{genome_index_name},
						'-o', $tophat_ouput_dir, $variables{genome_index}		);
	print "mapping R1_fastq and R2_fastq files of $sample_name to the reference genome using TopHat!\n\n";
	print "$tophat_file $tophat_options $sample_info{$sample_name}->{fastq_R1_files} $sample_info{$sample_name}->{fastq_R2_files}\n";
	system("$tophat_file $tophat_options $sample_info{$sample_name}->{fastq_R1_files} $sample_info{$sample_name}->{fastq_R2_files}");
	print "Alignment of $sample_name is done!\n";
	#remove annotation dir
	system("rm -R $tophat_annotation_dir") if -d $tophat_annotation_dir;
	
}
##############################################
sub mRNA_sam_analysis{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my $sample_name=$variables{sample_name};
	my $sample_dir=$variables{sample_dir};
	my $sample_log=$variables{sample_log};
	
	#print "Convert *.bam files into *.sam files\n";
	my $bam_file=$sample_dir.'/tophat/accepted_hits.bam';
	my $sam_file=$sample_dir.'/tophat/accepted_hits.sam';
	#system("samtools view $bam_file > $sam_file");
	$bam_file=$sample_dir.'/tophat/unmapped.bam';
	$sam_file=$sample_dir.'/tophat/unmapped.sam';
	#system("samtools view $bam_file > $sam_file");
	
	#
	
}
#################################################
#
sub mRNA_cufflinks_assembling{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my $sample_name=$variables{sample_name};
	my $sample_dir=$variables{sample_dir};
	my $sample_log=$variables{sample_log};
	
	print "\n\n Transcripts assembling of $sample_name:\n";
	#assemble expressed genes and transcripts using cufflinks
	my $cufflinks_file=$variables{dir_assembler}.'/cufflinks';
	my $cufflinks_output_dir=$sample_dir.'/cufflinks';
	my $cufflinks_options=join(' ', $variables{cufflinks_options}, '-o', $cufflinks_output_dir);
	my $bam_file=$sample_dir.'/tophat/accepted_hits.bam';
	print "###$cufflinks_file $cufflinks_options $bam_file\n";
	system("$cufflinks_file $cufflinks_options $bam_file");
		
	#compare assembled transcripts to a reference annotation
	my $cuffcompare_file=$variables{dir_assembler}.'/cuffcompare';
	my $cuff_gtf_file=$sample_dir.'/cufflinks/transcripts.gtf';
	my $output_prefix=$sample_dir.'/cufflinks/cuffcomp';
	my $cuffcompare_options=join(" ", $cuff_gtf_file, '-o', $output_prefix);
	print "###$cuffcompare_file $cuffcompare_options\n";
	system("$cuffcompare_file $cuffcompare_options");
	print "Transcripts assembling of $sample_name is done!\n";
	
	
}
#############################################
##########################################
#
sub tophat_cufflinks_FPKM_export {
	my ($variables_pointer, $type)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my @sample_names=split(",", $variables{sample_names});
	
	#
	my %FPKM;
	foreach my $sample_name(@sample_names){
		my $sample_dir=$sample_info{$sample_name}->{sample_dir};
		my $n=1;
		open my($IN), "<", $sample_dir.'/cufflinks/'.$type.'.fpkm_tracking';
		print "$sample_name:", $sample_dir.'/cufflinks/'.$type.'.fpkm_tracking', "\n";
		while(<$IN>){
			if ($n>1){
				chomp($_);
				my @items=split("\t", $_);
				my $tracking_id=$items[0];
				#my $gene_id=$items[3];
				my $fpkm=$items[9];
				$FPKM{$tracking_id}->{$sample_name}=$fpkm;
			}
			$n++;
		}
		close($IN);
	}
	my @id=keys %FPKM;
	@id=sort @id;
	
	#export
	my $fpkm_csv=$variables{dir_result}.'/'.$type.'_fpkm.csv';
	open my($OUT), ">", $fpkm_csv or die;
	print $OUT join(',', 'tracking_id', @sample_names), "\n";
	foreach my $tracking_id(@id){
		my @fpkm_values;
		foreach my $sample_name(@sample_names){
			my $value=(exists $FPKM{$tracking_id}->{$sample_name}) ? $FPKM{$tracking_id}->{$sample_name} : 0;
			push(@fpkm_values, $value);
		}
		print $OUT join(',', $tracking_id, @fpkm_values), "\n";
	}
	close($OUT);
}

############################################
#open parameter info file named "E_RNA_pipeline.info" needed for processing. 
sub pre_process_info{#1
	my %variables;  # save parameters for processing:
	open my ($INFO), "<", $_[0] or die "could not open the file of variables for E_RNA pipeline!";
	while (<$INFO>) {#2
		chomp($_);
		if($_=~/=/ and $_!~/#/){
			my ($name, $value) = split("=", $_); #split on the tabs
			$variables{$name}=$value;
		}
	}#2
	close($INFO);

	return(\%variables);
}#1

#############################################
sub refresh_R_script{
	my ($log_file, $r_name, $r_value)=@_;
	
	my @lines;
	#read old data
	open my($IN), "<", $log_file or die;
	while (<$IN>){
		chomp($_);
		if($_=~/^R\_/){
			my($name, $value)=split("=", $_);
			if($name eq $r_name){
				my $r_line=$name."=\"".$r_value."\"";
				push(@lines, $r_line);
			}
			else{	push(@lines, $_);	}
		}
		else{
			push(@lines, $_);
		}
	}
	close($IN);
	
	#export refreshed data
	open my($OUT), ">", $log_file or die;
	foreach (@lines){
		print $OUT "$_\n";
	}
	close($OUT);
}

#############################################
#type= 'replace' or 'add'
sub refresh_log{
	my ($log_file, $r_name, $r_value, $type)=@_;
	
	my %variables;
	#print "###$log_file, $r_name, $r_value\n";
	#read old data
	if(-f $log_file){
		open my($IN), "<", $log_file or die;
		while (<$IN>){
			chomp($_);
			if($_=~/=/){
				my($name, $value)=split("=", $_);
				$variables{$name}=$value;
			}
		}
		close($IN);
	}
	#refresh new data
	if (exists $variables{$r_name} and $type and $type eq 'add'){
		$variables{$r_name}=$variables{$r_name}+$r_value;
	}
	else{
		$variables{$r_name}=$r_value;
	}
	#export refreshed data
	open my($OUT), ">", $log_file or die;
	foreach my $key( sort {$a cmp $b} (keys %variables ) ){
		print $OUT "$key=$variables{$key}\n";
	}
	close($OUT);

}

#############################################
sub read_log{
	my ($log_file, $variable_name)=@_;
	
	my $variable_value='NA';
	open my($IN), "<", $log_file or die;
	while (<$IN>){
		if($_=~/=/){
			chomp($_);
			my($name, $value)=split("=", $_);
			$value=~s/\"//g;
			$variable_value=$value if $name eq $variable_name;
		}
	}
	close($IN);
	
	return($variable_value);
}
#######################################################################
#generate specie miRNA fasta file based on miRBase fasta files
sub specie_miRNA{
	my($precursor_fasta_file, $mature_fasta_file, $specie, $dir)=@_;
	$dir .= '/' unless $dir=~/\/$/;
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
#####################################################
#read Q values 
sub pre_sequencing_quality{ #1
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my $sample_name=$variables{sample_name};
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	
	my %hash=('fastq_R1_files'=>'.R1_QC', 'fastq_R2_files'=>'.R2_QC');
	foreach my $item(keys %hash){#2
		my $file_tail=$hash{$item};
		unless($sample_info{$sample_name}->{$item} eq 'NA'){#3
			my @fastq_files=split(',', $sample_info{$sample_name}->{$item});
			open my($QC), ">", $variables{sample_dir}.'/'.$sample_name.$file_tail or die;
			foreach my $fastq_file(@fastq_files){#4
				my @quals;
				my $num=0;
				print "\tCheck sequencing quality of $fastq_file.\n";
				my $in_obj = Bio::SeqIO->new(-file => $fastq_file, -format => 'fastq');
				while (my $seq = $in_obj->next_seq() ) {#5
					my $quality_pointer=$seq->Bio::Seq::Quality::qual;
					my @quality=@$quality_pointer;
					@quals= ($num==0) ? @quality : List::MoreUtils::pairwise {$a+$b} @quals, @quality;
					$num++;
					if ($num==$variables{QC_compression}){
						@quals=map{int($_/$num+0.5)} @quals;
						print $QC join("\t", @quals), "\n";
						$num=0;
						undef @quals;
					}
				}#5
				if($num>0){
					@quals=map{int($_/$num+0.5)} @quals;
					print $QC join("\t", @quals), "\n";
				}
			}#4
			close($QC);
		}#3
	}#2

	#return (\%variables);
}#1

#####################################################
#read all fastq files of one sample 
#generate *.RC files in the result file
sub ncRNA_read_rawdata{ #1
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my $sample_name=$variables{sample_name};
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	
	my $raw_reads_num=0;
	my @R1_fastq_files=split(',', $sample_info{$sample_name}->{fastq_R1_files});
	open my($OUT), ">", $variables{sample_dir}.'/'.$sample_name.".RC" or die; 
	foreach my $fastq_file(@R1_fastq_files){#2 fastq file circyling
		my %pairs; # grabs the FASTQ parser, specifies the Illumina variant
		
		#read R1 raw data
		print "\t read $fastq_file!\n";
		my $in_obj = Bio::SeqIO->new(-file => $fastq_file, -format => 'fastq');
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
		$fastq_file=~s/_R1_/_R2_/;
		if (-f $fastq_file){#3
			print "\t read $fastq_file!\n";
			my $in_obj = Bio::SeqIO->new(-file => $fastq_file, -format => 'fastq');
			while (my $seq_obj = $in_obj->next_seq() ) {#4
				my @display_id=split(/:/, $seq_obj->display_id());
				my $coordinate=$display_id[-2]."_".$display_id[-1];
				$pairs{$coordinate}->{R2_seq}=$seq_obj->seq();
				$raw_reads_num++;
				$pairs{$coordinate}->{R2_No}=$raw_reads_num;
			}#4
		}#3

		#export pair ends of read sequences of R1 and R2
		foreach (keys %pairs){#3
			print $OUT join("\t", $_, $pairs{$_}->{R1_No}, $pairs{$_}->{R1_seq},
						$pairs{$_}->{R2_No}, $pairs{$_}->{R2_seq}, $fastq_file), "\n";
		}#3
	}#2 fastq file circyling
	close($OUT);
	E_RNA::refresh_log($variables{sample_log}, "raw_reads_num", $raw_reads_num);

	#return (\%variables);
}#1

##################################################################
#function
sub reverse_complement {
	my $dna_seq = shift;

	# reverse the DNA sequence
	my $revcom = reverse($dna_seq);
	# complement the reversed DNA sequence
	$revcom =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;

	return ($revcom);
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


##########################################################
#adapter removal at the 3' end  
sub ncRNA_adapter_removal{#1
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my $sample_name=$variables{sample_name};
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	
	#generate the possible adapters for seq truncation
	my $adapters_pointer=adapters(\%variables);
	my %adapters=%$adapters_pointer;

	my (%insert, %insert_len);
	open my($RC), "<", $variables{sample_dir}.'/'.$sample_name.".RC" or die;
	open my($INSERT), ">", $variables{sample_dir}.'/'.$sample_name.".insert_RC" or die;
	while (<$RC>){#2
		chomp($_);
		my ($coordinate, $R1_seq_No, $R1_read_seq, $R2_seq_No, $R2_read_seq )=split("\t", $_);
		#truncate R1_seq
		my $R1_seq_trunc=adapter_3_truncation(\%variables, $R1_read_seq, \%adapters);
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
			my $R2_seq_trunc=adapter_5_truncation(\%variables, $R2_read_seq, \%adapters);
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
	open my($query), ">", $variables{sample_dir}.'/'.$sample_name."_query.fa" or die;
	open my($TXT), ">", $variables{sample_dir}.'/'.$sample_name.".query_RC" or die;
	delete $insert{'NA'};
	foreach my $seq(keys %insert){#2
		if (length($seq)>=$variables{query_length}){#4
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
	E_RNA::refresh_log($variables{sample_log}, "insert_reads_num", $insert_reads_num);
	E_RNA::refresh_log($variables{sample_log}, "nr_insert_reads_num", $nr_insert_reads_num);
	E_RNA::refresh_log($variables{sample_log}, "query_reads_num", $query_reads_num);
	E_RNA::refresh_log($variables{sample_log}, "nr_query_reads_num", $nr_query_reads_num);
	
	#export insert length distribution
	my $IL_file=$variables{sample_dir}.'/'.$sample_name.'.IL';
	open my($IL), ">", $IL_file or die;
	foreach my $len(sort {$a <=> $b} (keys %insert_len)){
		print $IL "$len\t", "$insert_len{$len}\n";
	}
	close($IL);

	return(\%variables);
}#1

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


####################################################
#read SAM file determined by bowtie1 and bowtie2
#1S12_GGCTAC:4316023_1   0       precursor_miRNA:hsa-mir-451a    17      255     18M     *       0       0       AGACCGTTACCATTACTG      IIIIIIIIIIIIIIIIII      XA:i:1  MD:Z:1A16       NM:i:1
sub read_bowtie1_sam{
	my ($version, $fasta_file, $alignment_output)=@_;

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
			my $flags=$items[1];
			my $ref_name=$items[2];
			my $ref_offset=$items[3];
			my $read_seq=$items[9];
			my ($sample_name, $query_str)=split(':', $query_name);
			my ($seq_No, $query_read_counts)=split('_', $query_str);
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
					$mappable_counts{$ref_name}->{seq_info} .= ';'.$seq_No.','.$read_seq;
					$mappable_counts{$ref_name}->{query_name_str} .= ','.$query_name;
				}
				else{	
					$mappable_counts{$ref_name}->{seq_info} =$seq_No.','.$read_seq;
					$mappable_counts{$ref_name}->{query_name_str} =$query_name;
				}
				#update %SD_counts
				$SD_counts{$seq_No}=$ref_name;
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
###################################################
#sequence alignment or genome mapping using bowtie 
#the software named bowtie-build and bowtie should be copied into the bowtie directory
sub ncRNA_seperate_alignment{#1
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my $sample_name=$variables{sample_name};
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	
	#sequence alignment using bowtie
	my @index_names=split(",", $variables{index_seperate}) ;
	foreach my $index_name(@index_names){#2  #index name circyling
		#get all read files with the name tail *_combined.fasta 
		my $bowtie_index=$variables{dir_bowtie1}.'/'.$index_name;
		my $fasta_file=$variables{sample_dir}.'/'.$sample_name."_query.fa";
		my $alignment_sam_output=$variables{sample_dir}.'/'.$sample_name."_".$index_name.'.sam';
		print "Seperate alignment: $variables{bowtie_options} $bowtie_index -f $fasta_file -S $alignment_sam_output\n\n";
		system("$variables{bowtie_options} $bowtie_index -f $fasta_file -S $alignment_sam_output ");
		my $alignment_out_pointer=E_RNA::read_bowtie1_sam('bowtie1', $fasta_file, $alignment_sam_output);
		my %alignment_out=%$alignment_out_pointer;
		
		#my $bam_file=$variables{sample_dir}.'/'.$sample_name."_".$index_name.'.bam';
		#my $sorted_bam_file=$variables{sample_dir}.'/'.$sample_name."_".$index_name.'.sorted';
		#print "generate bam file: $bam_file\n";
		#system("samtools view -bS $alignment_sam_output -o $bam_file");
		#system("samtools sort $bam_file $sorted_bam_file");
		#system("samtools index $sorted_bam_file.bam");
		
		#export miRNA quanlification (read counts of known miRNA versus mappable query read counts) 
		my $mappable_counts_pointer=$alignment_out{mappable_counts_pointer};
		my %mappable_counts=%$mappable_counts_pointer;
		open my($COU), ">", $variables{sample_dir}.'/'.$sample_name."_".$index_name.".alignment" or die; 
		foreach my $ref_name( sort (keys %mappable_counts) ){#3
			print $COU join("\t", $sample_name, $index_name, $ref_name, 
				$mappable_counts{$ref_name}->{top_RC}, $mappable_counts{$ref_name}->{middle_RC},
				$mappable_counts{$ref_name}->{bottom_RC}, $mappable_counts{$ref_name}->{seq_info}, ), "\n";
		}#3
		close($COU);
		
		#update sample log file
		my $num_out_pointer=$alignment_out{num_out_pointer};
		my %num_out=%$num_out_pointer;
		while( my($num_name, $value)=(each %num_out) ){
			E_RNA::refresh_log($variables{sample_log}, $index_name.':'.$num_name, $value);
		}
		print  "mappable redundant number: $num_out{nr_query_reads_num}\t";
		print  "$num_out{nr_mappable_reads_num}\t", "$num_out{nr_unmappable_reads_num}\n";
		
		#export unalignment and multiple alignment sequences
		my $un_alignment_pointer=$alignment_out{un_alignment_pointer};
		my %un_alignment=%$un_alignment_pointer;
		open my($UN), ">", $variables{sample_dir}.'/'.$sample_name."_".$index_name."_unalignment.fa" or die; 
		foreach my $query_name(sort {$un_alignment{$b}->{RC}<=>$un_alignment{$a}->{RC}} (keys %un_alignment)){
			print $UN ">$query_name\n", "$un_alignment{$query_name}->{seq}\n";
		}
		close($UN);
		my $multiple_alignment_pointer=$alignment_out{multiple_alignment_pointer};
		my %multiple_alignment=%$multiple_alignment_pointer;
		open my($MUL), ">", $variables{sample_dir}.'/'.$sample_name."_".$index_name."_multialignment.fa" or die; 
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
		open my($SD), ">", $variables{sample_dir}.'/'.$sample_name."_".$index_name.".SD" or die; 
		foreach my $seq_No(sort {$a<=>$b} keys(%SD_counts) ){
			my $ref_name=$SD_counts{$seq_No};
			push(@ref_names, $ref_name) unless List::Util::first {$ref_name eq $_} @ref_names;
			$ref_num=@ref_names;
			print $SD "$seq_No\t", "$ref_num\n";
		}
		my $seq_No_max=List::Util::max keys %SD_counts;
		my $raw_reads_num=E_RNA::read_log($variables{sample_log}, 'raw_reads_num');
		print $SD "$raw_reads_num\t", "$ref_num\n" if $raw_reads_num > $seq_No_max;
		close($SD);
	}#2 index name circyling end
  
	#return(\%variables);
}#1


##############################################################
#iterativesequence alignment using bowtie
sub ncRNA_iterative_alignment{#1
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
  	my $fasta_file=$variables{sample_dir}.'/'.$sample_name."_query.fa";
  	my $fasta_tmp=$fasta_file.'_tmp';
	system("cp $fasta_file $fasta_tmp");
	my $sam_tmp=$variables{sample_dir}.'/'.$sample_name.'_iterative.sam_tmp';
	
	#iterative alignment
	my (%total_num_out, %SD_counts); #total number counting
	my @index_names=split(",", $variables{index_iterative});
	open my($SAM), ">", $variables{sample_dir}.'/'.$sample_name.'_iterative_combined.sam' or die;
	open my($COU), ">", $variables{sample_dir}.'/'.$sample_name.'_iterative.alignment' or die;
	for (my $i=0; $i<@index_names; $i++){#2  index_name circyling
		#bowtie alignment
		my $index_name=$index_names[$i];
		my $bowtie_index=$variables{dir_bowtie1}.'/'.$index_name;
		print  "Iterative alignment of $sample_name against $bowtie_index!\n";
		system("$variables{bowtie_options} $bowtie_index -f $fasta_tmp -S $sam_tmp ");
		my $alignment_out_pointer=E_RNA::read_bowtie1_sam('bowtie1', $fasta_tmp, $sam_tmp);
		my %alignment_out=%$alignment_out_pointer;
		
		#export miRNA quanlification (read counts of known miRNA versus mappable query read counts) 
		my $mappable_counts_pointer=$alignment_out{mappable_counts_pointer};
		my %mappable_counts=%$mappable_counts_pointer;
		foreach my $ref_name( sort (keys %mappable_counts) ){#3
			print $COU join("\t", $sample_name, 'iterative_'.$i.'_'.$index_name, $ref_name, 
				$mappable_counts{$ref_name}->{top_RC}, $mappable_counts{$ref_name}->{middle_RC},
				$mappable_counts{$ref_name}->{bottom_RC}, $mappable_counts{$ref_name}->{seq_info}, ), "\n";
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
		open my($UN), ">", $variables{sample_dir}.'/'.$sample_name.'_iterative_'.$i.'_'.$index_name."_unalignment.fa" or die; 
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
		open my($MUL), ">", $variables{sample_dir}.'/'.$sample_name.'_iterative_'.$i.'_'.$index_name."_multialignment.fa" or die; 
		foreach my $query_name(sort {$multiple_alignment{$b}->{RC}<=>$multiple_alignment{$a}->{RC}} (keys %multiple_alignment) ){
			print $MUL ">$query_name\n", "$multiple_alignment{$query_name}->{seq}\n";
		}
		close($MUL);
		
		#update %SD_counts
		my $SD_counts_pointer=$alignment_out{SD_counts_pointer};
		my %tmp_SD=%$SD_counts_pointer;
		while (my($tmp_seq_No, $tmp_ref_name)=(each %tmp_SD)){
			$SD_counts{$tmp_seq_No}=$tmp_ref_name;
		}
		#update sample log file
		my $num_out_pointer=$alignment_out{num_out_pointer};
		my %num_out=%$num_out_pointer;
		while( my($num_name, $value)=(each %num_out) ){
			E_RNA::refresh_log($variables{sample_log}, 'iterative_'.$i.'_'.$index_name.':'.$num_name, $value);
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
	open my($SD), ">", $variables{sample_dir}.'/'.$sample_name."_iterative.SD" or die; 
	foreach my $seq_No(sort {$a<=>$b} keys(%SD_counts) ){
		my $ref_name=$SD_counts{$seq_No};
		push(@ref_names, $ref_name) unless List::Util::first {$ref_name eq $_} @ref_names;
		$ref_num=@ref_names;
		print $SD "$seq_No\t", "$ref_num\n";
	}
	my $seq_No_max=List::Util::max keys %SD_counts;
	my $raw_reads_num=E_RNA::read_log($variables{sample_log}, 'raw_reads_num');
	print $SD "$raw_reads_num\t", "$ref_num\n" if $raw_reads_num > $seq_No_max;
	close($SD);
		
	#total alignment analysis
	while( my($num_name, $value)=(each %total_num_out) ){
		E_RNA::refresh_log($variables{sample_log}, $num_name, $value);
	}
	
	#clear temporary files
	system("rm $sam_tmp $fasta_tmp");
	print  "\n Iterative sequences alignment of $sample_name using bowtie is done.\n\n\n";

}#1

##############################################################################
#
sub ncRNA_result_export{#1
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my @sample_names=split(",", $variables{sample_names});
	my $statistics_dir=$variables{dir_result}.'/'.'statistics';
	mkdir($statistics_dir, 0755) unless -d $statistics_dir;
	
	#get sample log information
	my (%statistics, @statistical_names);
	foreach my $sample_name(@sample_names){
		open my($IN), "<", $variables{dir_result}.'/'.$sample_name.'/'.$sample_name.'.log' or die;
		while(<$IN>){
			chomp($_);
			my($statistical_name, $value)=split("=", $_);
			$statistics{$sample_name}->{$statistical_name}=$value;
			push(@statistical_names, $statistical_name) unless List::Util::first {$_ eq $statistical_name} @statistical_names;
		}
		close($IN);	
	}
	
	print "\n\nIntegrate and export results into statistics.txt!\n\n\n";
	open my($STI), ">", $statistics_dir.'/statistics.txt' or die;
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
	my %data_frame;
	my $alignment_files_pointer=E_RNA::files_list($variables{dir_result}, 'incrusive_files', "alignment");
	my @alignment_files=@$alignment_files_pointer;
	foreach my $alignment_file(@alignment_files){#2
		print "$alignment_file\n";
		open my($IN), "<", $alignment_file or die;
		while(my $line=<$IN>){#4
			chomp($line);
			my ($sample_name, $index_name, $ref, $top_RC, $middle_RC, $bottom_RC, $seq_info)=split("\t", $line);
			$index_name='iterative' if $index_name=~/iterative/;
			$data_frame{$index_name}->{$ref}->{$sample_name}=$middle_RC;
		}#4
		close($IN);
	}#2
	
	foreach my $index_name(keys %data_frame){
		print "$index_name\n";
		my $pointer=$data_frame{$index_name};
		my %sample_df=%$pointer;
		
		#export %data_frame
		#my %hash=();
		my $frame_RC_file=$statistics_dir.'/df_totalRC_'.$index_name.'.txt';
		open my($RC), ">", $frame_RC_file or die;
			print $RC join("\t", 'Ref_names', @sample_names), "\n";
		my $frame_MRC_file=$statistics_dir.'/df_mappableRC_'.$index_name.'.txt';
		open my($MRC), ">", $frame_MRC_file or die;
			print $MRC join("\t", 'Ref_names', @sample_names), "\n";
		my $frame_IRC_file=$statistics_dir.'/df_insertRC_'.$index_name.'.txt';
		open my($IRC), ">", $frame_IRC_file or die;
			print $IRC join("\t", 'Ref_names', @sample_names), "\n";
		foreach my $ref_name(sort(keys %sample_df)){#3
			my (@RC, @MRC, @IRC);
			foreach my $sample_name(@sample_names){#4
				$sample_df{$ref_name}->{$sample_name}=0 unless exists $sample_df{$ref_name}->{$sample_name};
				my $RC=int($sample_df{$ref_name}->{$sample_name}+0.5);
				push(@RC, $RC );
				my $mappableRC=int( $RC * 1e6 / $statistics{$sample_name}->{$index_name.':mappable_reads_num'} +0.5 );
				push(@MRC, $mappableRC );
				my $insertRC=int( $RC * 1e6 / $statistics{$sample_name}->{'insert_reads_num'} +0.5 );
				push(@IRC, $insertRC );
			}#4
			print $RC  join("\t", $ref_name, @RC), "\n" if List::Util::max @RC>=$variables{read_counts_background};
			print $MRC  join("\t", $ref_name, @MRC), "\n" if List::Util::max @MRC>=$variables{read_counts_background};
			print $IRC  join("\t", $ref_name, @IRC), "\n" if List::Util::max @IRC>=$variables{read_counts_background};
		}#3
		close($RC);
		close($MRC);
		close($IRC);
	}#2


}#1

##################################################################
sub venn_diagram{
  my $variables_pointer=$_[0];
  my %variables=%$variables_pointer;
  my $pairs_pointer=$variables{venn_pairs_pointer};
  my @pairs=@$pairs_pointer;
  my $out_file=$variables{dir_result}.'/'.$_[1];

  # Create a communication bridge with R and start R
  my $R=Statistics::R->new();
  my @index_names=split(",", $variables{bowtie_index});
  my $index_num=@index_names;
  $R->set('index_num', $index_num);
  my $pairs_num=@pairs;
  $R->set('pairs_num', $pairs_num);
  $R->set("cutoff", $variables{venn_NormRC_cutoff});
  $R->set('out_file', $out_file);
  $R->run(qq`jpeg(out_file, units="mm", height=150*pairs_num, width=200*index_num,res=600, pointsize=10)`,
           q'par(mfrow=c(pairs_num,index_num), mar=c(2,2,2,2))',
           q'library(limma)',
         );

  foreach my $pointer(@pairs){#2 sample pair circyling
    my %pairs_hash=%$pointer;
    if ($pairs_hash{a} and $pairs_hash{b} and $pairs_hash{c}){#3 intersection of A, B and C
      foreach my $index_name(@index_names){#4
        my @A=split(",", $pairs_hash{a});
        my @B=split(",", $pairs_hash{b});
        my @C=split(",", $pairs_hash{c});
        $R->set("in_file", $variables{dir_result}.'/'.$index_name.".data_frame.txt");
        $R->set("index_name", $index_name);
        $R->set("A_samples", \@A);
        $R->set("B_samples", \@B);
        $R->set("C_samples", \@C);
        $R->set("A_name", $pairs_hash{a_name});
        $R->set("B_name", $pairs_hash{b_name});
        $R->set("C_name", $pairs_hash{c_name});
        $R->run(q'RNA<-read.delim(in_file, header=T)',
                q'A<-rowMeans(RNA[,A_samples])',
                q'A_name<-paste(sep="", A_name, "(", length(A_samples), ")=", length(A[A>=cutoff]) )',
                q'A<-(A>=cutoff)',
                q'B<-rowMeans(RNA[,B_samples])',
                q'B_name<-paste(sep="", B_name, "(", length(B_samples), ")=", length(B[B>=cutoff]) )',
                q'B<-(B>=cutoff)',
                q'C<-rowMeans(RNA[,C_samples])',
                q'C_name<-paste(sep="", C_name, "(", length(C_samples), ")=", length(C[C>=cutoff]) )',
                q'C<-(C>=cutoff)',
                q'ABC<-cbind(A, B, C)',
                q'ABC_counts<-vennCounts(ABC)',
                q'vennDiagram(ABC_counts, include="both", names=c(A_name, B_name, C_name), main=paste("Venn against", index_name, "mean of norm_RC>", cutoff), cex=2)',
             );
      }#4
    }#3
    elsif ($pairs_hash{a} and $pairs_hash{b}) {#3
      foreach my $index_name(@index_names){#4
        my @A=split(",", $pairs_hash{a});
        my @B=split(",", $pairs_hash{b});
        $R->set("in_file", $variables{dir_result}.'/'.$index_name.".data_frame.txt");
        $R->set("index_name", $index_name);
        $R->set("A_samples", \@A);
        $R->set("B_samples", \@B);
        $R->set("A_name", $pairs_hash{a_name});
        $R->set("B_name", $pairs_hash{b_name});
        $R->run(q'RNA<-read.delim(in_file, header=T)',
                q'A<-rowMeans(RNA[,A_samples])',
                q'A_name<-paste(sep="", A_name, "(", length(A_samples), ")=", length(A[A>=cutoff]) )',
                q'A<-(A>=cutoff)',
                q'B<-rowMeans(RNA[,B_samples])',
                q'B_name<-paste(sep="", B_name, "(", length(B_samples), ")=", length(B[B>=cutoff]) )',
                q'B<-(B>=cutoff)',
                q'AB<-cbind(A, B)',
                q'AB_counts<-vennCounts(AB)',
                q'vennDiagram(AB_counts, include="both", names=c(A_name, B_name), main=paste("Venn against", index_name, "mean of norm_RC>", cutoff), cex=2)',
             );
      }#4
    }#3
    else{ print "Wrong input in Venn diagramming!\n"; }

  }#2 sample pair circyling

  $R->run(q`dev.off()`);
  $R->stop();
  print "Intesection analysis among two or three samples is done!\n\n\n";
}#1


##############################################
sub clustering_analysis{
  my $variables_pointer=$_[0];
  my %variables=%$variables_pointer;
  my $index_name=$_[1];
  my $in_file=$_[1];
  my $out_file=$_[2];

  print "Clustering analysis begin:\n";
  my $R=Statistics::R->new();
  $R->set('in_file', $variables{dir_result}.'/'.$in_file);
  $R->set('out_file', $variables{dir_result}.'/'.$out_file);
  $R->run(q'df<-read.delim(in_file, header=T)',
          q'df_r<-t(df[,2:ncol(df)])',
          q'hc<-hclust(dist(df_r, method="euclidean"), method="ward")',
          qq'jpeg(out_file, units="mm", height=200, width=10*ncol(df),res=600, pointsize=10)',
          q'plot(hc)',
          q'plot(hc, hang = -1)',
          q'dev.off()',
         );
  $R->stop();
}

##############################################
sub RNA_categories{
  my $variables_pointer=$_[0];
  my %variables=%$variables_pointer;
  my $in_file=$variables{dir_result}.'/'.$_[1];
  my $out_file=$variables{dir_result}.'/'.$_[2];

  print "Graph pictures of RNA categories!\n";
  my $R=Statistics::R->new();
  my $sample_names_pointer=E_RNA::files_list($variables{dir_rawdata}, "sample_names");
  my @sample_names=@$sample_names_pointer;
  my $sample_num=@sample_names;
  $R->set('sample_num', $sample_num); 
  $R->set('in_file', $in_file);
  $R->set('out_file', $out_file);
  $R->run(qq`jpeg(out_file, units="mm", height=300, width=20*sample_num,res=600, pointsize=6)`,
           q'par(mfrow=c(4,1), mar=c(4,4,2,2))',
           q'df<-read.delim(in_file, header=T)',
           q'RNA_category_num<-(ncol(df)-10)/3',
           q'RNA_col<-RNA_category_num+5',
           q'barplot(as.matrix(df[,3:7]/1e6), legend=df[,2],  main="Relative distribution of RNAs in miRNA library", xlab="Reads", ylab="The number of reads (million)", beside=T)',
           q't_matrix<-t(as.matrix(df[,3:7]*100/df[,3]))',
           q'colnames(t_matrix)<-df[,2]',
           q'barplot(t_matrix, legend=names(df[,3:7]), main="Relative distribution of RNAs in miRNA library", ylim=c(0,100), xlab="sample names", ylab="The percent of reads versus the total number of raw reads %", beside=T)',
           q'barplot(as.matrix(df[,8:RNA_col]/1e6), legend=df[,2], main="Relative distribution of RNAs except human miRNA in miRNA library", ylab="The number of reads (million)", beside=T)',
           q't_matrix<-t(as.matrix(df[,8:RNA_col]*100/df[,3]))',
           q'colnames(t_matrix)<-df[,2]',
           q'barplot(t_matrix, legend=names(df[,8:RNA_col]), main="Relative distribution of RNAs except human miRNA in miRNA library", ylab="The percent of reads versus the total number of query reads %", beside=T)',
         );
  $R->run(q`dev.off()`);
  $R->stop();
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
#################################################################
#initiate and clear  the monitor.log and system_monitor.log
sub pre_initiate_monitor_log{
	my ($variables_pointer, $total_beginning_time_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my @total_beginning_time=@$total_beginning_time_pointer;
	
	#initiate and clear  the monitor.log
	open my($SYS), ">", $variables{file_system_monitor_log} or die;
		print $SYS join("\t", 'Time', 'Duration(s)', 'CPU_usage(%)', 'Memory_usage(%)'), "\n";
	close($SYS);

	open my($MON), ">", $variables{file_monitor_log} or die; 
	print $MON "Total:beginning_time=", join(",", @total_beginning_time), "\n";
	print $MON "Total:ending_time=NA\n";
	print $MON "Total:running_time=NA\n";
	print $MON "Total:supposed_time=NA\n";
	print $MON "Total:beginning_seconds=", time, "\n";
	foreach my $sample_name(sort (keys %sample_info) ) {#2
		my $supposed_time=int($sample_info{$sample_name}->{'files_size'}/228383); #unit is second
		print $MON $sample_name, ":supposed_time=$supposed_time\n";
		print $MON $sample_name, ":beginning_time=NA\n";
		print $MON $sample_name, ":ending_time=NA\n";
		print $MON $sample_name, ":running_time=NA\n";
	}#2
	close($MON);
}

#################################################################
#refresh  the monitor.log
sub pre_refresh_monitor_log{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	
	#get all sample_log_files
	my $log_files_pointer=E_RNA::files_list($variables{dir_log}, 'files');
	my @log_files=@$log_files_pointer;
	#refresh time_monitor_file
	foreach my $sample_log_file(@log_files){#2
		my @a=split('/', $sample_log_file);
		my $sample_name=$a[-1];
		$sample_name=~s/\.log//;
		open my($IN), "<", $sample_log_file or die;
		while(my $line=<$IN>){#3
			chomp($line);
			my($name, $value)=split("=", $line);
			$value=~s/\"//g;
			if($line=~/beginning_time|ending_time|running_time/){#4
				E_RNA::refresh_log($variables{file_monitor_log}, $sample_name.':'.$name, $value);
			}#4
		}#3
		close($IN);
	}#2

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



1;  # make sure the file returns true or require will not succeed!#






