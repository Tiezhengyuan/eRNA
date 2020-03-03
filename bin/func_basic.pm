#! /usr/bin/perl -w
use strict;
use warnings;
use Bio::SeqIO;
use List::Util;
use File::Find;

require "func_common.pm"; #sub_common::

####################################################
#

#get the directory of perl scripts involved in Pscore
#our $perl_dir=Cwd::getcwd();

#the file constains all subroutines required for running E_RNA_pipeline
package sub_basic;

#require $perl_dir."/functions_common.pm"; #sub_common::
######################################################

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
#############################################
#type= 'replace' or 'add' or 'append', the default is replace
sub refresh_log{
	my ($log_file, $r_name, $r_value, $type)=@_;
	$type='replace' unless $type;
	
	my %variables;
	#read old data
	if(-f $log_file){
		my $variables_pointer=sub_common::file_to_hash($log_file, "=");
		%variables=%$variables_pointer;
	}
	
	#refresh new data
	if ($variables{$r_name}){#2
		if ($type eq 'add'){ #add
				$variables{$r_name}=$variables{$r_name}+$r_value;
		}
		elsif ($type eq 'replace'){ #replace
				$variables{$r_name}=$r_value;
		}
		else{#append
				$variables{$r_name}=$variables{$r_name}.$type.$r_value;
		}
	}#2
	else{#2
		$variables{$r_name}=$r_value;
	}#2
	
	#export refreshed data
	sub_common::hash_to_file(\%variables, $log_file, "=");
	
}
#######################
sub refresh_log_hash{
	my ($log_file, $hash_pointer, $type)=@_;
	my %hash=%$hash_pointer;
	#
	while( my($name, $value)=(each %hash) ){
		sub_basic::refresh_log($log_file, $name, $value, $type);
	}
	
}
############################################
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


####################################33
#function: raw data files list
#list_type: incrusive_files
sub SM_rawfiles{#1
	my ($dir_str)=@_;
	my @rawdata_dir_arr=split(',', $dir_str);
	
	#print "The file list of raw data:\t";
	my @raw_files;
	for my $rawdata_dir(@rawdata_dir_arr){#2
		print "raw data directory: $rawdata_dir\n"; 
		#incrusive files and directories
		my $files_pointer=sub_common::files_list($rawdata_dir, 'incrusive_file');
		my @incrusive_files=@$files_pointer;
		#get raw files
		my @sub_files=grep(/\.fastq$|\.fq$/i, @incrusive_files);
		@raw_files=(@raw_files, @sub_files);
	}#2
	#
	my $raw_files_num=@raw_files;
	print "$raw_files_num files of raw data will be analyzed.\n\n";

	return(\@raw_files);
}#1

###########
sub SM_rawfiles_R1{#1
	my ($dir_str)=@_;
	
	#get all fastq files
	my $files_pointer=SM_rawfiles($dir_str);
	my @raw_files=@$files_pointer;
	
	#R1 files
	my @raw_files_input=grep(/R1_|_R1/, @raw_files) ;
	@raw_files_input=@raw_files if @raw_files_input==0;
	#print "@raw_files_input\n";
	my $raw_files_num=@raw_files_input;
	print "$raw_files_num files of raw data will be analyzed.\n\n";

	return(\@raw_files_input);
}#1



###################################################################
#function
#raw_file ~ raw file name
sub SM_sample_from_raw{
	my $raw_files_pointer=$_[0];
	my @raw_files=@$raw_files_pointer;
	#get only R1
	my @raw_files_input=grep(/R1_|_R1/, @raw_files) ;
	@raw_files_input=@raw_files if @raw_files_input==0;
	
	#
	my %samples;
	foreach my $raw_file(@raw_files_input){
		my $sample_name=sub_common::file_operation($raw_file, 'name_head');
		$samples{$raw_file}=$sample_name;
	}
	
	return(\%samples);
}
###################################################################
#function
sub SM_rawfile_sample{
	my $sample_info_pointer=$_[0];
	my %sample_info=%$sample_info_pointer;
	
	my %raw_to_sample;
	foreach my$sample_name(keys %sample_info){#2
		my $pointer=$sample_info{$sample_name}->{R1_files};
		my @sample_raw_files_input=@$pointer;
		foreach my $raw_file(@sample_raw_files_input){
			$raw_to_sample{$raw_file}=$sample_name;
		}
	}#2
	#
	return(\%raw_to_sample);
}
####################################
sub SM_assign_sample_dir{
	my($dir_result_array, $sample_name)=@_;
	my @result_dir_arr=split(',', $dir_result_array);
	
	my $sample_dir=$result_dir_arr[0];
	my $sample_dir_space=0; #free space
	my $sample_dir_old=0; # judge is have old sample_dir in a given sample_dir
	foreach my $dir(@result_dir_arr){#2
		$dir=sub_common::format_directory($dir);
		my $result_free_space=sub_common::free_space($dir); # by GB
		my $tmp_sample_dir=$dir.$sample_name;
		if(-d $tmp_sample_dir){
			$sample_dir=$tmp_sample_dir;
			last;
		}
		elsif ($result_free_space>$sample_dir_space){
			$sample_dir=$tmp_sample_dir;
			$sample_dir_space=0;
		}
	}#2
	#
	$sample_dir=sub_common::format_directory($sample_dir);
	return($sample_dir);
}
######################
#Status of log file: no, on, off
sub combine_log_files{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my @sample_names=split(",", $variables{sample_names});
	
	my $total_parts_status=0;
	my $found=1; #judge if running time;
	#print "Fresh sample log files in $variables{sample_log_dir}\n";
	foreach my $sample_name(@sample_names){#2
		my $log_files_pointer=sub_common::copy_files($variables{dir_result}.$sample_name, $variables{dir_log}.$sample_name, 'log');
		my @log_files=@$log_files_pointer;
		
		#initiate hash
		my %hash;
		$hash{sample_name}=$sample_name;
		$hash{status}='no';
		
		#refresh hash
		my @log_status;
		my $parts_status=0;
		foreach my $log_file(@log_files){#3
			my $log_variables_pointer=sub_common::file_to_hash($log_file, '=');
			my %log_variables=%$log_variables_pointer;
			$log_variables{status}='no' unless exists $log_variables{status};
			
			unless($log_variables{status} eq 'no'){#4
				foreach my $key(keys %log_variables){#5
					if(exists $hash{$key} and $key eq 'beginning_time'){
						$hash{$key}=get_time($hash{$key}, $log_variables{$key},'early');
					}
					elsif(exists $hash{$key} and $key eq 'ending_time'){
						$hash{$key}=get_time($hash{$key}, $log_variables{$key},'late');
					}
					elsif(exists $hash{$key} and $key=~/_num$/){
						$hash{$key} += $log_variables{$key};
					}
					elsif($key eq 'status'){
						push(@log_status, $log_variables{$key}) unless List::Util::first {$_ eq $log_variables{$key}} @log_status;
						$parts_status++ if $log_variables{$key} eq 'off';
					}
					else{
						$hash{$key}=$log_variables{$key};
					}
				}#5
			}#4
		}#3
		$hash{'supposed_time'}=$sample_info{$sample_name}->{'supposed_time'};
		$hash{'parallel_parts'}=$sample_info{$sample_name}->{'parallel_parts'};
		$hash{'parts_status'}=$parts_status;
		$total_parts_status += $parts_status;
		
		#status by samples
		if(@log_status==1 and $log_status[0] eq 'off'){
			$hash{status}='off';
			$found=0;
			$hash{'running_time'}=sub_common::get_time($hash{'beginning_time'}, $hash{'ending_time'},'duration');
		}
		elsif(List::Util::first {$_ eq 'on'} @log_status){
			$hash{status}='on';
			$hash{'ending_time'}='NA';
			$hash{'running_time'}='NA';
		}
		
		#export $hash;
		my $sample_log=$variables{dir_log}.$sample_name.'.log';
		sub_common::hash_to_file(\%hash, $sample_log, '=');
	}#2
	#refresh the file Total.log
	sub_basic::refresh_log($variables{file_total_log}, 'parts_status', $total_parts_status);
	
	return($found);
}

###################################################################
#the centerpiece function for sample management
sub SM_sample_info{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	
	#get files of raw data
	my $R1files_pointer=sub_basic::SM_rawfiles_R1($variables{dir_raw_data});
	#print @$R1files_pointer;
	#read sample_info file
	my $raw_to_sample_pointer=sub_common::file_to_hash($variables{file_sample_info}, ',');
	#sub_data::print_hash($raw_to_sample_pointer);
	#match fastq and sample name
	my $sample_rawfiles_pointer=sub_basic::SM_match_sample($R1files_pointer, $raw_to_sample_pointer);
	my %sample_rawfiles=%$sample_rawfiles_pointer;
	print "\n######The file list of raw data:\n";
	#sub_data::print_hash($sample_rawfiles_pointer);

	#sample info
	my $n=1;
	my %sample_info;
	while (my($sample_name, $R1files_str)=each (%sample_rawfiles) ){#2
		printf("[%d]  %s:\n", $n, $sample_name);
		#1: R1_files
		$sample_info{$sample_name}->{'R1_files'}=$R1files_str;
		my @sample_R1_files=split(',', $R1files_str);
		print "\t$_\n" for @sample_R1_files;
		if($variables{'sequencing_end'} eq 'double'){
			my @sample_R2_files=map {$_=~s/R1/R2/} @sample_R1_files;
			my $R2files_str=join(',', @sample_R2_files);
			$sample_info{$sample_name}->{'R2_files'}=$R2files_str;
			print "\t$_\n" for @sample_R2_files;
		}
		#2: file size and supposed time
		my $files_size=0;
		foreach my $raw_file(@sample_R1_files){
			$files_size += sub_common::file_operation($raw_file, 'file_size');
		}
		$sample_info{$sample_name}->{'files_size'}=$files_size;
		$sample_info{$sample_name}->{'supposed_time'}=int($files_size/228383); #unit is second
		$sample_info{$sample_name}->{'parallel_parts'}=@sample_R1_files;
		$sample_info{$sample_name}->{'sample_dir'}=$variables{'dir_result'}.$sample_name.'/';
		$n++;
	}
	print "######################\n\n";
	#
	my @sample_names=sort( keys %sample_info);
	$variables{sample_names}=join(",", @sample_names);
	$variables{sample_info_pointer}=\%sample_info;
	
	#raw data files to sample names
	$variables{raw_to_sample_pointer}=$raw_to_sample_pointer;
	#print %raw_to_sample;
	return(\%variables);
}


###################################################################
#sample name ~ R1 files
sub SM_match_sample{
	my ($raw_files_pointer, $raw_to_sample_pointer)=@_;
	my %raw_to_sample=%$raw_to_sample_pointer;
	my @raw_files=@$raw_files_pointer;
	my @raw_files_input=grep(/R1_|_R1/, @raw_files) ;
	@raw_files_input=@raw_files if @raw_files_input==0;
	
	#
	my %sample_rawfiles;
	foreach my$raw_file(@raw_files){#2
		#print "$raw_file\n";
		if(exists $raw_to_sample{$raw_file}){#3
			my $sample_name=$raw_to_sample{$raw_file};
			if (exists $sample_rawfiles{$sample_name}){
				$sample_rawfiles{$sample_name} .= ','.$raw_file;
			}
			else{
				$sample_rawfiles{$sample_name}=$raw_file;
			}
		}#3
	}#2
	return(\%sample_rawfiles);
}


##########################
sub auto_sample_names{#1
	my ($dir, $file_format)=@_;
	$dir .= '/' unless $dir=~/\/$/;

	my (%sample_size, @sample_names);
	if ($file_format eq 'FASTQ'){
		File::Find::find(sub {
				my $file=$File::Find::name;
				if (-f $file and $file=~/\.fastq$|\.fq$/){#4
					my @array1=split("/", $file);
					my $file_name=$array1[-1];
					$file_name=~s/\.fastq$|\.fq$//;
					my @array2=split("_", $file_name);
					my $sample_name;
					if (@array2==1){	$sample_name=$file_name;	}
					elsif (@array2==2){	$sample_name=$array2[0]."_".$array2[1];	}
					else{	$sample_name=$array2[0]."_".$array2[1]."_".$array2[2];	}
					my @args=stat($file);
					$sample_size{$sample_name} += $args[7];
				}#4
		}, $dir);
	}
	elsif ($file_format eq 'SAM'){
		File::Find::find(sub {
				my $file=$File::Find::name;
				if (-f $file and $file=~/\.sam$/i){#4
					my @array1=split("/", $file);
					my $file_name=$array1[-1];
					$file_name=~s/\.sam$//i;
					my @array2=split("_", $file_name);
					my $sample_name;
					if (@array2==1){	$sample_name=$file_name;	}
					elsif (@array2==2){	$sample_name=$array2[0]."_".$array2[1];	}
					else{	$sample_name=$array2[0]."_".$array2[1]."_".$array2[2];	}
					my @args=stat($file);
					$sample_size{$sample_name} += $args[7];
				}#4
		}, $dir);
	}
	elsif ($file_format eq 'BAM'){
		File::Find::find(sub {
				my $file=$File::Find::name;
				if (-f $file and $file=~/\.bam$/i){#4
					my @array1=split("/", $file);
					my $file_name=$array1[-1];
					$file_name=~s/\.bam$//i;
					my @array2=split("_", $file_name);
					my $sample_name;
					if (@array2==1){	$sample_name=$file_name;	}
					elsif (@array2==2){	$sample_name=$array2[0]."_".$array2[1];	}
					else{	$sample_name=$array2[0]."_".$array2[1]."_".$array2[2];	}
					my @args=stat($file);
					$sample_size{$sample_name} += $args[7];
				}#4
		}, $dir);
	}
	#order output by file size
	foreach my $sample_name(sort { $sample_size{$b}<=>$sample_size{$a} } (keys %sample_size) ){
		push(@sample_names, $sample_name);
	}
	return(\@sample_names);
}#1

############################################################
sub check_bowtie_index{
	my ($dir_bowtie, $version)=@_;
	#
	my %index_info;
	$index_info{bowtie1}->{tails}=join(',', '.1.ebwt', '.2.ebwt','.3.ebwt','.4.ebwt','.rev.1.ebwt','.rev.2.ebwt');
	$index_info{bowtie1}->{script}=$dir_bowtie.'bowtie-build';
	$index_info{bowtie2}->{tails}=join(',', '.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2');
	$index_info{bowtie2}->{script}=$dir_bowtie.'bowtie2-build';
	
	#get fasta files
	my $files_pointer=sub_common::files_list($dir_bowtie, 'file');
	my @files=@$files_pointer;
	my @fasta_files=grep(/\.fa$|\.fasta$/, @files); 
	#print @fasta_files;
	
	foreach my $fa_file(@fasta_files){#2
		my $bowtie_index_name=sub_common::file_operation($fa_file, 'name_head');
		my $bowtie_index=$dir_bowtie.$bowtie_index_name;
		my @tails=split(',', $index_info{$version}->{tails});
		#check bowtie index
		my $found=1;
		foreach (@tails){
			$found=0 unless -f $bowtie_index.$_;
		}
		#build bowtie index if no found
		my $script=$index_info{$version}->{script};
		system("$script $fa_file $bowtie_index") if $found==0;
	}#2

}


###########################################
#initiate the monitor log file, the system log file, and the sample log files
sub initiate_log_files{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my @sample_names=split(',', $variables{sample_names});
	my $total_beginning_time= join(",", localtime(time) );
	
	#initiate time_monitor.log
	my $total_supposed_time=0;
	my $total_parallel_parts=0;
	print "initiate and clear $variables{file_time_monitor_log}\n";
	open my($MON), ">", $variables{file_time_monitor_log} or die; 
	print $MON join("\t", 'sample_names', 'supposed_time', 'beginning_time', 'ending_time', 'running_time', 'parallel_parts', 'parts_status', 'status'), "\n";
	foreach my $sample_name(@sample_names) {
		print $MON join("\t", $sample_name, $sample_info{$sample_name}->{supposed_time}, 'NA', 'NA', 'NA', 
							$sample_info{$sample_name}->{parallel_parts}, 0, 'no'), "\n";
		$total_supposed_time += $sample_info{$sample_name}->{supposed_time};
		$total_parallel_parts += $sample_info{$sample_name}->{parallel_parts};
	}
	print $MON join("\t", 'Total', $total_supposed_time, $total_beginning_time, 'NA', 'NA', $total_parallel_parts, 0, 'no'), "\n";
	close($MON);
	
	print "initiate $variables{file_total_log}\n";
	sub_basic::refresh_log($variables{file_total_log}, 'beginning_time', $total_beginning_time);
	sub_basic::refresh_log($variables{file_total_log}, 'supposed_time', $total_supposed_time);
	sub_basic::refresh_log($variables{file_total_log}, 'parallel_parts', $total_parallel_parts);
	sub_basic::refresh_log($variables{file_total_log}, 'parts_status', 0);
	sub_basic::refresh_log($variables{file_total_log}, 'status', 'on');

	#initiate and clear the monitor.log
	print "initiate and clear $variables{file_system_monitor_log}\n";
	open my($SYS), ">", $variables{file_system_monitor_log} or die;
	print $SYS join("\t", 'Time', 'Duration(s)', 'CPU_usage(%)', 'Memory_usage(%)'), "\n";
	close($SYS);
	
	#initiate sample log files
	foreach my$sample_name(@sample_names){
		my $sample_dir=sub_common::format_directory($variables{dir_result}.$sample_name);
		my $files_pointer=sub_common::files_list($sample_dir,'file');
		my @files=@$files_pointer;
		my @log_files=grep(/\.log$/, @files);
		foreach my $log_file(@log_files){
			#refresh_log($log_file, 'status', 'no');
		}
	}
	#
	return($total_beginning_time);
}
###################
sub refresh_monitor_log{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my @sample_names=split(',', $variables{sample_names});
	my @monitor_items=('beginning_time', 'ending_time', 'running_time', 'parts_status', 'status');
	
	#read the file time_monitor.log
	my $time_pointer=sub_common::file_to_hash2($variables{file_time_monitor_log});
	my %time=%$time_pointer;
	
	#refresh by reading sample log
	foreach my $log_name(keys %time) {#2
		my $log_pointer=sub_common::file_to_hash($variables{dir_log}.$log_name.'.log', '=');
		my %log=%$log_pointer;
		$log{status}='no' unless exists $log{status};
		if($log{status} eq 'on' or $log{status} eq 'off'){#3
			foreach my $name(@monitor_items){
				$time{$log_name}->{$name}=$log{$name} if exists $log{$name};
			}
		}#3
	}#2
	
	#update
	sub_common::hash2_to_file(\%time, $variables{file_time_monitor_log});
	#
	return(\%time);
}

######################
#initiate start/end time and sample directory
sub initiate_starting_time{
	my ($variables_pointer, $raw_file)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my $raw_to_sample_pointer=$variables{raw_to_sample_pointer};
	my %raw_to_sample=%$raw_to_sample_pointer;
	
	#refresh sample_dir
	$variables{raw_file}=$raw_file;
	$variables{sample_name}=$raw_to_sample{$raw_file}; #sample name
	$variables{sample_dir}=$sample_info{$variables{sample_name}}->{sample_dir};
	#Note: The result directory is temporarily changed to $variables{sample_dir} for $variables{sample_name}

	#refresh sample_log
	$variables{beginning_time}=join(',', localtime(time) );
	my $raw_file_name_head=sub_common::file_operation($raw_file, 'name_head');
	$variables{sample_out_file}=$variables{sample_dir}.$raw_file_name_head;
	$variables{sample_log}=$variables{sample_out_file}.'.log';
	$variables{supposed_time}=int($sample_info{$variables{sample_name} }->{'supposed_time'}/$sample_info{$variables{sample_name} }->{'parallel_parts'}); #unit is second
	refresh_log($variables{sample_log}, "sample_name", $variables{sample_name});
	refresh_log($variables{sample_log}, "beginning_time" , $variables{beginning_time});
	refresh_log($variables{sample_log}, "supposed_time" , $variables{supposed_time});
	refresh_log($variables{sample_log}, "status" , 'on');
	
	return(\%variables);
}

##########
sub initiate_ending_time{
	my ($log_file, $beginning_time_str)=@_;

	my $ending_time_str=join(',', localtime(time) );
	my $running_time=sub_common::get_time($beginning_time_str, $ending_time_str, 'duration');
	refresh_log($log_file, "ending_time" , $ending_time_str);
	refresh_log($log_file, "running_time", $running_time);
	refresh_log($log_file, "status", 'off');
}


###############
sub scaling_normalization{
	my($hash2_pointer, $raw_reads_num)=@_;
	my %hash2=%$hash2_pointer;
	
	my %norm_hash2;
	foreach my $row(keys %hash2){
		my $hash_pointer=$hash2{$row};
		my %hash=%$hash_pointer;
		foreach my $col(keys %hash){
			my $counts=$hash2{$row}->{$col};
			$norm_hash2{$row}->{$col}=int(($counts*1e6)/$raw_reads_num + 0.5);
		}
	}

	return(\%norm_hash2);
}


################################
#both $a and $b are gene_ids from RNA-seq
#first is intersections, and second is union if no intersections
sub most_intersections{
	my($a_str, $b_str)=@_;
	
	#intersections
	my $ids;
	if ($a_str eq 'NA'){
		$ids=$b_str;
	}
	else{
		my @a_arr=split(',', $a_str);
		my @b_arr=split(',', $b_str);
		my @most;
		foreach my $a(@a_arr){
			push(@most, $a) if List::Util::first {$_ eq $a} @b_arr;
		}
		#
		@most=(@a_arr, @b_arr) if @most==0;
		$ids=join(',', @most);
	}
	return($ids);
}
################################
#both $a and $b are gene_ids from RNA-seq
#first is union, and second is 'NA';
sub most_unions{
	my($pointer)=@_;
	my @arr=@$pointer;
	
	my $references='NA';
	if(@arr>0){
		my @refs;
		my @ids=split(',', join(',', @arr));
		foreach my $id(@ids){
			push(@refs, $id) unless List::Util::first {$_ eq $id} @refs;
		}
		$references=join(',', @refs);
	}
	
	return($references);
}

################
#-- fetch file from url
sub download_url_file{
		my($url, $file_name, $local_dir)=@_;
		$url .= '/' unless $url=~/\/$/;
		$local_dir .= '/' unless $local_dir=~/\/$/;
  
		#fetch file
		my $url_file= $url.$file_name;
		my $local_file=$local_dir.$file_name;
		if (-f $local_file){#2
				my $status = getstore($url_file, $local_file);
				#report status
				if ( is_success($status) ){#3
						print "file downloaded correctly\n";
						#unzip file
						system("gunzip $local_file") if $local_file=~/\.gz$/;
				}#3
				else{#3
						print "error downloading file: $status\n";
				}#3
		}#2
		else{
			print "$local_file exists. So skip downloading.\n";
		}
  #return
}

#############################
1;  # make sure the file returns true or require will not succeed!#

