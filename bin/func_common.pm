#! /usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use File::Find;
use File::Fetch;
use LWP::Simple;
#use constant false => 0;
#use constant true  => 1;

####################################################
#
#the file constains all subroutines required for running E_RNA_pipeline
package sub_common;

##############################


########################
#
#key vs value must be one on one or one on multiple
sub file_to_hash{
	my ($file, $sep, $key_index, $value_index)=@_;
	$sep='=' unless $sep;
	
	my %hash;
	if(-f $file){#2
		open my ($INFO), "<", $file or die;
		while (<$INFO>) {#3
			chomp($_);
			#print"$key:$value\n";
			if ($_=~/$sep/){#4
				#seperate line
				my($key,$value);
				if ($key_index and $value_index){
					my @items = split(/$sep/, $_); #split on the tabs
					$key=$items[$key_index]; #
					$value=$items[$value_index];
				}else{
					($key, $value) = split(/$sep/, $_, 2); #split into two parts only
				}
				#update hash
				if(exists $hash{$key}){ $hash{$key} .= ','.$value ; }
				else{ $hash{$key} =$value ;  }
			}#4
		}#3
		close($INFO);
	}#2
	
	return(\%hash);
}

##########################
#read file as the hash with two levels
sub file_to_hash2{
	my($file, $sep)=@_;
	$sep="\t" unless $sep;
	
	#
	my %hash2;
	open my($OUT), "<", $file or die;
	my $header=<$OUT>;
	chomp($header);
	#print "####$header###\n";
	my @col_names=split(/$sep/,$header);
	#
	while (<$OUT>){
		chomp($_);
		my @items=split(/$sep/, $_);
		my $row_name=$items[0];
		for(my $i=1;$i<@items; $i++){
			my $col_name=$col_names[$i];
			$hash2{$row_name}->{$col_name}=$items[$i];
			#print "$row_name:$col_name=$items[$i]\n";
		}
	}
	close($OUT);
	
	return(\%hash2);
}

##########################
#read file as the hash with two levels
#A file has at least two columns, no header no rownames
#records would be multiple vs. multiple
sub fileflat_to_hash2{
	my($file, $sep)=@_;
	$sep="\t" unless $sep;
	
	#
	my %hash2;
	open my($OUT), "<", $file or die;
	while (<$OUT>){
		chomp($_);
		my @items=split(/$sep/, $_);
		my $key1=shift @items;
		my $key2=shift @items;
		$hash2{$key1}->{$key2}=(@items>0) ? join($sep, @items) : "";
	}
	close($OUT);
	
	return(\%hash2);
}

########################
sub file_to_multihash{
	my($file, $sep)=@_;
	$sep="\t" unless $sep;
	#
	my %hash2;
	open my($OUT), "<", $file or die "can't open $file\n";
	my $header=<$OUT>;
	chomp($header);
	#print "####$header###\n";
	my @names=split(/$sep/,$header);
	#read lines
	while (<$OUT>){
		chomp($_);
		my @items=split(/$sep/, $_);
		my $key1_name=$items[0];
		my $key2_name=$items[1];
		for(my $i=2;$i<@items; $i++){
			my $key_name=$names[$i];
			$hash2{$key1_name}->{$key2_name}->{$key_name}=$items[$i];
		}
	}
	close($OUT);
	
	return(\%hash2);
}
#############################
#export the hash with two levels into text file
sub hash2_to_file{
	my($hash2_pointer, $statistics_file, $sep, $first_col_name)=@_;
	my %hash2=%$hash2_pointer;
	$sep="\t" unless defined $sep;
	$first_col_name ='names' unless defined $first_col_name;
	
	#get keys1 and key2
	my @row_names= sort keys %hash2;
	my %names;
	my @col_names=map { keys %{$hash2{$_}} } @row_names;
	@col_names= sort sub_data::array_unique(@col_names);
	
	open my($OUT), ">", $statistics_file or die;
	print $OUT join($sep, $first_col_name, @col_names), "\n";
	foreach my $row_name(@row_names){#2
		my @counting_num;
		foreach my $col_name(@col_names){
			my $value= (exists $hash2{$row_name}->{$col_name}) ? $hash2{$row_name}->{$col_name} : 0;
			push(@counting_num, $value);
		}
		print $OUT join($sep, $row_name, @counting_num), "\n";
	}#2
	close($OUT);
	#
}
########################################################################
#function: the specific files list in the specific direcotry
#list_type: files, incrusive_file, incrusive_file_name, file_name, dir, incrusive_dir 
sub files_list{#1
	my ($dir, $list_type)=@_;
	$dir .= '/' unless $dir=~/\/$/;
	my (@file_names, @dir_names, @incrusive_dir, @incrusive_files);
	#print "$dir\n";
	
	# relative file name and directory name
	opendir(DIR, $dir);
	my @all=readdir(DIR);
	closedir(DIR);
	foreach (@all){
		unless ($_ =~ /^\.|~/){
			if (-d $dir.$_) {	push(@dir_names, $_) ;		}
			else {	push(@file_names, $_) ;		}
		}
	}
	#incrusive files and directories
	File::Find::find(sub {
			my $name=$File::Find::name;
			push(@incrusive_files, $name) if -f $name;
			push(@incrusive_dir, $name) if -d $name;
	}, $dir);

	my @files;
	if ($list_type eq 'dir_name'){
		@files=@dir_names;
		#print "$list_type===@files\n";
	}
	elsif ($list_type eq 'dir'){
		@files=map{$dir.$_} @dir_names;
	}
	elsif ($list_type eq 'incrusive_dir'){
		@files=@incrusive_dir;
	}
	elsif ($list_type eq 'file_name'){
		@files=@file_names;
	}
	elsif ($list_type eq 'file'){
		@files=map {$dir.$_} @file_names;
	}
	elsif ($list_type eq 'incrusive_file'){
		@files=@incrusive_files;
	}
	elsif ($list_type eq 'incrusive_file_name'){
		foreach (@incrusive_files){
			my @array=split("/", $_);
			my $name=$array[-1];
			push(@files, $name);
		}
		#unique files
		@files = sub_data::array_unique(@files);
	}

	return(\@files);
}#1

##############################
#
sub file_operation{
	my ($file, $type)=@_;

	#
	my @array=split("/", $file);
	my $file_name=$array[-1];
	my @items=split(/\./, $file_name);
	my $file_name_head=$file_name;
	my $file_name_tail="";
	if (@items>1){
		$file_name_tail=pop @items;
		$file_name_head=join('.', @items);
	}
	
	my $out;
	if($type eq 'file_name'){# file name
		$out=$file_name;
	}
	elsif($type eq 'name_head'){# file name head
		$out=$file_name_head;
	}
	elsif($type eq 'name_tail'){# file name tail
		$out=$file_name_tail;
	}
	elsif($type eq 'file_head'){# file head
		pop @array;
		$out=join('/', @array, "").$file_name_head;
	}
	elsif($type eq 'file_size'){# file size
		my @stats=stat($file);
		$out=$stats[7];
	}
	elsif($type eq 'directory'){# directory
		pop @array;
		$out=join('/', @array, "");
	}
	return($out);
}

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


#########################
sub copy_files{
	my ($in_dir, $out_dir, $file_tail)=@_;
	$in_dir .= '/' unless $in_dir=~/\/$/;
	$out_dir .= '/' unless $out_dir=~/\/$/;
	mkdir($out_dir, 0755) unless -d $out_dir;
	
	my $file_names_pointer=files_list($in_dir,'file_name');
	my @file_names=@$file_names_pointer;
	my @sub_file_names=grep(/\.$file_tail$/, @file_names);
	#
	my @out_files;
	if(@sub_file_names>0){#2
		foreach my$file_name(@sub_file_names){#3
			my $in_file=$in_dir.$file_name;
			my $out_file=$out_dir.$file_name;
			system("cp $in_file $out_file"); 
			push(@out_files, $out_file);
		}#3
	}#2
	
	return(\@out_files);
}
############################################
#store data from array into a file
#######
sub array_to_file{
	my($arr_pointer, $outfile, $sep)=@_;
	my @arr=@$arr_pointer;
	$sep="\n" unless $sep; #default is = as split
	
	#
	open my($OUT), ">", $outfile or die;
	foreach my $line(@arr){
		print $OUT $line, $sep; 
	}
	close($OUT);

}

############################################
#######
sub hash_to_file{
	my($hash_pointer, $out_file, $sep)=@_;
	my %hash=%$hash_pointer;
	$sep="=" unless $sep; #default is = as split
	
	open my($OUT), ">", $out_file or die "can't open $out_file\n";
	foreach my$key(sort (keys %hash)){
		print $OUT join("", $key, $sep, $hash{$key}), "\n";
	}
	close($OUT);
}

#################################
#calculate free space of a given directory
sub free_space{
	my($dir)=@_;
	#
	open my ($SPA), "df $dir | " or die;
	my $result_free_space=0;
	my($l1, $l2);
	while( defined($l1=<$SPA>) && defined($l2=<$SPA>) ){
		chomp($l1, $l2);
		my ($disk_id, $total_space, $used_space, $free_space, $space_ratio, $mount_dir)=split(" ", $l2);
		$result_free_space=int($free_space/(1024**2)+0.5);
		#print "$result_free_space\n"; #return GB
	}
	close($SPA);
	return($result_free_space);
}
########################
sub check_free_space {
	my ($raw_data_dir, $result_dir)=@_;
	
	my (@incrusive_raw_data_files);
	my $raw_data_size=0;
	#incrusive files
	File::Find::find(sub {
			my $name=$File::Find::name;
			push(@incrusive_raw_data_files, $name) if -f $name;
		}, $raw_data_dir);
	#raw data size
	foreach my $raw_file(@incrusive_raw_data_files){
		my @stats=stat($raw_file);
		$raw_data_size+= $stats[7];
	}
	my $supposed_result_space=int($raw_data_size*3/(1024**3)+0.5);
	#
	my $result_free_space=sub_common::free_space($result_dir);
	my $space_info="Free space of $result_dir is $result_free_space GB, and $supposed_result_space GB might be used.\n";
	return($space_info);
}
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

#################
sub system_loading{
	
	#return load average
	open(my $LOAD, "uptime |") or die;
		my $line=<$LOAD>;
		chomp($line);
		my @array=split(" ", $line);
		my $load=$array[-1];
	close($LOAD);
	
	#return cpu usage and memory usage
	my $tmp_cpu=0;
	my $tmp_memory=0;
	open(my $PROC, "ps u |") or die;
	while (my $line=<$PROC>){
		chomp($line);
		my @array=split(" ", $line);
		if ($array[2] cmp '%CPU' and $array[2]>20){
			$tmp_cpu +=$array[2];
			$tmp_memory +=$array[3];
		}
	}
	close($PROC);
	
	my $out=join(',', $load, $tmp_cpu, $tmp_memory);
	return($out);
}

######################
#######3
sub combine_files{
	my ($sample_dir, $file_tail, $out_file_head)=@_;
    
    #get files before combination
	my $files_pointer=sub_common::files_list($sample_dir, 'file');
	my @files=@$files_pointer;
	my @sub_files=grep(/$file_tail$/, @files);
	
	#combine files
	my $in_sub_files_str=join (' ', @sub_files);
	my $out_file=$out_file_head.'.'.$file_tail;
	system("cat $in_sub_files_str > $out_file") if @sub_files>0; #combine
	
}

################################################################
#get date and time
sub get_time{
	my ($begin_time_str, $end_time_str, $type)=@_;
	$end_time_str='0,0,0,0,0,0,0,0,0' if $end_time_str eq '0';
	$type='duration' unless $type;
	
	my $my_time=0;
	my ($sec1,$min1,$hour1,$monthday1,$month1,$year1,$weekday1,$yearday1,$isdaylight1)=split(',', $begin_time_str);
	$year1 += 1900;
	my $begin_time=$year1.'/'.$month1.'/'.$monthday1.', '.$hour1.':'.$min1.':'.$sec1;	
	my ($sec2,$min2,$hour2,$monthday2,$month2,$year2,$weekday2,$yearday2,$isdaylight2)=split(',', $end_time_str);
	$year2 += 1900;
	my $end_time=$year2.'/'.$month2.'/'.$monthday2.', '.$hour2.':'.$min2.':'.$sec2;
	my $duration_m=($year2-$year1)*365*24*60 + ($yearday2-$yearday1)*24*60 + ($hour2-$hour1)*60 + ($min2-$min1);
	#
	if($type eq 'duration'){
		$my_time = $begin_time.'-----'.$end_time.'. Duration: '.$duration_m.'min' if $duration_m > 0;
	}
	elsif($type eq 'early'){
		$my_time= ($duration_m>=0) ? $begin_time_str : $end_time_str ;
	}
	elsif($type eq 'late'){
		$my_time= ($duration_m<=0) ? $begin_time_str : $end_time_str ;
	}
	elsif($type eq 'seconds'){
		$my_time = $duration_m*60 if $duration_m>0;
	}
	elsif($type eq 'readable'){
		$my_time = $hour1.':'.$min1.', '.($month1+1).'/'.$monthday1.'/'.$year1;
	}
	#
	return($my_time);
}

###################
sub format_directory{
	my $indir=$_[0];
	
	#make directory if it doesn't exists
	my @dir_levels=split('/', $indir);
	my $outdir=$dir_levels[0];
	for (my $i=1; $i<@dir_levels; $i++){
		$outdir=$outdir.'/'.$dir_levels[$i];
		#print "$outdir\n";
		mkdir($outdir, 0755) unless -d $outdir;
	}
	#add '/' at the end of directory
	$outdir .= '/' unless $outdir=~/\/$/;
	return($outdir);
}

###########################
#judge if two fragments are overlapped.
sub overlapping{
	my($mate_start, $mate_end, $start, $end)=@_;
	
	my $judging;
	if ($end<=$mate_start ){#2
		$judging='upstream';
	}#2
	elsif($start>=$mate_end){#2
		$judging='downstream';
	}#2
	else{#2
		if($start<$mate_start){#3
			$judging= ($end<$mate_end)? 'cross' : 'cover';
		}#3
		elsif($start==$mate_start){#3
			if ($end<$mate_end){	$judging='in';		}
			elsif ($end==$mate_end){	$judging='same';		}
			else{	$judging='cover';		}
		}#3
		else{#3
			$judging= ($end<=$mate_end)? 'in' : 'cross';
		}#3
	}#2
	
	return($judging);
}

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

	return(\%files_hash);  
}#1

########################
#at least 3 parameters, usually 4 or more
sub multi_threading{#1
	my $func=shift @_;  # the passed subFunction
	my $threads_num=shift @_; # number of multi-threads
	my $inputs_pointer=shift @_; # names of samples
	my @inputs=@$inputs_pointer;
	
	my $num=@inputs;
	if ($threads_num<2){
		foreach my $input(@inputs){
			printf( "Single processing: \t %s %s\n", $num, $input);
			$func->($input, @_);
			$num--;
		}
	}
	else{#2
		while(1){#3
			if(threads->list() < $threads_num and @inputs>0 ){#4
				my $input=shift @inputs;
				printf( "Parellel processing: \t %s %s\n", $num, $input);
				threads->create(\&$func, $input, @_); 
				$num--;
			}#4
			#recover all threads
			foreach my $sub_thread( threads->list() ){#4
				$sub_thread->join() if $sub_thread->is_joinable();
			}#=4
			last if threads->list()==0 and @inputs==0;
		}#3
	}#2
	#
}#1

##########################
#get file names with input url
sub web_list_files{
	my ($url, $out_type, $pattern)=@_;
	$url=~s/\/$//;
	#
	my @file_names;
	my @html=split(/\n/, LWP::Simple::get $url);
	foreach my $web_line(@html){
		my @items=split(/ +/, $web_line);
		push(@file_names, $items[-1]);
		#print "$items[-1]\n";
	}
	#
	if ($pattern){
		@file_names=grep {$_=~/$pattern/} @file_names;
	}
	#
	if ($out_type eq 'hash'){
		my %web_files;
		foreach (@file_names){
			$web_files{$_}=$url.'/'.$_;
		}
		return(\%web_files);
	}
	elsif ($out_type eq 'files'){
		my @files=map {$url.'/'.$_} @file_names;
		return(\@files);
	}
	elsif ($out_type eq 'file'){
		my @files=map {$url.'/'.$_} @file_names;
		my $first_file=shift @files;
		return($first_file);
	}
	else{# file_names
		return(\@file_names);
	}
	#
}

###################3
sub web_download_file{
	my($file_url, $out_dir, $unpack)=@_;
	$file_url=~s/\/$//;
	$out_dir=format_directory($out_dir);
	
	#download file
	my $ff = File::Fetch->new(uri => $file_url);
	$ff->fetch(to=>$out_dir) or die $ff->error;
	#return local file
	my $local_file=$out_dir.file_operation($file_url, 'file_name');
	#uncompress gz file
	if ($unpack==1){
		#system("rm $local_file") if -f $local_file;
		print "uncompress $local_file\n";
		system("gunzip $local_file");
		$local_file=~s/\.gz$//;
	}
	printf( "Download a file from %s and save it into %s\n", $file_url, $local_file);
	return($local_file);
}

##############
sub check_duplicates{
	my ($infile, $index, $sep)=@_;
	my %hash;
	open my($IN), "<", $infile or die;
	while(<$IN>){
		chomp($_);
		my @items=split($sep, $_);
		my $key=$items[$index];
		if($hash{$key}){
			printf("%s\n, %s\n\n", $hash{$key}, $_);
		}
		else{
			$hash{$key}=$_;
		}
		print "$_\n" unless @items==33;
	}
	close($IN);
}
########################
#select one or multiple separated by comma
sub stdin_select{
	my($pointer, $label)=@_;
	my @arr=@$pointer;
	$label //= "Select";
	printf("\n###%s###\n", $label);
	#
	my $out='';
	while($out eq ''){#2
		#print "$out\n";
		#display
		for(my$i=0; $i<@arr; $i++){
			printf ("\t[%d] %s\n", $i+1, $arr[$i]);
		}
		print"Select:";
		my $input=<STDIN>;
		chomp($input);
		if($input eq ''){
			$out=shift @arr;
			#print "$out\n";
		}else{
			my @index_arr=sub_data::array_unique( split(',', $input) );
			my @select_arr=map {$arr[($_-1)]} @index_arr;
			$out=join(',', @select_arr);
		}
		print "$out\n";
	}#2
	print "##########\n";
	return($out);
}


##########################3
sub stdin_dir{
	my($default, $label)=@_;
	$label //= "Enter directory:"; #default
	printf("###%s(default: %s):\n", $label, $default);
	#
	my $out='';
	while($out eq ''){#2
		#print "$out\n";
		#display
		my $input=<STDIN>;
		chomp($input);
		if($input eq ''){
			$out=$default;
		}elsif(-d $input){
			$input .= '/' unless $input=~/\/$/;
			$out=$input;
		}else{
			$input=Cwd::getcwd().'/'.$input unless $input=~/^\//;
			printf("%s doesn't exits, so create it.\n", $input);
			$out=format_directory($input);
		}
		print "$out\n";
	}#2
	print "\n";
	return($out);
}
##########################3
sub stdin_file{
	my($default, $label)=@_;
	$label //= "Enter file name with it path:"; #default
	printf("###%s:\n", $label);
	#
	my $out='';
	while($out eq ''){#2
		#print "$out\n";
		#display
		my $input=<STDIN>;
		chomp($input);
		$out=$input if -f $input;
		print "$out\n";
	}#2
	print "\n";
	return($out);
}
########################
sub stdin_DNA{
	my($default, $label)=@_;
	$label //= "Enter DNA sequence:"; #default
	printf("###%s:\n", $label);
	#
	my $out='';
	while($out eq ''){#2
		my $input=<STDIN>;
		chomp($input);
		if($input eq ''){
			$out=$default;
		}else{
			$out=$input;
		}
		print "$out\n";
	}#2
	print "\n";
	return($out);
}
###

#############################
1;  # make sure the file returns true or require will not succeed!#
