#! /usr/bin/perl -w
use strict;
use warnings;
use threads; 
use Cwd;
use List::Util;
use List::MoreUtils;


################################################################################################################
#main programe
################################################################################################################
print "Main running begins:\n\n\n";
#get the directory of perl scripts involved in Pscore
my $perl_dir=Cwd::getcwd();

#get subroutines 
require $perl_dir."/eRNA_subroutines.pm";

#get E_RNA variables from the file E_RNA_pipeline.info
#it should be saved in the same directory with E_RNA.pl
my $var_file=$ARGV[0];
#my $var_file='/home/yuan/mysql_pre/eRNA/result/variables.txt';
my $variables_pointer=E_RNA::process_info($var_file);
my %variables=%$variables_pointer;
foreach(sort (keys %variables) ){
	print "$_=$variables{$_}\n";
}
my @sample_names=split(",", $variables{sample_names});

#initiate envirinmental variables
$ENV{'PATH'} .= ":$variables{alignment_dir}";
print $ENV{'PATH'}, "\n";


#initiate and clear  the monitor.log
open my($MON), ">", $variables{monitor_log_file} or die; 
#get *.fastq files list
my $files_pointer=E_RNA::files_list($variables{rawdata_dir}, 'incrusive_files', 'fastq');
my @total_fastq_files=@$files_pointer;
foreach my $sample_name(@sample_names){
	my @fastq_files;
	if(-f $variables{sample_info_file}){
		my $sample_info_pointer=E_RNA::read_sample_info($variables{sample_info_file}, 'hash');
		my %sample_info=%$sample_info_pointer;
		my @fastq_names=split(',', $sample_info{$sample_name}->{'fastq_names'});
		foreach my $name(@fastq_names){
			foreach my $file(@total_fastq_files){
				push(@fastq_files, $file) if $file=~/\/$name\_/;
			}
		}
	}
	else{	@fastq_files=grep(/$sample_name/, @total_fastq_files);	}
	#return the size of all fastq files
	my $file_size=0;
	foreach(@fastq_files){
		my @args=stat($_);
		$file_size +=$args[7];
	}
	my $supposed_time=int($file_size/228383);
	print $MON $sample_name, ":supposed_time=$supposed_time\n";
	print $MON $sample_name, ":beginning_time=NA\n";
	print $MON $sample_name, ":ending_time=NA\n";
	print $MON $sample_name, ":running_time=NA\n";
}
close($MON);


#get references information
$variables_pointer=E_RNA::references_info(\%variables);
%variables=%$variables_pointer;


#############
#main running

#DEG analysis
if ($variables{R_targets_method} eq 'DEseq' and $variables{R_DEG_analysis} eq 'yes'){
	print "\n\nDifferential expression profiling analysis using DEseq in R environment!\n\n";
	system("Rscript $perl_dir/R_DEseq.R > $variables{result_dir}/R.log");
}
if ($variables{R_targets_method} eq 'RF' and $variables{R_RF_analysis} eq 'yes'){
	print "\n\nTarget genes prediction using RandomForest in R environment!\n\n";
	system("Rscript $perl_dir/R_RF.R > $variables{result_dir}/R.log");
}



print "\n\nAnalysis is done! GOOD LUCK!\n\n\n";
















