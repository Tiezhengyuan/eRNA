#! /usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use threads; 
use List::Util;
use List::MoreUtils;


################################################################################################################
#main programe
################################################################################################################



#subroutines:
require $perl_dir."/functions_ncrna.pm"; #sub_ncrna::
require $perl_dir."/functions_rna.pm"; #sub_rna::
require $perl_dir."/functions_basic.pm"; #sub_basic::
require $perl_dir."/functions_common.pm"; #sub_common::

#get the directory of perl scripts involved in Pscore
our $perl_dir=Cwd::abs_path(Cwd::getcwd().'/..' );

#get E_RNA variables from the file E_RNA_pipeline.info
#it should be saved in the same directory with E_RNA.pl
my $var_file=$ARGV[0];
my $variables_pointer=sub_common::file_to_hash($var_file, '=');
my %variables=%$variables_pointer;
$variables{'file_var'}=$var_file;

#check directories
$variables_pointer=sub_rna::initiate_directories(\%variables, $perl_dir);
%variables=%$variables_pointer;
sub_common::print_hash(\%variables);#print %variables


#initiate envirinmental variables
$ENV{'PATH'} .= ":$variables{dir_bowtie1}";
print $ENV{'PATH'}, "\n";

#associate sample name with *.fastq files list
$variables_pointer=sub_basic::SM_sample_info(\%variables);
%variables=%$variables_pointer;

#initiate and clear  the monitor.log and system_monitor.log
my @total_beginning_time=localtime(time); #beginning time record
sub_basic::initiate_log_files(\%variables);

#get references information
print "\n\nGet references information.\n";
my $ref_info_pointer=sub_basic::pre_references_info(\%variables);
$variables{$ref_info_pointer}=$ref_info_pointer;

print "\n\ncheck bowtie index.\n";
sub_basic::check_bowtie_index(\%variables, 'bowtie1');



####################################################
#main running
print "\n\n\n###########################\n";
print "Main running begins:\n\n\n";
if($variables{ncRNA_read_rawdata} eq 'yes' or $variables{ncRNA_adapter_removal} eq 'yes' or 
		$variables{ncRNA_seperate_alignment} eq 'yes' or $variables{ncRNA_iterative_alignment} eq 'yes'){#1
	my @sample_names=split(',', $variables{sample_names});
	#print "$variables{sample_names}\n";
	while(1){#2
		my $threads_num=sub_basic::read_log($var_file, 'threads_num');
		$variables{threads_num}=$threads_num unless length($threads_num)==0 or $threads_num=~/[^0-9]/;
		if(threads->list() < $variables{threads_num} and @sample_names>0){#3
			my $sample_name=shift @sample_names;
			threads->create(\&sub_ncrna::main_running, \%variables, $sample_name);
		}#3
		#recover all threads
		foreach my $sub_thread( threads->list() ){#3
			if ( $sub_thread->is_joinable() ){
				$sub_thread->join() ;
			}
		}#3
		last if threads->list()==0 and @sample_names==0;
		sleep 1;
	}#2
	
}#1

print "\n\nMain running is done!\n\n\n";

#generate data frame based on alignment
if ($variables{ncRNA_counting} eq 'yes'){
	print  "\n\n Export statistical results and data frames of transcriptional level\n";
	sub_ncrna::ncRNA_S5_counting(\%variables) ;
}

my @total_ending_time=localtime(time);
sub_basic::refresh_log($variables{file_time_monitor_log}, "Total:ending_time" , join(",", @total_ending_time) );
my $total_running_time=sub_common::get_time(join(',', @total_beginning_time), join(',',@total_ending_time) );
sub_basic::refresh_log($variables{file_time_monitor_log}, "Total:running_time" , $total_running_time );

print "\n\nAnalysis is done! GOOD LUCK!\n\n\n";
















