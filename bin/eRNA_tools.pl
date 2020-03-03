#! /usr/bin/perl

use strict;
use warnings;
use List::Util;
use List::MoreUtils;

#get subroutines
require "func_data.pm";
require "func_gui.pm";
require "func_basic.pm";
require "func_common.pm";
require "func_rna.pm";
require "func_mrna.pm";
require "func_ncrna.pm";
require "func_bioseq.pm";

#get the directory of exe scripts involved in
my $dir_home=Cwd::abs_path(Cwd::getcwd().'/..' ).'/';


#first select operations
my @operations=('Split FASTQ files', 'Analyze protein *.gbK', 'Analyze RNA *.gbK to fa',);
my $operation=sub_common::stdin_select(\@operations, 'Select operations');


######################################
if ($operation eq 'Split FASTQ files'){
	my $in_dir=sub_common::stdin_dir($dir_home, 'Enter a directory of storing FASTQ files');
	my $out_dir=sub_common::stdin_dir($in_dir, 'Enter a directory for storing splited files');
	
	#split fastq file
	sub_bioseq::split_fastq_file($in_dir, $out_dir);
}

######################################
#extract annotation save into a table txt
if ($operation eq 'Analyze protein *.gbk'){
	my $gbk_file=sub_common::stdin_file("", 'Enter a GBK file');
	my $out_file=sub_common::file_operation($gbk_file, 'file_head').'.txt';
	sub_bioperl::protein_gbk($gbk_file, $out_file);
}

######################################
if ($operation eq 'Analyze RNA *.gbK to fa'){
	my $gbk_file=sub_common::stdin_file("", 'Enter a GBK file');
	my $result_dir=sub_common::file_operation($gbk_file, 'directory').'/RNA_annot';
	print "Enter specie or annotation name:";
	my $specie=<STDIN>;
	#
	sub_bioperl::gbk_to_fa($gbk_file, $result_dir, $specie);
}
#####################################
print "\n\nGreat! It is done!\n\n";