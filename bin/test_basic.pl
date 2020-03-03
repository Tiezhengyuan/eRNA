#! /usr/bin/perl -w
use strict;
use warnings;
use Bio::SeqIO;
use List::MoreUtils;
use List::Util;
use File::Find;


#################

require "func_basic.pm"; #sub_basic::


######
my %variables=('dir_raw_data'=>'/home/yuan/eRNA/raw_data/', 'file_sample_info'=>'/home/yuan/eRNA/result/sample_info.csv');
#sub_basic::SM_raw_files(\%variables);
sub_basic::SM_sample_info(\%variables);


print "ok\n";