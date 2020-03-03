#!/usr/bin/perl

use strict;
use LWP::Simple;
require "func_download.pm"; #sub_download::

=head1
print STDERR "\n
###############################################################################################
#
# This is the ### installer. 
# It will work under a bash,csh and ksh shell.
# It will try to download all necessary files and install them. 
# Please restart your shell to make changes take effect
#
###############################################################################################


";
=cut

###########################
#current directory
my $dir_bin= `pwd 2>&1`; 
chomp $dir_bin;
#printf("###The bin directory: \n", $dir_bin);
#home directory
my $dir=$dir_bin;
$dir=~s/bin$//;



#####################
#2
#check default software 
#sub_download::check_gcc();
#sub_download::check_wget();
#sub_download::check_gzip();

########################
#3
#sub_download::get_bowtie1($dir);
#sub_download::get_bowtie2($dir);
#sub_download::get_tophat($dir);
#sub_download::get_cufflinks($dir);
#sub_download::get_samtools($dir);

print "ok\n";
