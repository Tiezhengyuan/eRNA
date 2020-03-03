use warnings;
use strict;
use Cwd;

use func_rna;
use func_common; #sub_common::
use func_basic; #sub_basic::
use func_bioseq;


my $a = "attcggcagcttatttctgt";
my $b = "actcggccatgtatcgctgt";
my $c = "gctcggccatgtatcgctgt";
#print sub_bioseq::shared_seq2($a,$b), "\n";
my @seq_arr=($a,$b, $c);
#print sub_bioseq::shared_seq(\@seq_arr), "\n";

#print sub_bioseq::seq_mismatch(\@seq_arr), "\n";

#my @cluster=qw/1,2;3 2,3;1 2,1,3 3;4 5/;
#print "@cluster\n";
#sub_bioseq::remove_same_comb(\@cluster);

#split fastq file
sub_bioseq::split_fastq_file('/home/yuan/raw_data/test_rawdata', '/home/yuan/raw_data/test_rawdata/t');


print "ok\n";
