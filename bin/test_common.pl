#! /usr/bin/perl -w
use strict;
use warnings;
use List::Util;
use File::Find;


############
#import personal modules
use func_basic;
use func_common;
use func_bioperl;
use func_bioseq;
use func_data;


#############
#sub_common::check_duplicates('/home/yuan/phip/ref_seq/T7Pep2_human_annot.txt', 0, "\t");

my @a=qw/1 2 3 5/;
my @b=qw/5 1 2 /;
my @c=qw/6 7 8 1 5/;
#sub_common::array_diff(\@a,\@b);
#sub_common::array_fdiff(\@a,\@b);
#sub_common::array_share2(\@a,\@b);
#sub_common::array_share(\@a,\@b, \@c);
#sub_common::array_merge(\@a,\@b, \@c);


my @strings = qw(A B C D E F G H I J K L M N O P Q R S T U V W X Y Z);
#print "@strings\n";

my @add=qw/1 2 3/;
sub_common::array_to_file(\@strings, '/home/yuan/test.txt');

#my $perm_pointer=sub_common::array_permutation(\@strings, 6);
#my @perm=@$perm_pointer;
#print "@$_\n" for @perm;

@strings = qw(1 2 3 4 5 6);
#width first
#my $comb_pointer=sub_common::array_combination(\@strings, 4);
#depth first
#my $comb_pointer=sub_common::array_combination2(\@strings, 4);
#my @comb=@$comb_pointer;
#print "@$_\n" for @comb;

@a=qw/1 2 3 7/;
@b=qw/6 12 21/;
#print sub_common::array_identical(\@a, \@b);
#my $pointer=sub_common::array_fdiff(\@a, \@b);
#print @$pointer, "\n";

@a=qw/a,j,b c,d f/;
#print sub_common::unique_combination(\@a);

#@a=('a,b', 'z','a,b,c',  'c,d', 'e', 'e,h', 'h,j', 'x', 'a','b','c','d','e','h','j');
#sub_common::array_cluster(\@a, 2);

#std
#print sub_common::stdin_select(\@a, 'dsa');
#sub_common::stdin_dir('/home/yuan/', 'dsa');
print "OK\n";