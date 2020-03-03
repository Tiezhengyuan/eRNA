use warnings;
use strict;
use constant false => 0;
use constant true  => 1;
use Cwd;

use func_rna;
use func_common;    #sub_common::
use func_basic;     #sub_basic::
use func_gui;
use func_bioseq;
use func_data;

my @a1 = qw/5 4 10 5 6 0 7 8/;
my @a2 = qw/3 8cf 3 5fa 4fa 6fa 7ba 8cf 3 ' '' '/;
my %d1 = ( 1 => 2, 2 => 3 );
my %d2;
$d2{'a'}->{'stat'}=30;
$d2{'c'}->{'stat'}=1;
$d2{'d'}->{'stat'}=8;
$d2{'e'}->{'stat'}=-1;
#my %d3=%{$d2{'a'}};

my $s='tr|D6RHA5|D6RHA5_MOUSE|CTERM|STOP|C-PADDED-16';
if ($s=~/CTERM|STOP/){ print '4';}


