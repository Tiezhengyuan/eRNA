#! /usr/bin/perl -w
use strict;
use warnings;
use List::MoreUtils;
use List::Util;
use File::Find;

####################################################
#
#the file constains all subroutines required for running E_RNA_pipeline
package sub_data;

##############################
#list unique elements given an array
sub array_unique(@){
	#equalvalent to List::MoreUtils::uniq
	#List::MoreUtils::uniq @_;
	
	#
	my %seen;
    my @uniq = grep {!$seen{$_}++} @_;
    #print @uniq;
    return(@uniq);
}

##########################
#find shared element between two arrays
sub array_share2{
	my ($arr1_pointer, $arr2_pointer) =@_;
	my @arr1=@$arr1_pointer;
	my @arr2=@$arr2_pointer;
	
	my @share;
	foreach my $a(@arr1){
		push(@share, $a) if List::Util::first {$a eq $_} @arr2;
	}
	#print "@share\n";
	return(\@share);
}
##################################
#multiple arrays
sub array_share{
	my $pointer=shift @_;
	my @share=@$pointer;
	while(@_>0){
		my $arr_pointer=shift @_;
		my @arr=@$arr_pointer;
		my $share_pointer=sub_common::array_share2(\@share,\@arr);
		@share=@$share_pointer;
		last if @share==0;
	}
	print "@share\n";
	return(\@share);
}
##########################
#find different elemets betwee two array
sub array_diff{
	my ($arr1_pointer, $arr2_pointer) =@_;
	my @arr1=@$arr1_pointer;
	my @arr2=@$arr2_pointer;
	
	my @diff;
	foreach my $a(@arr1){
		push(@diff, $a) unless List::Util::first {$a eq $_} @arr2;
	}
	foreach my $b(@arr2){
		push(@diff, $b) unless List::Util::first {$b eq $_} @arr1;
	}
	#print "@diff\n";
	return(\@diff);
}
#####################################
#judge difference
sub array_identical{
	my ($arr1_pointer, $arr2_pointer) =@_;
	my @arr1=@$arr1_pointer;
	my @arr2=@$arr2_pointer;
	
	my %h;
	foreach my $a(@arr1){
		if(exists $h{$a}){
			$h{$a} +=1;
		}
		else{
			$h{$a}=1;
		}
	}
	foreach my $b(@arr2){
		if(exists $h{$b}){
			$h{$b}+=1;
		}
		else{
			$h{$b}=1;
		}
	}
	my $tag='equal'; #default is unequal
	foreach (keys %h){
		if($h{$_}==1){
			$tag='unequal';
			last;
		}
	}
	return($tag);
}

#########################
sub unique_combination{
	my $pointer=shift @_;
	my @arr=@$pointer;
	
	my %h;
	my $n=0;
	foreach my $a_str(@arr){
		my @a=split(',', $a_str); #elements should be seperated by comma
		foreach my $b(@a){
			if (exists $h{$b}){
				return(0);
			}else{
				$h{$b}=1;
				$n++;
			}
		}
	}
	return($n);
}
#############################
#find the elemets in the first array that doesn't exist in the second one
sub array_fdiff{
	my ($arr1_pointer, $arr2_pointer) =@_;
	my @arr1=@$arr1_pointer;
	my @arr2=@$arr2_pointer;
	
	my @diff;
	foreach my $a(@arr1){
		push(@diff, $a) unless List::Util::first {$a eq $_} @arr2;
	}
	#print "@diff\n";
	return(\@diff);
}
##################################
#merge multiple arrays and removing duplicates
sub array_merge{
	my @merge;
	while(@_>0){
		my $arr_pointer=shift @_;
		my @arr=@$arr_pointer;
		foreach my $a(@arr){
			push(@merge, $a) unless List::Util::first {$a eq $_} @merge;
		}
	}
	#print "@merge\n";
	return(\@merge);
}
############################
#
sub array_permutation{
	my ($list_pointer, $n) = @_;
	my @list=@$list_pointer;
	my @cache;
	#outliers
	$n=@list if $n > @list;
	foreach (@list){
		my @a=($_,);
		push(@cache, \@a);
	}
	#
	my $tag=1;
	my @comb;
	while (@cache>0){#2
		my $pointer=shift @cache;
		my @a=@$pointer;
		if($n>@a){#3
			my $diff_pointer = array_fdiff(\@list, \@a);
			foreach (@$diff_pointer){#4
				my @new_a=@a;
				push(@new_a, $_);
				push(@cache, \@new_a);
				#print "@a\n";
			}#4
		}else{
			push(@comb, $pointer);
		}#3
	}#2
	return \@comb;
}

################
#width first
sub array_combination{
	my ($pool_pointer, $k) = @_;
	my @pool=@$pool_pointer;
	
	#init list of index
	my @arr; # take @arr as stack, kick off element at the head and add element at the end
	foreach my $index(0..(@pool-$k)){
		my @item_arr=($index);
		#print "@item_arr\t";
		push(@arr, \@item_arr);
	}
	
	#combination of index
	foreach my $tail_index( (@pool-$k+1)..(@pool-1) ){#2
		print "#$tail_index\t";
		my $len=@arr-1;
		foreach (0..$len){#3
			my $first=shift @arr; #kick off the first element
			#print @$first, ": ";
			my $head_index=@$first[-1]+1;
			foreach my $add_index($head_index..$tail_index){#4
				my @current=@$first;
				push(@current, $add_index);
				push(@arr, \@current);
				#print @current, "\t";
			}#4
			#print "\n";
		}#3
	}#2

	#translate index into elements
	my @comb;
	foreach my $pointer(@arr){
		my @a_values=map {$pool[$_]} @$pointer;
		#print "@a_values\n";
		push(@comb, \@a_values);
	}	
	return(\@comb);
}

################
#depth first
sub array_combination2{
	my ($pool_pointer, $k) = @_;
	my @pool=@$pool_pointer;
	
	#init list of index
	my @arr; # take @arr as stack, kick off element at the head and add element at the end
	foreach my $index(0..(@pool-$k)){
		my @item_arr=($index);
		#print "@item_arr\t";
		push(@arr, \@item_arr);
	}
	#print @arr, "\n";
	return(\@arr) if $k<=1; #outlier
	
	#combination of index
	my %tail_index_hash;
	my $current_index=1;
	foreach my $tail_index( (@pool-$k+1)..(@pool-1) ){#2
		$tail_index_hash{$current_index}=$tail_index;
		#print $current_index, "\t", $tail_index, "==\n";
		$current_index++;
	}
	
	#
	my @comb;
	while( @arr>0 ){#2
		my $first=pop @arr; #kick off the last element
		#print @$first, ": ";
		my $head_index=@$first[-1]+1;
		my $hash_index=@$first;
		#print $hash_index, "==\n";
		my $tail_index=$tail_index_hash{$hash_index};
		#printf("%s..%s\n", $head_index, $tail_index);
		foreach my $add_index($head_index..$tail_index){#4
			my @current=@$first;
			push(@current, $add_index);
			#print @current, "\t";
			if (@current<$k) {
				push(@arr, \@current); #add unfinished element
			}else{
				my @values=map {$pool[$_]} @current;
				#print @values, "\n";
				push(@comb, \@values);
			}
		}#4
		#print "\n";
	}#2

	return(\@comb);
}

#######################
#all combinations
sub array_all_combination{
	my ($list_pointer) = @_;
	my @list=@$list_pointer;
	
	my %comb_hash;
	my $max_num=@list;
	foreach my $n(1..$max_num){
		print "$n\n";
		my $comb_pointer=array_combination($list_pointer, $n);
		$comb_hash{$n}=$comb_pointer;
	}
	return(\%comb_hash);
}
######################
#cluster elements
sub array_cluster{
	my ($list_pointer, $cluster_min) = @_;
	my @list=@$list_pointer;
	my $total=@list;
		
	my @cluster;
	while(@list>0){#2
		my $l=@list;
		#print "$l\n";
		my $tag=0;
		my $first_str=pop @list;
		my @first_arr=split(/,|;/, $first_str); 
		@first_arr=sub_data::unique(@first_arr);
		foreach my $a_str(@list){#3
			my @a_arr=split(/,/, $a_str); 
			my $shared_pointer=array_share2(\@first_arr, \@a_arr);
			if (@$shared_pointer>0){
				@list= grep {$_ cmp $a_str} @list;
				#print "$first_str=$a_str\n";
				push(@list, $first_str.';'.$a_str);
				$tag=1;
				last;
			}
		}#3
		if ($tag==0){
			push(@cluster, $first_str) ;
			#print "$first_str\n";
		}
	}#2
	
	#combine small blocks
	my (@new, @other);
	foreach my $ss(@cluster){
		my @aa=split(';', $ss);
		if (@aa<=$cluster_min){ 	push(@other, $ss);	}
		else{	push(@new, $ss);	}
	}
	
	#export
	#my $n=1;
	#foreach my $s(@new){
	#	print "#$total:$n##\n";
	#	print "$_\n" for split(';', $s);
	#	$n++;
	#}
	return(\@new, \@other);
}

###########
#@arr=('a,b,c' 'b,c', 'd') => @arr=(a b c d)
sub array_flat{
	my($pointer, $sep)=@_;
	
	my @flat;
	foreach (@$pointer){
		my @a=split(/$sep/, $_);
		my $flat_pointer=array_merge(\@flat, \@a);
		@flat=@$flat_pointer;
	}
	return(\@flat);
}

##intersection of two arrays
sub two_intersection{#1
  	my ($a_pointer, $b_pointer)=@_;
	my @a=@$a_pointer;
	my @b=@$b_pointer;

	my %out=(a_num=>0, b_num=>0, a_only_num=>0, b_only_num=>0, ab_only_num=>0,  total_num=>0); 
	$out{a_num}=@a;
	$out{b_num}=@b;
	my @ab=(@a, @b);
	@ab = sub_data::unique( grep($_, @ab) );
	$out{total_num}=@ab;
	foreach my $item(@ab) {#3
		if (List::Util::first {$item eq $_} @a and List::MoreUtils::none {$item eq $_} @b){ $out{a_only_num}++; }
		if (List::Util::first {$item eq $_} @a and List::Util::first {$item eq $_} @b){ $out{ab_only_num}++; }
		else{  $out{b_only_num}++; }
	}#3
	
	return(\%out);
}#1
#############
#intersection of three arrays
sub three_intersection{#1
	my ($a_pointer, $b_pointer, $c_pointer)=@_;
	my @a=@$a_pointer;
	my @b=@$b_pointer;
	my @c=@$c_pointer;
	my %out=(a_num=>0, b_num=>0, c_num=>0, a_only_num=>0, b_only_num=>0, c_only_num=>0, 
           ab_only_num=>0, ac_only_num=>0, bc_only_num=>0, abc_only_num=>0, total_num=>0); 

	$out{a_num}=@a;
	$out{b_num}=@b;
	$out{c_num}=@c;
	my @abc=(@a, @b, @c);
	@abc = sub_data::unique(grep($_, @abc));
	$out{total_num}=@abc;
	foreach my $item(@abc) {#3
		if (List::Util::first {$item eq $_} @a) {#4
			if (List::Util::first {$item eq $_} @b and List::Util::first {$item eq $_} @c){
				$out{abc_only_num}++;
			}
			elsif (List::Util::first {$item eq $_} @b and List::MoreUtils::none {$item eq $_} @c){
				$out{ab_only_num}++;
			}
			elsif (List::Util::first {$item eq $_} @c and List::MoreUtils::none {$item eq $_} @b){
				$out{ac_only_num}++;
			}
			else{  $out{a_only_num}++; }
		}#4
		elsif (List::Util::first {$item eq $_} @b) {#4
			if ( List::Util::first {$item eq $_} @c){	$out{bc_only_num}++;	}
			else{  $out{b_only_num}++; }
		}
		else{  
			$out{c_only_num}++; 
			#print "$item\n";
		}
	}#3

	return(\%out);
}#1
###################################################
#the first two arrays
sub intersect{
		my ($arr1_pointer, $arr2_pointer)=@_;
		my @arr1=@$arr1_pointer;
		my @arr2=@$arr2_pointer;
		#print "@arr1\n @arr2\n\n" if @arr1==0 or @arr2==0;
		
		my @intersect;
		foreach my $a(@arr1){
			push(@intersect, $a) if List::Util::first {$a eq $_} @arr2;
		}
		return(\@intersect);
}
#intersections of arrays
sub intersect_arrays{
	#
	my @intersect;
	if (@_<2) {#2
		print "Error input\n";
	}#2
	else{#2
		my $arr1_pointer=shift @_;
		my $arr2_pointer=shift @_;
		#the first two arrays
		my $intersect_pointer=intersect($arr1_pointer, $arr2_pointer);
		#the other arrays
		while(@_>0){
			my $arr3_pointer=shift @_;
			$intersect_pointer=intersect($intersect_pointer, $arr3_pointer);
		}
		@intersect=@$intersect_pointer;
		#
		#print "==@intersect==\n";
	}#2
	return(\@intersect);
}

##################################
#############
sub combine_hash2_2{
		my($pointer1, $pointer2)=@_;
		my %hash1=%$pointer1;
		my %hash2=%$pointer2;
		
		#add hash2 into hash1 by pro_id
		foreach my $pro_id(keys %hash1){#2
				my $pointer1=$hash1{$pro_id};
				my %subhash1=%$pointer1;
				if (exists $hash2{$pro_id}){#3
						my $subpointer=$hash2{$pro_id};
						my %subhash2=%$subpointer;
						foreach my $key2(keys %subhash2){#4
								foreach my $pep_id(keys %subhash1){
										$hash1{$pro_id}->{$pep_id}->{$key2}=$subhash2{$key2};
								}
						}#4
				}#3
				else{#3
						print "$pro_id: no keys of the first hash2\n";
				}#3
		}#2
		return(\%hash1);
}

###################
sub combine_hash2{
	my($h1_pointer, $h2_pointer)=@_;
	my %h1=%$h1_pointer;
	my %h2=%$h2_pointer;
	
	foreach my $row(keys %h2){#2
		my $hash_pointer=$h2{$row};
		my %hash=%$hash_pointer;
		foreach my $col(keys %hash){
			if(exists $h1{$row}->{$col}){
				$h1{$row}->{$col} += $h2{$row}->{$col};
			}
			else{
				$h1{$row}->{$col}=$h2{$row}->{$col};
			}
		}
	}#2
	
	return(\%h1);
}

#######################################################################################
#distance of two strings for comparison
#algorithm is Needleman-Wunsch
sub NW_dist{#2
  my $str1 = $_[0];
  my $str2 = $_[1];
  my $len1 = length($str1);
  my $len2 = length($str2);

  # We have to be a little careful with the index of the table we are using
  # Since we do not want to reconstruct the alignment, only two rows of the
  # matrix are sufficient. Obviously we will use pointers to these rows, so
  # we can exchange them freely without copying elements.

  # Initialize old_row as reference
  my $old_row;
  for (my $i = 0; $i <= $len1; $i++){#2
    $old_row->[$i] = $i;
  }#2

  my $new_row;
  for (my $i = 1; $i <= $len2; $i++) {#2
    $new_row = [];
    $new_row->[0] = $i;
    for (my $j = 1; $j <= $len1; $j++) {
      $new_row->[$j]=$old_row->[$j]+1;
      my $b=$new_row->[$j-1]+1;
      my $c=$old_row->[$j-1] + ((substr($str1, $j-1, 1) eq substr($str2, $i-1, 1)) ? 0 : 1);
      $new_row->[$j]=$b if $new_row->[$j] > $b;
      $new_row->[$j]=$c if $new_row->[$j] > $c;
    }
    $old_row = $new_row;
  }#2
  return($new_row->[$len1]);
}#2


#######################
sub print_hash{
	my $pointer=shift @_;
	my %hash=%$pointer;
	
	my $n=1;
	foreach my $key (sort { $hash{$a} cmp $hash{$b} }  keys %hash){
		printf("%s: %s\t%s\n", $n, $key, $hash{$key});
		$n++;
	}
}

#########################
1;