use strict;
use warnings;
use Bio::Perl;
use Bio::SeqIO;
use Bio::DB::GenBank;
use Bio::DB::Query::GenBank;
use Bio::DB::SwissProt;
use Bio::DB::EntrezGene;
use Bio::DB::Taxonomy;
#use IO::String;
use Bio::SeqFeatureI;
use threads;
use List::MoreUtils;
use List::Util;

############
#import personal modules
use lib "/home/yuan/eRNA/bin";
use func_basic;
use func_common;
use func_data;
use func_bioperl;
use func_bioseq;
#require "/home/yuan/eRNA/bin/func_common.pm"; #sub_common::
#require "/home/yuan/eRNA/bin/func_basic.pm"; #sub_basic::
#require "/home/yuan/eRNA/bin/func_bioperl.pm"; #sub_bioperl::

####################################
package sub_phip;


################
sub C_protein_motifs{
	my($attr_pointer)=@_;
	my %attr=%$attr_pointer;
	my %hash=%{$attr{annot_pointer}};
	my $out_dir=sub_common::format_directory($attr{annot_dir}.$attr{lib}.'_ProMotifs/');
	
	#get motif info
	my $motif_file=$attr{annot_dir}.'prosite_protein_motifs.txt';
	sub_bioperl::PROSITE_protein_motifs($motif_file) unless -f $motif_file;
	my $motifs_pointer=sub_common::file_to_hash2($motif_file, "\t");
	my %motifs=%$motifs_pointer;
	
	#
	my $n=keys %hash;
	foreach my $pro_id(keys %hash){#2
		printf( "%s:%s\n", $n, $pro_id);
		my %pep_motif;
		my $pointer2=$hash{$pro_id};
		my %h2=%$pointer2;
		foreach my $pep_id(keys %h2){#3
			my $key='pro_motifs:'.$pep_id;
			$pep_motif{$key}='none';
			my $pep_aa=$h2{$pep_id}->{'pep_aa'};
			foreach my $motif_acc(keys %motifs){
				my $pattern=$motifs{$motif_acc}->{'regular_pattern'};
				#print "$pattern\n";
				if ($pep_aa=~/$pattern/){
					#my $description=$motifs{$motif_acc}->{'description'};
					my $value=$motif_acc.','.$&;
					$pep_motif{$key}= ($pep_motif{$key} eq 'none') ? $value : $pep_motif{$key}.';'.$value;
				}
			}
			#print "$pep_motif{$pep_id}\n";
		}#3
		#export
		sub_common::hash_to_file(\%pep_motif, $out_dir.$pro_id, '=');
		$n--;
	}#2
}
#####################################################
#hash{pro_id}->{pep_id}->{XX}=XX
sub D_export_fa{
		my($pointer, $pro_cols_pointer, $file_head, $pros_pointer)=@_;
		my %h=%$pointer;
		my @pro_cols=@$pro_cols_pointer;
		my $fa_file=$file_head.'.fa'; #DNA fasta file
		my $AA_file=$file_head.'_AA.fa'; #AA fasta file
		my $annot_file=$file_head.'_annot.txt'; #annotation txt file
		my @pros= ($pros_pointer) ? @$pros_pointer : sort (keys %h);
		
		#export into text
		my $pro=0;
		my $pep=0;
		open my($OUT1), ">", $fa_file or die;
		open my($OUT2), ">", $annot_file or die;
		open my($OUT3), ">", $AA_file or die;
		printf ("\nExport annotations into the file %s, %s and %s\n\n: ", $fa_file, $AA_file, $annot_file);
		print $OUT2 join("\t", 'pep_id', 'pro_id', @pro_cols), "\n";
		foreach my $pro_id ( @pros){#2
				my %subh=%{$h{$pro_id}};
				foreach my $pep_id(sort (keys %subh)){ #3
						#export to fasta file
						print $OUT1 ">$pep_id\n", "$h{$pro_id}->{$pep_id}->{phip_seq}\n";
						print $OUT3 ">$pep_id\n", "$h{$pro_id}->{$pep_id}->{pep_aa}\n";
						#print ">$pep_id\n", "$h{$pro_id}->{$pep_id}->{taxon_phip}\n";
						#export to annotation file
						my @a=($pep_id, $pro_id);
						#print "$h{$key}->{pep_aa}\n";
						for(@pro_cols){
							my $v= (exists $h{$pro_id}->{$pep_id}->{$_}) ? $h{$pro_id}->{$pep_id}->{$_} : 'NA';
							push(@a, $v);
							#print "##$_:$h{$pro_id}->{$pep_id}->{$_}==$v==\n" if $pep_id=='9540';
						}
						print $OUT2 join("\t", @a), "\n";
						$pep++;
						#printf("\n%s=%s=%s=\n", $pep, $pro_id, $pep_id) if $pro_id eq 'P08318';
				}#3
				$pro++;
		}#2
		close($OUT1);
		close($OUT2);
		close($OUT3);
		printf("%s proteins and %s peptides in fa.\n", $pro, $pep);
}

#hash{pro_id}->{pep_id}->{XX}=XX
sub D_update_fa{
		my($pointer, $pro_cols_pointer, $file_head)=@_;
		my %h=%$pointer;
		my @pro_cols=@$pro_cols_pointer;
		my $fa_file=$file_head.'.fa'; #DNA fasta file
		my $AA_file=$file_head.'_AA.fa'; #AA fasta file
		my $annot_file=$file_head.'_annot.txt'; #annotation txt file
		
		#export into text
		my $pro=0;
		my $pep=0;
		open my($OUT1), ">>", $fa_file or die;
		open my($OUT2), ">>", $annot_file or die;
		open my($OUT3), ">>", $AA_file or die;
		printf ("\nAppend annotations into the file %s, %s and %s\n\n: ", $fa_file, $AA_file, $annot_file);
		foreach my $pro_id (sort (keys %h) ){#2
				my $pointer=$h{$pro_id};
				my %subh=%$pointer;
				foreach my $pep_id(sort (keys %subh)){ #3
						#export to fasta file
						print $OUT1 ">$pep_id\n", "$h{$pro_id}->{$pep_id}->{phip_seq}\n";
						print $OUT3 ">$pep_id\n", "$h{$pro_id}->{$pep_id}->{pep_aa}\n";
						#export to annotation file
						my @a=($pep_id, $pro_id);
						#print "$h{$key}->{pep_aa}\n";
						for(@pro_cols){
							my $v= exists $h{$pro_id}->{$pep_id}->{$_} ? $h{$pro_id}->{$pep_id}->{$_} : 'NA';
							push(@a, $v);
							#print "##$_:$h{$pro_id}->{$pep_id}->{$_}==$v==\n" if $pep_id=='9540';
						}
						print $OUT2 join("\t", @a), "\n";
						$pep++;
				}#3
				$pro++;
		}#2
		close($OUT1);
		close($OUT2);
		close($OUT3);
		printf("%s proteins and %s peptides in fa.\n", $pro, $pep);
}

##############################################
####
#hash{pro_id}->{pep_id}->{XX}=XX
sub Ztrash_E_add_human{
	my($infile, $pro_cols_pointer, $file_head)=@_;
	
	#read infile
	my %h;
	open my ($IN), "<", $infile or die;
	my $first=<$IN>;
	while ( my $L=<$IN>) {#2
		chomp($L);
		my($order, $pep_id,  $pep_dna, $pep_aa, $UniProt_acc, $pep_position, $gene_symbol, 
		$product, $gene_synonym, $pro_len, $UniProt_keywords, $UniProt_features)=split("\t", $L);
		$UniProt_acc=~s/^UniRef100_//;
		$h{$UniProt_acc}->{$pep_id}->{'UniProt_acc'}=$UniProt_acc ;
		$h{$UniProt_acc}->{$pep_id}->{'pro_id'}=$UniProt_acc ;
		$h{$UniProt_acc}->{$pep_id}->{'phip_seq'}=$pep_dna;
		$h{$UniProt_acc}->{$pep_id}->{'pep_dna'}=$pep_dna;
		$h{$UniProt_acc}->{$pep_id}->{'pep_aa'}=$pep_aa;
		$h{$UniProt_acc}->{$pep_id}->{'pep_position'}=$pep_position; 
		my($start,$end)=split('-', $pep_position); 
		$h{$UniProt_acc}->{$pep_id}->{'pep_rank'}=$start;
		$h{$UniProt_acc}->{$pep_id}->{'gene_symbol'}=$gene_symbol;
		$h{$UniProt_acc}->{$pep_id}->{'product'}=$product;
		$h{$UniProt_acc}->{$pep_id}->{'gene_synonym'}=$gene_synonym, 
		$h{$UniProt_acc}->{$pep_id}->{'pro_len'}=$pro_len;
		#printf( "##%s##\n", $L);
		#printf( "##%s##\n", $L) if length($h{$UniProt_acc}->{$pep_id}->{'pro_id'})<4;
		#last if length($h{$UniProt_acc}->{$pep_id}->{'pro_id'})<4;
	}
	close($IN);
	
	#export into text
	D_update_fa(\%h, $pro_cols_pointer, $file_head);
}

###############
sub extract_control_pep{
	my ($infile)=@_;
	
	my %control_pep;
	open my ($IN), "<", $infile or die;
	my($L1, $seq);
	while ( defined ($L1=<$IN>) && defined ($seq=<$IN>)) {#2
		chomp($L1, $seq);
		if($L1 eq '>Rnl2_SPIKEIN_'){#3
			$seq=~tr/a-z/A-Z/;
			my $fseq=substr($seq, 0, 51);
			my $aa_seq=sub_bioseq::translate_DNA($seq);
			$control_pep{'Rnl2_SPIKEIN'}->{'Rnl2_SPIKEIN'}->{pep_id}='Rnl2_SPIKEIN';
			$control_pep{'Rnl2_SPIKEIN'}->{'Rnl2_SPIKEIN'}->{pro_id}='Rnl2_SPIKEIN';
			$control_pep{'Rnl2_SPIKEIN'}->{'Rnl2_SPIKEIN'}->{phip_seq}=$fseq;
			$control_pep{'Rnl2_SPIKEIN'}->{'Rnl2_SPIKEIN'}->{pep_aa}=$aa_seq;
			$control_pep{'Rnl2_SPIKEIN'}->{'Rnl2_SPIKEIN'}->{pep_dna}=$seq;
		}#3
	}#2
	close ($IN);
	return(\%control_pep);
	
}

##################################
sub export_annot_table{
	my ($pointer, $indir, $outdir)=@_;
	my @colnames=@$pointer;
	my %taxon;
	
	#get taxon str
	my $pro_ids_pointer=sub_common::files_list($indir, 'file_name');
	my @pro_ids=@$pro_ids_pointer;
	for my $pro_id(@pro_ids){#2
		#get annotions 
		my $pointer=sub_common::file_to_hash($indir.$pro_id, "=");
		my %annot=%$pointer;
		my @tmp_taxon;
		for my $col(@colnames){#3
			if (exists $annot{$col}){ 
				push(@tmp_taxon, $annot{$col});
			}else{
				push(@tmp_taxon, 'NA');
			}
		}#3
		my $taxon_str=join(',', @tmp_taxon);
		$taxon{$pro_id} = $taxon_str;
	}#2
	
	#export
	sub_common::hash_to_file(\%taxon, $outdir.'taxonmy_list.txt', "\t");
}
###############
sub A_extract_human{
	my ($attr_pointer)=@_;
	my %attr=%$attr_pointer;
	my $pep=0;
	my (%annot, %pep_pro);
	
	#extract from infile1
	open my ($IN1), "<", $attr{phip_fa} or die;
	my($L1, $seq);
	while ( defined ($L1=<$IN1>) && defined ($seq=<$IN1>)) {#2
		chomp($L1, $seq);
		#print "$L1\n$seq\n";
		if ($L1=~/^>hb/){#3 human only
			#sequence
			$seq=~tr/a-z/A-Z/;
			my $fseq=substr($seq, 0, 51);
			my $pep_aa=sub_bioseq::translate_DNA($seq);
			#print length($fseq), $fseq, "\n";

			#name
			my ($pep_id, $pro_id, $pep_rank,$pep_rank_type, $pro_symb)=split(/\|/, $L1);
			$pep_id=~s/^>|"//g;
			unless (exists $pep_pro{$pep_id}){#4 remove duplicates in fasta
				$annot{$pro_id}->{$pep_id}->{name}=$L1;
				$annot{$pro_id}->{$pep_id}->{phip_seq}=$fseq;
				#print "$pro_id\n";
				$annot{$pro_id}->{$pep_id}->{pep_rank}=$pep_rank;
				$annot{$pro_id}->{$pep_id}->{pep_dna}=$seq;
				$annot{$pro_id}->{$pep_id}->{pep_aa}=$pep_aa;
				$pep_pro{$pep_id}=$pro_id;
				$pep++;
			}#4
		}#3
		else{#3
			print "$L1\n" unless $L1=~/^>[0-9]|^">[0-9]/;
		}#3
	}#2
	close ($IN1);
	my @pros=sort(keys %annot);
	
	#extract from the infile2
	my @pros_add;
	open my ($IN2), "<", $attr{phip_human1} or die;
	my $first=<$IN2>;
	while ( my $L=<$IN2>) {#2
		chomp($L);
		my($order, $pep_id,  $pep_dna, $pep_aa, $UniProt_acc, $pep_position, $gene_symbol, 
		$product, $gene_synonym, $pro_len, $UniProt_keywords, $UniProt_features)=split("\t", $L);
		$UniProt_acc=~s/^UniRef100_//;
		$annot{$UniProt_acc}->{$pep_id}->{'phip_seq'}=$pep_dna;
		$annot{$UniProt_acc}->{$pep_id}->{'pep_dna'}=$pep_dna;
		$pep_aa=~s/\*$//;
		$annot{$UniProt_acc}->{$pep_id}->{'pep_aa'}=$pep_aa;
		$annot{$UniProt_acc}->{$pep_id}->{'pep_position'}=$pep_position; 
		my($start,$end)=split('-', $pep_position); 
		$annot{$UniProt_acc}->{$pep_id}->{'pep_rank'}=$start;
		$annot{$UniProt_acc}->{$pep_id}->{'gene_symbol'}=$gene_symbol;
		$annot{$UniProt_acc}->{$pep_id}->{'product'}=$product;
		$annot{$UniProt_acc}->{$pep_id}->{'gene_synonym'}=$gene_synonym, 
		$annot{$UniProt_acc}->{$pep_id}->{'pro_len'}=$pro_len;
		#printf( "##%s##\n", $L);
		#printf( "##%s##\n", $L) if length($annot{$UniProt_acc}->{$pep_id}->{'pro_id'})<4;
		#last if length($h{$UniProt_acc}->{$pep_id}->{'pro_id'})<4;
		#aa stretches
		$annot{$UniProt_acc}->{$pep_id}->{aa_arr_pointer}=aa_overlapped($pep_aa);
		$pep_pro{$pep_id}=$UniProt_acc;
		$pep++;
		push(@pros_add, $UniProt_acc) unless List::Util::first {$_ eq $UniProt_acc} @pros_add;
	}
	close($IN2);
	push(@pros, sort @pros_add); #append additional protein ids
	
	#export
	my $pro_num=@pros;
	printf("Proteins: %s, peptides: %s\n", $pro_num, $pep);
	$attr{annot_pointer}=\%annot;
	$attr{pep_pro_pointer}=\%pep_pro;
	$attr{pros_pointer}=\@pros;
	return(\%attr);
}


########################
sub Ztrash_humanS2_SwissProt_query{#1
	my ($hash_pointer, $annot_dir, $type)=@_;
	my %hash=%$hash_pointer;
		
	#print "###add uniprot accession number first\n";
	print "read ref_uniprotkb file\n";
	my $refseq_uniprot_file='/home/yuan/data_preparation/ftp.ncbi.nlm.nih.gov/refseq/uniprotkb/gene_refseq_uniprotkb_collab';
	my $pu_pointer=sub_common::file_to_hash($refseq_uniprot_file, "\t");
	my %pu=%$pu_pointer;
	
	#
	my (@uniprot_accs, %nu);
	my $GB_dir=$annot_dir.'human_GenBank/';
	foreach my $acc(keys %hash){#2
		my $acc_file=$GB_dir.$acc;
		if (-f $acc_file){#3
			my $pro_id=sub_basic::read_log($acc_file, 'protein_id');
			$pro_id=~s/\.[0-9]*//;
			if (exists $pu{$pro_id}){
				my @a=split(',', $pu{$pro_id});
				@uniprot_accs=(@uniprot_accs, @a);
				$nu{$acc}->{first}=$a[0];
				$nu{$acc}->{all}=$pu{$pro_id};
				#printf("GenBank gene:%s \t protein:%s\t SwissProt %s\n", $acc, $pro_id, $pu{$pro_id} );
			}
			else{
				#print "############no swissprot annotation for gene $acc\n";
			}
		}#3
	}#2
	
	#parrallel processing
	#threads num: 8
	my $tmp_dir=$annot_dir.'human_annot_SwissProt_tmp/';
	mkdir($tmp_dir, 0755) unless -d $tmp_dir;
	#query accession
	my $query_accs_pointer=judge_query_acc(\@uniprot_accs, $tmp_dir, $type);
	my @query_accs=@$query_accs_pointer;
	if (@query_accs>0){
		my $num=@query_accs;
		printf( "\nThe total numer of proteins for SwissProt searching is %s\n\n", $num);
		#
		sub_common::multi_threading(\&sub_bioperl::SwissProt_query, 8, \@query_accs, $tmp_dir);
	}
	#
	print "update annotation by GenBank accession\n";
	my $out_dir=$annot_dir.'human_annot_SwissProt/';
	mkdir($out_dir, 0755) unless -d $out_dir;
	foreach my $acc(keys %nu){#2
		my $first_uniprot_acc=$nu{$acc}->{first};
		my $pointer=sub_common::file_to_hash($tmp_dir.$first_uniprot_acc, '=');
		my %uniprot_hash=%$pointer;
		$uniprot_hash{'uniprot_accs'}=$nu{$acc}->{all};
		sub_common::hash_to_file(\%uniprot_hash, $out_dir.$acc, '=');
	}

}#1

##################
sub human_RefGenome_query{
	my ($accs_pointer, $out_dir, $gff_file)=@_;
	my @accs=@$accs_pointer;
	
	my $n=0;
	my @out_features=('pos_start','pos_end', 'HGNC', 'MIM');
	open my($IN), "<", $gff_file or die;
	while (<$IN>){#2
		chomp($_);
		my @items=split("\t", $_);
		if ($_!~/^#/ and $items[2] eq 'mRNA'){#3 only mRNA line
			#get annotations of mRNA
			my %features;
			$features{'pos_start'}=$items[3];
			$features{'pos_end'}=$items[4];
			my @item_9=split(';', $items[8]);
			foreach my $fea(@item_9){#4
				my ($a, $b)=split('=', $fea);
				if ($a eq 'Dbxref'){
						my @sub_features=split(',', $b);
						foreach my $sub_fea(@sub_features){
							my @aa=split(':', $sub_fea);
							my $f_name=shift @aa;
							my $f_value=join(':', @aa);
							$features{$f_name}=$f_value;
						}
				}
				else{
					$features{$a}=$b;
				}
			}#4
			#
			$features{Genbank}=~s/\.[0-9]*//;
			if (List::Util::first {$features{Genbank} eq $_} @accs){#4
				$n++;
				@accs=grep {$_ cmp $features{Genbank}} @accs;
				printf("%d:%s\n", $n, $features{Genbank}); 
				my $acc_file=$out_dir.$features{Genbank};
				sub_common::hash_to_file(\%features, $acc_file, '=');
			}#4
		}#3
	}#2
	close($IN);
}

##################
sub human_PPI_query{
	my ($attr_pointer)=@_;
	my %attr=%$attr_pointer;
	my %human=%{$attr{annot_pointer}};
	my $PPI_file=$attr{annot_dir}.'BIOGRID-ALL-3.4.140.tab2.txt';
	
	#get Entrez GeneID
	my %geneid;
	foreach my $acc(keys %human){#2
		my $annot_file=$attr{annot_dir}.$attr{lib}.'_GenBank/'.$acc;
		#get annotations
		if (-f $annot_file){#4
			my $annot_pointer=sub_common::file_to_hash($annot_file, '=');
			my %annot=%$annot_pointer;
			if ($annot{'GeneID'}){
				$annot{'GeneID'}=~s/HGNC$//;
				$geneid{$annot{'GeneID'}}=$acc;
			}
		}
	}#2
	
	#get PPI 
	#my $n=0;
	open my($IN), "<", $PPI_file or die;
	#my @first_lines=split("\t", <$IN>);
	#foreach my $f(@first_lines){
	#	print "$n:$f\n";
	#	$n++;
	#}
	my %ppi;
	while (<$IN>){#2
		chomp($_);
		my @items=split("\t", $_);
		my $A_geneid=$items[1];
		my $B_geneid=$items[2];
		my $A_taxon=$items[15];
		my $B_taxon=$items[16];
		if ($A_taxon eq '9606' and $B_taxon eq '9606'){
			if (exists $geneid{$A_geneid} and exists $geneid{$B_geneid}){
				my $A=$geneid{$A_geneid};
				my $B=$geneid{$B_geneid};
				$ppi{$A}->{$B}=1;
				$ppi{$B}->{$A}=1;
				#print "$A_geneid---$B_geneid\n";
			}
			else{
				print "NO transcript id found: $A_geneid---$B_geneid\n";
			}
		}
	}#2
	close($IN);
	

	#export
	foreach my $A(keys %ppi){#3
		my $pointer=$ppi{$A};
		my %B_hash=%$pointer;
		my $B_str=join(',', sort keys %B_hash);
		#
		my $annot_file=$attr{task_dir}.$A;
		sub_basic::refresh_log($annot_file, 'PPI', $B_str);
	}
}

##################
sub human_autoantigen{
	my ($attr_pointer)=@_;
	my %attr=%$attr_pointer;
	my @acc_arr=@{$attr{pros_pointer}};
	
	#get connection between aagatlas id, ensembl
	my %disease_hash;
	open my($IN1), "<", $attr{annot_dir}."autoantigen_disease.txt" or die 'fail to open $_';
	while(<$IN1>){
		chomp($_);
		my @items=split("\t", $_);
		my $disease=($items[2] eq '\N') ? 'autoantigen' : $items[2];
		my $annot_str=$items[-1];
		if ($annot_str=~/list_uids=[0-9]* /){
			my $id=$&;
			$id=~s/list_uids=|\s//g;
			$disease_hash{$id} = ['autoantigen'] unless exists $disease_hash{$id};
			push($disease_hash{$id}, $disease) unless List::Util::first {$disease eq $_} @{$disease_hash{$id}};
			#if ($id== 4988){
			#	my $pointer=$disease_hash{$id};
			#	printf("%s:%s\t%40s\t%s\n\n", $id, $items[1], $disease, join(',',@$pointer)) ;
			#}
			#printf("%s\t%s\n", $id, $disease_hash{$id});
		}
	}
	close($IN1);

	#get Ensembl accession
	print "update annotation files\n";
	my %geneid;
	my $m=0;
	my $n=0;
	foreach my $acc(@acc_arr){#
		my $annot_file=$attr{annot_dir}.'human_GenBank/'.$acc;
		#get annotations
		my $annot_pointer=sub_common::file_to_hash($annot_file, '=');
		my %annot=%$annot_pointer;
		if (exists $annot{'GeneID'} ){
			my $annot_id=$annot{'GeneID'} ;
			$annot_id=~s/[A-Z]+$//;
			#print "$annot{'transcript_id'}:$annot_id:$disease_hash{4988}\n" if $annot_id==4988;
			#update annot_file
			if( exists $disease_hash{$annot_id}){
				my $pointer=$disease_hash{$annot_id};
				#my $disease_str='autoantigen';
				#if (@$pointer>1){
				#	shift @$pointer; #remove the first element 'autoantigen'
				#	$disease_str=join(',', @$pointer);
				#}
				my $disease_str=join(',', @$pointer);
				$annot{'autoantigen'}=$disease_str;
				sub_common::hash_to_file(\%annot, $annot_file, '=');
				printf("%s/%s\t=%s=\n", $acc, $annot{'GeneID'}, $annot{'autoantigen'});
				$m++ unless $annot{'autoantigen'} eq 'autoantigen';
				$n++;
				delete $disease_hash{$annot_id};
			}else{
				#print "No autoantigen:=$annot_id=\n";
			}
		}
	}#2
	printf("Total autoantigen=%s\tdisease autoantigen=%s\n", $n, $m);
	print "unknown:\n";
	foreach my $a(keys %disease_hash){
		print "$a\n";
	}

}
####################
sub human_combine_annot{
	my ($hash_pointer, $annot_dirs_pointer)=@_;
	my %hash=%$hash_pointer;
	my @annot_dirs=@$annot_dirs_pointer;
	#
	foreach my $annot_dir(@annot_dirs){
		printf ("Annotation directory: %s\n", $annot_dir);
		foreach my $acc(keys %hash){#3
				my $annot_file=$annot_dir.$acc;
				#get annotations
				if (-f $annot_file){#4
					my $annot_pointer=sub_common::file_to_hash($annot_file, '=');
					my %annot=%$annot_pointer;
					foreach my $key(keys %annot){#5
						if($key=~/^pro_motifs:/){#6
							my ($tag, $tag_pep_id)=split(':', $key);
							my @motifs=split(';', $annot{$key});
							my @long_motifs;
							foreach my $motif(@motifs){
								my($a, $b)=split(',', $motif);
								push(@long_motifs, $motif) if length($b)>2;
							}
							$hash{$acc}->{$tag_pep_id}->{$tag}=join(';', @long_motifs);
						}#6
						else{#6
							#enter nested hash
							my $hash2_pointer=$hash{$acc};
							my %hash2=%$hash2_pointer;
							foreach my $pep_id(keys %hash2){
								$hash{$acc}->{$pep_id}->{$key}=$annot{$key};
								#$hash{$acc}->{$pep_id}->{'pro_len'}=int($annot{$key}/3)-1 if $key eq 'rna_len';
							}
						}#6
					}#5
				}#4
				else{
					#print "$acc have not an annotation file\n";
				}
		}#3
	}#2
	return \%hash;
}


####################
sub D_combine_annot{
	my $attr_pointer=shift @_;
	my %attr=%$attr_pointer;
	my $hash_pointer=$attr{'annot_pointer'};
	my %hash=%$hash_pointer;
	
	#read annotation files
	while(@_>0){#2
		my $annot_name=shift @_;
		my $annot_dir=$attr{annot_dir}.$attr{lib}.'_'.$annot_name.'/';
		printf("%s annotation: %s\n", $attr{lib}, $annot_dir);
		#
		foreach my $acc(keys %hash){#3
			#enter nested hash
			my $hash2_pointer=$hash{$acc};
			my %hash2=%$hash2_pointer;
			#get annotations
			#print "$acc\n";
			my $annot_file=$annot_dir.$acc;
			if (-f $annot_file){#4
				#print "$annot_file\n";
				my $annot_pointer=sub_common::file_to_hash($annot_file, '=');
				my %annot=%$annot_pointer;
				#sub_data::print_hash(\%annot);
				foreach my $key(keys %annot){#5
					#inter-specie or inter-genus
					if($key=~/^InterSpecie:|^InterGenus:|^InterTaxon:/){#6
						my ($tag, $tag_pep_id)=split(':', $key);
						if ($annot{$key} eq 'none'){#7
							$hash{$acc}->{$tag_pep_id}->{$tag}=$annot{$key};
						}#7
						else{#7
							my @overlaps=split(';', $annot{$key});
							my @overlap_peps;
							foreach my $overlap(@overlaps){
								my @items=split(',', $overlap);
								my $overlap_pep_id=$items[2];
								push(@overlap_peps, $overlap_pep_id);
							}	
							$hash{$acc}->{$tag_pep_id}->{$tag}=join(',', @overlap_peps);
						}#7
						#printf("%s\t%s\t%s\n",  $acc, $tag_pep_id, $annot{$key});
					}#6
					elsif($key=~/^pro_motifs:/){#6
						my ($tag, $tag_pep_id)=split(':', $key);
						my @motifs=split(';', $annot{$key});
						my @long_motifs;
						foreach my $motif(@motifs){
							my($a, $b)=split(',', $motif);
							push(@long_motifs, $motif) if length($b)>2;
						}
						$hash{$acc}->{$tag_pep_id}->{$tag}=join(';', @long_motifs);
					}#6
					else{#6
						foreach my $pep_id(keys %hash2){
							$hash{$acc}->{$pep_id}->{$key}=$annot{$key};		
						}
					}#6
				}#5
				#print "$hash{$acc}->{'taxon_phip'}\n" if exists $hash{$acc}->{'taxon_phip'};
			}#4
			else{#4
				printf( "No annotation file in %s\n", $annot_file);
			}#4
		}#3
	}#2
	
	#export taxon
	return \%hash;
}

##############################################
#sub_phip::E_add_virus($annot_dir.'VirScanSupplement_annot.txt', \@pro_cols, $out_dir.'VirScan');
####
#hash{pro_id}->{pep_id}->{XX}=XX
sub E_add_virus{
	my($infile, $pro_cols_pointer, $file_head)=@_;
	
	#read infile
	my %h;
	open my ($IN), "<", $infile or die;
	my $first=<$IN>;
	while ( my $L=<$IN>) {#2
		chomp($L);
		my ($pep_id, $UniProt_acc, $pep_rank, $pep_start, $pep_end, $pep_dna, $pep_aa, $status, $product, $gene_symbol, $organism, 
			$pro_len, $organism_id, $taxonomic_lineage, $gene_synonym, $pro_aa, $annotation, $keywords, $features)=split("\t", $L);
		$h{$UniProt_acc}->{$pep_id}->{'UniProt_acc'}=$UniProt_acc ;
		$h{$UniProt_acc}->{$pep_id}->{'pro_id'}=$UniProt_acc ;
		$h{$UniProt_acc}->{$pep_id}->{'phip_seq'}=$pep_dna;
		$h{$UniProt_acc}->{$pep_id}->{'pep_dna'}=$pep_dna;
		$h{$UniProt_acc}->{$pep_id}->{'pep_aa'}=$pep_aa;
		$h{$UniProt_acc}->{$pep_id}->{'pep_pos'}=$pep_start.'-'.$pep_end; 
		$h{$UniProt_acc}->{$pep_id}->{'pep_start'}=$pep_start; 
		$h{$UniProt_acc}->{$pep_id}->{'pep_end'}=$pep_end; 
		$h{$UniProt_acc}->{$pep_id}->{'pep_rank'}=$pep_rank;
		$h{$UniProt_acc}->{$pep_id}->{'gene_symbol'}=$gene_symbol;
		$h{$UniProt_acc}->{$pep_id}->{'product'}=$product;
		$h{$UniProt_acc}->{$pep_id}->{'gene_synonym'}=$gene_synonym, 
		$h{$UniProt_acc}->{$pep_id}->{'pro_len'}=$pro_len;
		#printf( "##%s##\n", $L);
		#printf( "##%s##\n", $L) if length($h{$UniProt_acc}->{$pep_id}->{'pro_id'})<4;
		#last if length($h{$UniProt_acc}->{$pep_id}->{'pro_id'})<4;
	}
	close($IN);
	
	#export into text
	D_update_fa(\%h, $pro_cols_pointer, $file_head);
}

########################
sub Ztrash_SwissProt_query{#1
		my ($attr_pointer, $type)=@_;
		my %attr=%$attr_pointer;
		my %annot=%{$attr{'annot_pointer'}};
		my %pep_pro=%{$attr{'pep_pro_pointer'}};
		my $out_dir=$attr{annot_dir}.$attr{lib}.'_UniProt/';
		mkdir($out_dir, 0755) unless -d $out_dir;
		
		print "select uniprotkb_id for query\n";
		my @accs=keys %annot;
		my $accs_pointer=judge_query_acc(\@accs, $out_dir, $type);
		my @uniprot_accs=@$accs_pointer;
		
		my $p=@uniprot_accs;
		while(1){#2
			if(threads->list() < 16 and @uniprot_accs>0 ){#3
				my @sub_accs=splice @uniprot_accs, 0, 49;
				threads->create(\&C_SwissProt_query_a, $p, \@sub_accs, $out_dir);
				$p=$p-50;
			}#3
			#recover all threads
			foreach my $sub_thread( threads->list() ){#3
				$sub_thread->join() if $sub_thread->is_joinable();
			}#3
			last if threads->list()==0 and @uniprot_accs==0;
			#print @uniprot_accs;
		}#2
		
		#check the annotation files
		foreach my $pep_id(keys %pep_pro){
			my @a=split(',', $pep_pro{$pep_id});
			my $found=0;
			foreach my $acc(@a){
				$found=1 if -f $out_dir.$acc;
			}
			printf("No swissprot annotation of %s, %s\n", $pep_pro{$pep_id}, $pep_id) if $found==0;
		}
}#1
#############
sub Ztrash_SwissProt_query_a{
	my ($p, $accs_pointer, $out_dir)=@_;
	my @uniprot_accs=@$accs_pointer;
	my @taxon_arr=qw/genus family superkingdom/;
	#
	my$sp_obj = Bio::DB::SwissProt->new();
	foreach my $acc(@uniprot_accs){#2
		eval{#3
			my %pro_info;
			my $seq_obj;
			if ($seq_obj = $sp_obj->get_Seq_by_acc($acc) ){#4
				$pro_info{UniProt_acc}=$seq_obj->accession_number;
				$pro_info{display_id}=$seq_obj->display_id;
				$pro_info{display_name}=$seq_obj->display_name;
				$pro_info{pro_len}=$seq_obj->length;
				$pro_info{description}=$seq_obj->desc;
				#print $pro_info{description}, "\n";
				#specie name
				$pro_info{taxon_specie}=$seq_obj->species->species();
				$pro_info{taxon_id}=$seq_obj->species->taxon->id;
				$pro_info{taxon_scientific_name}=$seq_obj->species->taxon->scientific_name;
				#other taxonomy name
				foreach my $taxon_name(@taxon_arr){#5
					eval{
						$pro_info{'taxon_'.$taxon_name}= sub_bioperl::taxon_id($pro_info{taxon_id}, $taxon_name);
					};
					#printf("No %s of %s: %s throught the taxonomic tree\n", $taxon_name, $pro_info{taxon_id}, $pro_info{uniprot_acc}) if $@;
				}#5
				#get annotations
				my $anno_collection=$seq_obj->annotation;
				foreach my $key ( $anno_collection->get_all_annotation_keys ) { #4 features circyling
					#print "$key\n";
					my @annotations=$anno_collection->get_Annotations($key);
					#print "@annotations\n";
					for my $value(@annotations){#5
						#printf("\t%s###%s\n", $value->tagname, $value->display_text);
						if ($value->tagname eq 'dblink'){#6
							my@a=split(":", $value->display_text);
							my $f_tag=shift @a;
							my $f_value=join(':', @a);
							#printf("\t%s###%s\t%s\n", $value->tagname, $f_tag, $f_value);
							$pro_info{$f_tag}= $pro_info{$f_tag} ? $pro_info{$f_tag}.','.$f_value : $f_value;
						}#6
					}#5
				}#4
				#export
				sub_common::hash_to_file(\%pro_info, $out_dir.$acc, '=') unless keys %pro_info==0;
				printf("%d:%s,%s\n", $p, $pro_info{'taxon_superkingdom'}, $acc);
			}#4
		};#3
		printf("No annotations of %s found.\n",$acc) if $@; #return error when query
		$p--;
	}#2
	
}
###########################
#Old code and no use now
#taxon rank would be specie, genus or family
sub virus_update_taxon{
	my($virus_pointer, $taxon_file, $taxon_rank)=@_;
	my %virus=%$virus_pointer;
	
	#update
	my $new_taxon_pointer=sub_common::file_to_hash2($taxon_file, ',');
	my %new_taxon=%$new_taxon_pointer;
	#
	my $n=0;
	my %virus_taxon;
	foreach my $acc(keys %virus){#2
		my$pep_pointer=$virus{$acc};
		my %pep=%$pep_pointer;
		foreach my $pep_id(keys %pep){ #3
			my $uniprot_specie=$virus{$acc}->{$pep_id}->{'uniprot_specie'};
			#print "$uniprot_specie\n";
			my $out_taxon=$new_taxon{$uniprot_specie}->{$taxon_rank};
			if ($out_taxon eq 'NA' or $out_taxon=~/_exclude/){#4
				print "Skip specie: $uniprot_specie\n";
			}
			else{
				#printf("Accession: %s \t Taxonomy: %s, %s\n", $acc, $uniprot_specie, $out_taxon);
				#update %virus_taxon
				$virus_taxon{$acc}=$out_taxon;
				$n++;
			}#4
			last;
		}#3
	}#2
	delete $virus_taxon{'control_pep'} if $virus_taxon{'control_pep'};
	printf("\nTotal number of %d proteins by taxonomy %s\n\n ",$n, $taxon_rank);
	return(\%virus_taxon);
}
############3#
sub Ztrash_personal_taxonomy{
	my($virus_pointer, $annot_dir)=@_;
	my %hash=%$virus_pointer;
	my $uniprot_dir=$annot_dir.'virus_UniProt/';
	
	#phip_genus vs  phip_taxon
	my $taxon_file=$annot_dir.'20150414_FamilyConversion.csv';
	my $pointer=sub_common::file_to_hash($taxon_file, ',');
	my %taxon=%$pointer;
	#unirpot_genus vs phip_genus
	$taxon_file=$annot_dir.'20150323_GenusConversion.csv';
	$pointer=sub_common::file_to_hash($taxon_file, ',');
	my %genus=%$pointer;
	#uniprot_specie vs phip_specie
	$taxon_file=$annot_dir.'20150319_SpeciesConversion.csv';
	$pointer=sub_common::file_to_hash($taxon_file, ',');
	my %specie=%$pointer;
	#
	my %taxons;
	my $n=keys %hash;
	#print "$n\n";
	foreach my $acc(keys %hash){
		print" $acc\n" if $acc=~/SPIK/;
		my $hash2_pointer=sub_common::file_to_hash($uniprot_dir.$acc, '=');
		my %hash2=%$hash2_pointer;
		my $uniprot_specie=$hash2{'uniprot_specie'} ? $hash2{'uniprot_specie'} : 'NA';
		my $uniprot_genus=$hash2{'uniprot_genus'} ? $hash2{'uniprot_genus'} : 'NA';
		#print "$uniprot_specie:$uniprot_genus\n";
		$taxons{$uniprot_specie}->{'uniprot_genus'}=$uniprot_genus;
		$taxons{$uniprot_specie}->{'uniprot_subfamily'}=$hash2{'uniprot_subfamily'} ? $hash2{'uniprot_subfamily'} : 'NA';
		$taxons{$uniprot_specie}->{'uniprot_family'}=$hash2{'uniprot_family'} ? $hash2{'uniprot_family'} : 'NA';
		$taxons{$uniprot_specie}->{'phip_specie'}= $specie{ $uniprot_specie} ? $specie{ $uniprot_specie} : $uniprot_specie;
		$taxons{$uniprot_specie}->{'phip_genus'}=$genus{ $uniprot_genus} ? $genus{ $uniprot_genus} : $uniprot_genus;
		$taxons{$uniprot_specie}->{'phip_taxon'}=$taxon{$taxons{$uniprot_specie}->{'phip_genus'}} ? $taxon{$taxons{$uniprot_specie}->{'phip_genus'}}: $taxons{$uniprot_specie}->{'phip_genus'};
		#printf("%s, %s,%s, %s\n", $acc, $taxons{$uniprot_specie}->{'phip_specie'}, $taxons{$uniprot_specie}->{'phip_genus'},  $taxons{$uniprot_specie}->{'phip_taxon'});
	}
	my $taxon_outfile=$annot_dir.'virus_personal_taxonomy.csv';
	sub_common::hash2_to_file(\%taxons, $taxon_outfile, ',', 'uniprot_specie');
	#
}
############3#
sub personal_taxonomy_v2{
	my($annot_dir)=@_;
	#printf("###%s###\n", $annot_dir);
	
	#phip_genus vs  phip_taxon
	my $taxon_file=$annot_dir.'20171222_FamilyConversion.csv';
	my $pointer=sub_common::file_to_hash($taxon_file, ',');
	my %taxon=%$pointer;

	#
	foreach my $acc(keys %taxon){
		my $infile=$annot_dir.'virus_UniProt/'.$acc;
		#printf("###%s###\n", $infile);
		#update phip_taxon
		 if (-f $infile){
		 	#printf("%s:%s\n", $acc, $taxon{$acc});
			sub_basic::refresh_log($infile,'taxon_phip', $taxon{$acc}, '=');
		}else{
			printf($infile);
		}
	}
	#
}

########################
sub virusS3_interpep_query{#1
	my($attr_pointer, $taxon_rank, $out_tag, $accs_pointer)=@_;
	my %attr=%$attr_pointer;
	my $virus_pointer=$attr{annot_pointer};
	my @pro_ids=@$accs_pointer;
	
	#update annotations from swissprot
	my $virus_taxon_pointer=D_combine_annot($attr_pointer, 'UniProt');
	my %virus_taxon=%$virus_taxon_pointer;
	#update taxon rank and taxon name  #old code#
	#my $virus_taxon_pointer=virus_update_taxon($virus_pointer, $taxon_file, $taxon_rank);
	
	print "read depedent peptids\n";
	my %dependent;
	open my($IN), "<", $attr{dependent_file} or die;
	while(<$IN>){
		chomp($_);
		my ($pro_id1, $pep_id1, $pro_id2, $pep_id2, $pep) =split(',', $_);
		#print "$_\n";
		$dependent{$pep_id1}->{$pep_id2}=join(',', $pro_id2, $pep_id2, $pep);
	}
	close($IN);
	
	#multi-threading
	sub_common::multi_threading(\&virusS3_interpep_query_a, 8, \@pro_ids, $virus_pointer, \%virus_taxon, \%dependent, $taxon_rank, $out_tag, $attr{task_dir});
	#
}#1
######################
sub virusS3_interpep_query_a{
	my($in_pro_id, $virus_pointer, $virus_taxon_pointer, $dependent_pointer, $new_taxon_rank, $out_tag, $out_dir)=@_;
	my %dependent=%$dependent_pointer;
	#target virus hash
	my %virus_taxon=%$virus_taxon_pointer;
	my $in_taxon=$virus_taxon{$in_pro_id};
	#print "####$in_pro_id:$virus_taxon{$in_pro_id}\n";
	#
	my %inter_taxon;
	$inter_taxon{$new_taxon_rank}=$in_taxon;
	#the pool of pro_ids
	my @out_pro_ids;
	foreach my $out_pro_id(keys %virus_taxon){
		push(@out_pro_ids, $out_pro_id) unless $virus_taxon{$out_pro_id} eq $in_taxon;
	}
	
	#
	my %virus=%$virus_pointer;
	#get peptides belonging to $in_pro_id
	my $in_peps_pointer=$virus{$in_pro_id};
	my %in_peps=%$in_peps_pointer;
	my $pep_num=keys %in_peps;
	printf("\t%s: %s specific of %s peptides: %s", $new_taxon_rank, $in_taxon, $pep_num, join(',', keys %in_peps) );
	#inter-specie searching
	foreach my $in_pep_id(keys %in_peps){#2
		my @dependent_arr;
		#my $found=0;
		foreach my $out_pro_id(@out_pro_ids){#3
			my $out_peps_pointer=$virus{$out_pro_id};
			my %out_peps=%$out_peps_pointer;
			foreach my $out_pep_id(keys %out_peps){#4
				 if (exists $dependent{$in_pep_id}->{$out_pep_id}){#5
				 	my $overlap=join(',', $out_peps{$out_pep_id}->{taxon_specie}, $dependent{$in_pep_id}->{$out_pep_id});
				 	push(@dependent_arr, $overlap);
				 	#$found=1;
				 }#5
				#last if $found==1;
			}#4
			#last if $found==1;
		}#3
		#export
		my $key=$out_tag.':'.$in_pep_id;
		$inter_taxon{$key}= (@dependent_arr==0) ? 'none': join(';', @dependent_arr);
	}#2
	#export
	my $annot_file=$out_dir.$in_pro_id;
	sub_common::hash_to_file(\%inter_taxon, $annot_file, '=');
	printf("\tSave into %s\n", $annot_file);
	#
}

#####################################333
sub virusS4_interpep_stat{
	my($attr_pointer, $taxon_rank, $inter_name)=@_;
	my %attr=%$attr_pointer;
	my $virus_pointer=$attr{annot_pointer};
	my %virus=%$virus_pointer;
	#combine annotations
	$virus_pointer=D_combine_annot($attr_pointer, 'UniProt');

	print "read depedent peptids\n";
	my %dependent;
	open my($IN), "<", $attr{dependent_file} or die;
	while(<$IN>){
		chomp($_);
		my ($pro_id1, $pep_id1, $pro_id2, $pep_id2, $pep) =split(',', $_);
		#print "$_\n";
		$dependent{$pep_id1}->{$pep_id2}=join(',', $pro_id2, $pep_id2, $pep);
	}
	close($IN);
	
	#print keys %virus;
	#get inter counting
	my (%pep_taxon, %taxon_pep);
	foreach my $pro_id (sort (keys %virus) ){#2
		my $pointer=$virus{$pro_id};
		my %subh=%$pointer;
		foreach my $pep_id(sort (keys %subh)){ #3
			my $taxon_rank_name=$subh{$pep_id}->{$taxon_rank} ? $subh{$pep_id}->{$taxon_rank} : 'NA';
			if($taxon_rank_name eq 'NA' or $taxon_rank_name eq ''){
				printf("No taxonomy %s, %s, %s\n", $pro_id, $pep_id, $taxon_rank_name);
			}
			else{
				$pep_taxon{$pep_id}=$taxon_rank_name;
				$taxon_pep{$taxon_rank_name}= $taxon_pep{$taxon_rank_name} ? $taxon_pep{$taxon_rank_name}.','.$pep_id : $pep_id;
			}
		}#3
	}#2
	my @taxons=keys %taxon_pep;
		
	#
	my (%taxon_counts,%per_taxon_counts);
	foreach my $row(@taxons){#2
		print "$row\n";
		my @row_peps=split(',', $taxon_pep{$row});
		foreach my $col(@taxons){#3
			my @col_peps=split(',', $taxon_pep{$col});
			my $num=0;
			foreach my $a(@row_peps){#4
				foreach my $b(@col_peps){#5
					if ($dependent{$a}->{$b} or $b eq $a){
						$num++;
						last;
					}
				}#5
			}#4
			$taxon_counts{$row}->{$col}=$num;
			$per_taxon_counts{$row}->{$col}=$num/@row_peps;
		}#3
	}#2

	my $outfile=$attr{out_dir}.$inter_name.'_peptides_counting.csv';
	printf ("\nExport annotations into the file %s\n\n: ", $outfile);
	sub_common::hash2_to_file(\%taxon_counts, $outfile, ',');
	$outfile=$attr{out_dir}.$inter_name.'_peptides_counting_percentage.csv';
	printf ("\nExport annotations into the file %s\n\n: ", $outfile);
	sub_common::hash2_to_file(\%per_taxon_counts, $outfile, ',');
	
}

#########################
sub A_judge_query_acc{
	my ($attr_pointer)=@_;
	my %attr=%$attr_pointer;
	
	my (@query_accs, @skip);
	#remove unknown proteins
	my @accs=(exists $attr{pros_pointer}) ? @{$attr{pros_pointer}} : keys %{$attr{annot_pointer}};
	@accs=grep {$_ !~ /\\|\/|;|^chr|^TDP|^PTB|^k99|^gi_/} @accs;
	foreach my $acc(@accs){#3
		my $annot_file=$attr{task_dir}.$acc;
		if (-f $annot_file){
			push(@skip, $acc);
			#printf( "\tskip %s\n", $acc);
		}else{	
			push(@query_accs, $acc);	
		}
	}#3

	my $sk=@skip;
	printf( "\n\nskipped proteins: %s\n", $sk);
	my $total=@query_accs;
	printf( "Querys: %s\n\n", $total);

	return(\@query_accs);
}

#######################
sub aa_overlapped{
	my($pep_seq)=@_;
	
	my @aa_arr;
	if (length($pep_seq)>6){#4
		for (my $i=0; $i<length($pep_seq)-6; $i++){
			my $epi=substr($pep_seq, $i, 7);
			unless (List::Util::first {$epi eq $_} @aa_arr){# remove duplicated epi-aa
				push(@aa_arr, $epi) unless $epi=~/[^A-Za-z]/; #remove wrong aa like AT*TT
			}
		}
	}
	else{
		push(@aa_arr, $pep_seq );
	}
	return(\@aa_arr);
}
######################
sub aa_matching{
	my($virus_pointer, $virus_species_pointer, $pep_ids_pointer, $aa_arr_pointer)=@_;
	my %virus=%$virus_pointer;
	my %virus_species=%$virus_species_pointer;
	my @pep_ids=@$pep_ids_pointer;
	my @aa_arr=@$aa_arr_pointer;
	
	#
	my $inter_out='none';
	foreach my $aa_stretch(@aa_arr){#2
		foreach my $pro_id(@pep_ids){#3
			my $specie=$virus_species{$pro_id};
			my $pep_hash_pointer=$virus{$pro_id};
			my %pep_hash=%$pep_hash_pointer;
			foreach my $pep_id(keys %pep_hash){#4
				my $pep_aa=$pep_hash{$pep_id}->{pep_aa};
				if ($pep_aa=~/$aa_stretch/){#5
					$inter_out=$specie.','.$pro_id.','.$pep_id.','.$aa_stretch;
					#printf("%s:\t%s\n", $specie, $inter_out);
					last;
				}#5
			}#4
			last unless $inter_out eq 'none';
		}#3
		last unless $inter_out eq 'none';
	}#2
	
	return ($inter_out);
}

###############
sub A_extract_virus{
	my ($attr_pointer, $infile2)=@_;
	my %attr=%$attr_pointer;
	my(%virus, %virus_pep_pro, %virus_pro_taxon);
	
	#read the old file first
	my $n=0;
	open my ($IN), "<", $attr{'phip_fa'} or die;
	my($L1, $seq);
	while ( defined ($L1=<$IN>) && defined ($seq=<$IN>)) {#1
		chomp($L1, $seq);
		$L1=~s/^"|"$//g;
		#name
		my @items=split(/\|/, $L1);
		#print join("\t", $pep_id, $pro_id, $pep_order, $pep_type, $pro_name), "\n";
		#sequence
		$seq=~tr/a-z/A-Z/;
		my $fseq=substr($seq, 0, 51); # Note: only the first 50nt !!!
		#
		if ($L1=~/^>[1-9]/){#2
			my ($pep_id, $specie, $name, $pro_id, $pep_pos,$pep_aa)=@items;
			$pep_id=~s/^>|^">//g;
			unless($virus_pep_pro{$pep_id}){#3
				#print("$pep_pos\t");
				$virus{$pro_id}->{$pep_id}->{name}=$L1;
				$virus{$pro_id}->{$pep_id}->{phip_seq}=$fseq;
				$virus{$pro_id}->{$pep_id}->{pep_pos}=$pep_pos;
				my($start, $end)=split('_', $pep_pos);
				$virus{$pro_id}->{$pep_id}->{pep_rank}=$start;
				$virus{$pro_id}->{$pep_id}->{pep_start}=$start;
				$virus{$pro_id}->{$pep_id}->{pep_end}=$end;
				$virus{$pro_id}->{$pep_id}->{pep_aa}=$pep_aa;
				$virus{$pro_id}->{$pep_id}->{pep_dna}=$seq;
				$virus{$pro_id}->{$pep_id}->{product}=$name;
				$virus{$pro_id}->{$pep_id}->{specie}=$specie;
				#print "$virus{$pro_id}->{$pep_id}->{pep_rank}\n" if $pep_id=='9540';
				# #one pep_id vs. multiple pro_id
				$virus_pep_pro{$pep_id}=$pro_id;
				$n++;
			}#3
		}#2
		else{#2
			print "$L1\n" unless $L1=~/^>hb/;
		}#2
	}#1
	close ($IN);
	my @pros=sort(keys %virus); #protein ids 
		
	#read infile2
	my (@pros_add, %tmp);
	open my ($IN2), "<", $attr{'phip_virus1'} or die;
	my $first=<$IN2>;
	while ( my $L=<$IN2>) {#2
		chomp($L);
		my ($pep_id, $UniProt_acc, $pep_rank, $pep_start, $pep_end, $pep_dna, $pep_aa, $status, $product, $gene_symbol, $organism, 
			$pro_len, $organism_id, $taxonomic_lineage, $gene_synonym, $pro_aa, $annotation, $keywords, $features)=split("\t", $L);
		$virus{$UniProt_acc}->{$pep_id}->{'UniProt_acc'}=$UniProt_acc ;
		$virus{$UniProt_acc}->{$pep_id}->{'pro_id'}=$UniProt_acc ;
		$virus{$UniProt_acc}->{$pep_id}->{'phip_seq'}=$pep_dna;
		$virus{$UniProt_acc}->{$pep_id}->{'pep_dna'}=$pep_dna;
		$virus{$UniProt_acc}->{$pep_id}->{'pep_aa'}=$pep_aa;
		$virus{$UniProt_acc}->{$pep_id}->{'pep_pos'}=$pep_start.'-'.$pep_end; 
		$virus{$UniProt_acc}->{$pep_id}->{'pep_start'}=$pep_start; 
		$virus{$UniProt_acc}->{$pep_id}->{'pep_end'}=$pep_end; 
		$virus{$UniProt_acc}->{$pep_id}->{'pep_rank'}=$pep_rank;
		$virus{$UniProt_acc}->{$pep_id}->{'gene_symbol'}=$gene_symbol;
		$virus{$UniProt_acc}->{$pep_id}->{'product'}=$product;
		$virus{$UniProt_acc}->{$pep_id}->{'gene_synonym'}=$gene_synonym;
		$virus{$UniProt_acc}->{$pep_id}->{'pro_len'}=$pro_len;
		#aa stretches
		$virus{$UniProt_acc}->{$pep_id}->{aa_arr_pointer}=aa_overlapped($pep_aa);
		if(exists $virus_pep_pro{$pep_id}){
			printf("%s=%s\n", $UniProt_acc, $pep_id);
		}else{
			$virus_pep_pro{$pep_id}=$UniProt_acc;
		}
		if (List::Util::first {$_ eq $UniProt_acc} @pros){
			#printf("%s=%s\n", $UniProt_acc, $pep_id);
			$tmp{$UniProt_acc}='';
		}else{
			push(@pros_add, $UniProt_acc) unless List::Util::first {$_ eq $UniProt_acc} @pros_add;
		}
	}
	close($IN2);
	push(@pros, sort @pros_add); #append additional viral protein ids
	#sub_data::print_hash(\%tmp);
	sub_common::array_to_file(\@pros_add, '/home/yuan/Downloads/VirScan_add.txt' );
	
	##########export
	my $m=0;
	foreach my $pro_id(sort(keys %virus)){
		my $p=$virus{$pro_id};
		my %pp=%$p;
		my $p_n=keys %pp;
		$m +=$p_n;
		#print "=$pro_id=,";
	}
	my $pro_num=keys %virus;
	my $pep_n=keys %virus_pep_pro;
	printf("Total %d protiens and %d (%d, %d) peptides were achieved from fasta file\n\n", $pro_num, $n, $m, $pep_n);
	$attr{annot_pointer}=\%virus;
	$attr{pep_pro_pointer}=\%virus_pep_pro;
	$attr{pros_pointer}=\@pros;
	return(\%attr);
}
#######################
sub B_pep_overlapping{
	my($attr_pointer)=@_;
	my %attr=%$attr_pointer;
	my %annot=%{$attr{'annot_pointer'}};
	
	#directories
	#extract pro_ids, get 7-aa, 
	my @pro_ids;
	foreach my $pro_id(keys %annot){#2
		my %hash2=%{$annot{$pro_id}};
		push(@pro_ids, $pro_id) unless -f $attr{task_dir}.$pro_id;
		#2: aa stretches
		foreach my $pep_id(keys %hash2){#3
			my $pep_aa=$annot{$pro_id}->{$pep_id}->{pep_aa};
			$annot{$pro_id}->{$pep_id}->{aa_arr_pointer}=aa_overlapped($pep_aa);
		}#3
	}#2
	
	#parrallel processing
	#threads num: 8
	if (@pro_ids>0){
		my $num=@pro_ids;
		printf( "\nThe total numer of proteins for searching overlapped peptides is %s\n\n", $num);
		#
		sub_common::multi_threading(\&B_pep_overlapping_a, 8, \@pro_ids, \%annot, $attr{task_dir});
	}
	
	print "combine and export dependent files\n";
	my $dependent_file=$attr{out_dir}.$attr{lib}.'_dependent_peptides.csv';
	system("cat $attr{task_dir}/*  > $dependent_file");
}
#######################
sub B_pep_overlapping_a{
	my($pro_id1, $annot_pointer, $dependent_dir)=@_;
	my %annot=%$annot_pointer;
	my %hash1=%{$annot{$pro_id1}};
	
	#
	my @dependent;
	foreach my $pep_id1(keys %hash1){#2
		my @aa1=@{$hash1{$pep_id1}->{'aa_arr_pointer'}};
		#screen the whole peptides
		foreach my $pro_id2(keys %annot){#3
			my %hash2=%{$annot{$pro_id2}};
			my @pep_id2_arr= keys %hash2 ;
			#remove $pep_id1
			@pep_id2_arr= grep {$_ cmp $pep_id1} @pep_id2_arr if $pro_id1 eq $pro_id2;
			foreach my $pep_id2(@pep_id2_arr){#4
				#print "$pep_id2\n";
				my $pep_aa2=$hash2{$pep_id2}->{'pep_aa'};
				#search overlapped peptides
				my @overlap=grep {$pep_aa2=~/$_/} @aa1;
				##append overlapped peptides to array
				unless (@overlap==0){
					my $overlap_str=join(',', $pro_id1, $pep_id1, $pro_id2, $pep_id2, @overlap);
					push(@dependent, $overlap_str); 
				}
			}#4
		}#3
	}#2
	#export
	sub_common::array_to_file(\@dependent, $dependent_dir.$pro_id1) if @dependent>0;
}

##################
sub read_pep_order{
	my ($RC_file)=@_;
	
	print "read RC_file: $RC_file \n";
	my @peps;
	open my ($IN), "<", $RC_file or die;
	while(<$IN>){
		chomp($_);
		my @items=split("\t", $_);
		my @a=split(/\|/, $items[0]);
		my $pep_name=$a[0];
		#printf("%s\n", $pep_name);
		push(@peps, $pep_name);
	}
	close($IN);
	return(\@peps);
}
##################
sub order_pep_fa{
	my ($RC_peps_pointer, $fa_file)=@_;
	my @RC_peps=@$RC_peps_pointer;
	
	#
	print "Read $fa_file\n";
	my %fa_hash;
	open my ($IN), "<", $fa_file or die;
	my($L, $seq);
	while ( defined ($L=<$IN>) && defined ($seq=<$IN>)) {
		chomp($L, $seq);
		$L=~s/^>//;
		$fa_hash{$L}=$seq;
	}
	close($IN);
	
	#export 
	open my ($OUT), ">", $fa_file.'_old' or die;
	foreach my $pep (@RC_peps ) {
		if($fa_hash{$pep}){
			print $OUT ">$pep\n", "$fa_hash{$pep}\n";
		}
		else{
			print "no $pep in fasta\n";
		}
	}
	close($OUT);
}
############
##################
sub order_pep_annot{
	my ($RC_peps_pointer, $annot_file)=@_;
	my @RC_peps=@$RC_peps_pointer;
	
	#
	print "Read $annot_file\n";
	my %hash;
	open my ($IN), "<", $annot_file or die;
	my $first=<$IN>;
	while (<$IN>) {
		chomp($_);
		my @items=split("\t", $_);
		$hash{$items[0]}=$_;
	}
	close($IN);
	
	#export 
	open my ($OUT), ">", $annot_file.'_old' or die;
	print $OUT $first;
	foreach my $pep (@RC_peps ) {
		if($hash{$pep}){
			print $OUT "$hash{$pep}\n";
		}
		else{
			print "no $pep in $annot_file\n";
		}
	}
	close($OUT);
}

###########3
sub treat_intra_peptides{
	my ($in_file, $out_file)=@_;
	
	my $num=1;
	my %hash;
	open my($IN), "<", $in_file or die;
	my $first_line=<$IN>;
	chomp($first_line);
	while (<$IN>){
		chomp($_);
		$_=~s/'//g;
		my @items=split(',', $_);
		my $pep_id=$items[2];
		my $pep_seq=$items[6];
		print "$pep_id:$items[6]\n";
		#aa stretches
		my @aa_arr;
		if (length($pep_seq)>6){
			for (my $i=0; $i<length($pep_seq)-6; $i++){
				push(@aa_arr, substr($pep_seq, $i, 7) );
			}
		}
		else{		push(@aa_arr, $pep_seq );		}
		$hash{$num}->{'line_pointer'}=\@items;
		$hash{$num}->{'pep_id'}=$pep_id;
		$hash{$num}->{'pep_seq'}=$pep_seq;
		@aa_arr=List::MoreUtils::uniq @aa_arr;
		$hash{$num}->{'aa_pointer'}=\@aa_arr;
		$num++;
	}
	close($IN);
	
	open my($OUT), ">", $out_file or die;
	print $OUT join("\t", split(',', $first_line), 'Overlapping','Dependency'), "\n";
	my $end_num=keys %hash;
	$end_num --;
	foreach my $m(1..$end_num){#2
		my $A_pep_id=$hash{$m}->{'pep_id'};
		my $A_pointer=$hash{$m}->{'aa_pointer'};
		my @A_aa=@$A_pointer;
		my %A_aa=map {$_=>1} @A_aa;
		print "Rank $m: $A_pep_id\n";
		#
		my @overlapping;
		for (my $n=1; $n<$m; $n++){#3
				my $B_pep_id=$hash{$n}->{'pep_id'};
				my $B_pointer=$hash{$n}->{'aa_pointer'};
				my @B_aa=@$B_pointer;
				# the intersection of @females and @simpsons:
				my @intersects = grep( $A_aa{$_}, @B_aa );
				push(@overlapping, join(',', $B_pep_id, @intersects) )if @intersects>0;
		}#3
		#export
		my $pointer=$hash{$m}->{'line_pointer'};
		my @line=@$pointer;
		if (@overlapping==0){ 		push(@line, 'NA', 'independent');		}
		else{		push(@line, join(';', @overlapping), 'dependent' );		}
		print $OUT join("\t", @line), "\n";
	}#2
	close($OUT);
	
}
#############
sub collapse_matrix{
	my($infile, $outfile1, $outfile2)=@_;
	
	#read hash
	open my($IN), "<", $infile or die;
	my $header=<$IN>;
	chomp($header);
	my @proteins=split("\t", $header);
	shift @proteins;
	my %pro_specie;
	foreach my $pro(@proteins){
		my @items=split('_', $pro);
		my $specie=$items[0].'_'.$items[1];
		$pro_specie{$pro}=$specie;
		#print "$pro\t$specie\n";
	}
	#
	my %hash;
	open my($OUT), ">", $outfile2 or die;
	while (my$line=<$IN>){#2
		chomp($line);
		my @items=split("\t", $line);
		my $pep_id=shift @items;
		my $value_sum=List::Util::sum @items;
		#print $pep_id, "\t",  $value_sum, "\n"  if  $value_sum ==0;
		print $OUT "$pep_id\n" if  $value_sum ==0;
		#
		for(my $i=0; $i<@items; $i++){
			my $perc=$items[$i];
			my $pro=$proteins[$i];
			my $spe=$pro_specie{$pro};
			if ( exists $hash{$pep_id}->{$spe} ) {
				$hash{$pep_id}->{$spe} = $perc if $perc > $hash{$pep_id}->{$spe};
			}
			else{
				$hash{$pep_id}->{$spe}=$perc;
			}
		}
		#print "$line\n";
	}#2
	close($IN);
	close($OUT);
	printf ("export to %s", $outfile1);
	sub_common::hash2_to_file(\%hash, $outfile1, "\t", 'pep_id');
}
####################3

sub A_extract_toxome{
	my($attr_pointer)=@_;
	my %attr=%$attr_pointer;
	
	#read infile
	my (%annot, %pep_pro, %pro);
	my $n=0;
	open my ($IN1), "<", $attr{phip_toxome1} or die;
	my $first=<$IN1>;
	while ( my $L=<$IN1>) {#2
		chomp($L);
		my($order, $UniProt_acc, $pep_position, $pep_id, $pep_dna, $pep_aa,  $gene_symbol, 
			$status, $product, $gene_names, $organism, $pro_len, $organism_id, $taxonomic_lineage, 
			$gene_synonym, $pro_seq,$annotation,$keywords)=split("\t", $L);
		$annot{$UniProt_acc}->{$pep_id}->{'phip_seq'}=$pep_dna if $pep_dna;
		$annot{$UniProt_acc}->{$pep_id}->{'pep_dna'}=$pep_dna if $pep_dna;
		$annot{$UniProt_acc}->{$pep_id}->{'seq_len'}=length($pep_dna);
		$annot{$UniProt_acc}->{$pep_id}->{'pep_aa'}=$pep_aa if $pep_aa;
		$annot{$UniProt_acc}->{$pep_id}->{'pep_pos'}=$pep_position; 
		if (exists $pro{$UniProt_acc}){
			$pro{$UniProt_acc} ++;
			$annot{$UniProt_acc}->{$pep_id}->{'pep_rank'}=$pro{$UniProt_acc};
		}else{
			$annot{$UniProt_acc}->{$pep_id}->{'pep_rank'}=1; 
			$pro{$UniProt_acc}=1;
		}
		$annot{$UniProt_acc}->{$pep_id}->{'gene_symbol'}=$gene_symbol if $gene_symbol;
		$annot{$UniProt_acc}->{$pep_id}->{'product'}=$product if $product;
		$annot{$UniProt_acc}->{$pep_id}->{'gene_synonym'}=$gene_synonym if $gene_synonym;
		$annot{$UniProt_acc}->{$pep_id}->{'pro_len'}=$pro_len if $pro_len;
		#printf( "##%s:%s:%s\n", $ll, $order, $L) if length($annot{$UniProt_acc}->{$pep_id}->{pep_aa})<10;
		$pep_pro{$pep_id}=$UniProt_acc;
		#last if $n==5;
		$n++;
	}
	close($IN1);
	
	open my ($IN2), "<", $attr{phip_toxome2} or die;
	$first=<$IN2>;
	while ( my $L=<$IN2>) {#2
		chomp($L);
		my($pep_id, $UniProt_acc, $pep_rank, $pep_start, $pep_end, $pep_dna, $pep_aa, $status, 
			$product, $gene_names, $organism, $pro_len, $organism_id, $taxonomic_lineage, 
			$gene_synonym, $pro_seq,$annotation,$keywords, $features)=split("\t", $L);
		$annot{$UniProt_acc}->{$pep_id}->{'phip_seq'}=$pep_dna if $pep_dna;
		$annot{$UniProt_acc}->{$pep_id}->{'pep_dna'}=$pep_dna if $pep_dna;
		$annot{$UniProt_acc}->{$pep_id}->{'seq_len'}=length($pep_dna);
		$annot{$UniProt_acc}->{$pep_id}->{'pep_aa'}=$pep_aa if $pep_aa;
		$annot{$UniProt_acc}->{$pep_id}->{'pep_rank'}=$pep_rank;
		$annot{$UniProt_acc}->{$pep_id}->{'pep_start'}=$pep_start; 
		$annot{$UniProt_acc}->{$pep_id}->{'pep_end'}=$pep_end; 
		$annot{$UniProt_acc}->{$pep_id}->{'pep_pos'}=$pep_start.'-'.$pep_end; 
		$annot{$UniProt_acc}->{$pep_id}->{'product'}=$product if $product;
		$annot{$UniProt_acc}->{$pep_id}->{'gene_synonym'}=$gene_synonym if $gene_synonym;
		$annot{$UniProt_acc}->{$pep_id}->{'pro_len'}=$pro_len if $pro_len;
		#printf( "##%s:%s:%s\n", $ll, $order, $L) if length($annot{$UniProt_acc}->{$pep_id}->{pep_aa})<10;
		$pep_pro{$pep_id}=$UniProt_acc;
		#last if $n==5;
		$n++;
	}
	close($IN2);
	#
	$attr{'annot_pointer'}=\%annot;
	$attr{'pep_pro_pointer'}=\%pep_pro;
	return(\%attr);
}

####################3

sub A_extract_allergome{
	my($attr_pointer)=@_;
	my %attr=%$attr_pointer;
	
	#read infile
	my (%annot, %pep_pro, %pro);
	my $n=0;
	open my ($IN1), "<", $attr{phip_allergome1} or die;
	my $first=<$IN1>;
	while ( my $L=<$IN1>) {#2
		chomp($L);
		my($order, $UniProt_acc, $pep_position, $pep_id, $pep_dna, $pep_aa,  $gene_symbol, 
			$status, $product, $gene_names, $organism, $pro_len, $organism_id, $taxonomic_lineage, 
			$gene_synonym, $pro_seq,$annotation,$keywords, $features)=split("\t", $L);
		$annot{$UniProt_acc}->{$pep_id}->{'phip_seq'}=$pep_dna if $pep_dna;
		$annot{$UniProt_acc}->{$pep_id}->{'pep_dna'}=$pep_dna if $pep_dna;
		$annot{$UniProt_acc}->{$pep_id}->{'seq_len'}=length($pep_dna);
		$annot{$UniProt_acc}->{$pep_id}->{'pep_aa'}=$pep_aa if $pep_aa;
		$annot{$UniProt_acc}->{$pep_id}->{'pep_pos'}=$pep_position; 
		if (exists $pro{$UniProt_acc}){
			$pro{$UniProt_acc} ++;
			$annot{$UniProt_acc}->{$pep_id}->{'pep_rank'}=$pro{$UniProt_acc};
		}else{
			$annot{$UniProt_acc}->{$pep_id}->{'pep_rank'}=1; 
			$pro{$UniProt_acc}=1;
		}
		$annot{$UniProt_acc}->{$pep_id}->{'gene_symbol'}=$gene_symbol if $gene_symbol;
		$annot{$UniProt_acc}->{$pep_id}->{'product'}=$product if $product;
		$annot{$UniProt_acc}->{$pep_id}->{'gene_synonym'}=$gene_synonym if $gene_synonym;
		$annot{$UniProt_acc}->{$pep_id}->{'pro_len'}=$pro_len if $pro_len;
		#printf( "##%s:%s:%s\n", $ll, $order, $L) if length($annot{$UniProt_acc}->{$pep_id}->{pep_aa})<10;
		$pep_pro{$pep_id}=$UniProt_acc;
		#last if $n==5;
		$n++;
	}
	close($IN1);
	#
	open my ($IN2), "<", $attr{phip_allergome2} or die;
	$first=<$IN2>;
	while ( my $L=<$IN2>) {#2
		chomp($L);
		my($pep_id, $UniProt_acc, $pep_rank, $pep_start, $pep_end, $pep_dna, $pep_aa, $status, 
			$product, $gene_names, $organism, $pro_len, $organism_id, $taxonomic_lineage, 
			$gene_synonym, $pro_seq,$annotation,$keywords, $features)=split("\t", $L);
		$annot{$UniProt_acc}->{$pep_id}->{'phip_seq'}=$pep_dna if $pep_dna;
		$annot{$UniProt_acc}->{$pep_id}->{'pep_dna'}=$pep_dna if $pep_dna;
		$annot{$UniProt_acc}->{$pep_id}->{'seq_len'}=length($pep_dna);
		$annot{$UniProt_acc}->{$pep_id}->{'pep_aa'}=$pep_aa if $pep_aa;
		$annot{$UniProt_acc}->{$pep_id}->{'pep_rank'}=$pep_rank;
		$annot{$UniProt_acc}->{$pep_id}->{'pep_start'}=$pep_start; 
		$annot{$UniProt_acc}->{$pep_id}->{'pep_end'}=$pep_end; 
		$annot{$UniProt_acc}->{$pep_id}->{'pep_pos'}=$pep_start.'-'.$pep_end; 
		$annot{$UniProt_acc}->{$pep_id}->{'product'}=$product if $product;
		$annot{$UniProt_acc}->{$pep_id}->{'gene_synonym'}=$gene_synonym if $gene_synonym;
		$annot{$UniProt_acc}->{$pep_id}->{'pro_len'}=$pro_len if $pro_len;
		#printf( "##%s:%s:%s\n", $ll, $order, $L) if length($annot{$UniProt_acc}->{$pep_id}->{pep_aa})<10;
		$pep_pro{$pep_id}=$UniProt_acc;
		#last if $n==5;
		$n++;
	}
	close($IN2);
	#
	$attr{'annot_pointer'}=\%annot;
	$attr{'pep_pro_pointer'}=\%pep_pro;
	return(\%attr);
}


####################3
sub A_extract_provirome{
	my($attr_pointer)=@_;
	my %attr=%$attr_pointer;
	
	#read infile
	my (%annot, %pro, %pep_pro);
	my $n=0;
	open my ($IN), "<", $attr{phip_provirome} or die;
	my $first=<$IN>;
	while ( my $L=<$IN>) {#2
		chomp($L);
		my($order, $pep_id, $pep_dna, $pep_aa, $pro_id, $pep_position, $UniProt_acc, 
			$product, $gene_symbol, $organism, $pro_len) = split(",", $L);
		$pro_id=~s/^UniRef90_//;
		$annot{$pro_id}->{$pep_id}->{'pro_id'}=$pro_id;
		$annot{$pro_id}->{$pep_id}->{'phip_seq'}=$pep_dna if $pep_dna;
		$annot{$pro_id}->{$pep_id}->{'pep_dna'}=$pep_dna if $pep_dna;
		$annot{$pro_id}->{$pep_id}->{'seq_len'}=length($pep_dna);
		$annot{$pro_id}->{$pep_id}->{'pep_aa'}=$pep_aa if $pep_aa;
		$annot{$pro_id}->{$pep_id}->{'pep_pos'}=$pep_position; 
		if (exists $pro{$pro_id}){
			$pro{$pro_id} ++;
			$annot{$pro_id}->{$pep_id}->{'pep_rank'}=$pro{$pro_id};
		}else{
			$annot{$pro_id}->{$pep_id}->{'pep_rank'}=1; 
			$pro{$pro_id}=1;
		}
		$annot{$pro_id}->{$pep_id}->{'gene_symbol'}=$gene_symbol if $gene_symbol;
		$annot{$pro_id}->{$pep_id}->{'product'}=$product if $product;
		$annot{$pro_id}->{$pep_id}->{'pro_len'}=$pro_len if $pro_len;
		#printf( "##%s:%s:%s\n", $ll, $order, $L) if length($annot{$pro_id}->{$pep_id}->{pep_aa})<10;
		#printf("%s:%s\n", $pro_id, $pep_id);
		$pep_pro{$pep_id}=$pro_id;
		#last if $n==5;
		$n++;
	}
	close($IN);
	#export
	$attr{annot_pointer}=\%annot;
	$attr{pep_pro_pointer}=\%pep_pro;
	return(\%attr);
}


####################3

sub A_extract_mouse{
    my($attr_pointer)=@_;
    my %attr=%$attr_pointer;
    
    #read infile
    my (%annot, %pep_pro);
    my $n=0;
    open my ($IN1), "<", $attr{phip_mouse} or die;
    my $first=<$IN1>;
    while ( my $L=<$IN1>) {#2
        chomp($L);
        my($acc, $pep_pos, $pep_id, $pep_dna, $pep_aa, $product, $gene_symbol)=split("\t", $L);
        $annot{$acc}->{$pep_id}->{'phip_seq'}=$pep_dna if $pep_dna;
        $annot{$acc}->{$pep_id}->{'pep_dna'}=$pep_dna if $pep_dna;
        $annot{$acc}->{$pep_id}->{'seq_len'}=length($pep_dna);
        $annot{$acc}->{$pep_id}->{'pep_aa'}=$pep_aa if $pep_aa;
        if($pep_id =~ /\|([0-9]*)-([0-9]*)/){
        	$pep_pos=join('-', $1, $2);
        	$annot{$acc}->{$pep_id}->{'pep_start'}=$1;
        	$annot{$acc}->{$pep_id}->{'pep_end'}=$2;
        	#if ($acc eq 'D6RI83'){
        	#   printf("%s:%s-%s\n", $pep_id, $annot{$acc}->{$pep_id}->{'pep_start'},$annot{$acc}->{$pep_id}->{'pep_end'});
        	#}
        }elsif($pep_id =~ /CTERM|STOP/){
        	$annot{$acc}->{$pep_id}->{'pep_start'}=-1;
            $annot{$acc}->{$pep_id}->{'pep_end'}=-1;
            #printf("%s:%s\n", $pep_id, $annot{$acc}->{$pep_id}->{'pep_start'}) if $acc eq 'D6RI83';
        }else{
        	#print "$L\n";
        }
        #printf("%s:#%s#\n", $pep_id, $annot{$acc}->{$pep_id}->{'pep_start'}) if $acc eq 'D6RI83';
        $annot{$acc}->{$pep_id}->{'pep_pos'}=$pep_pos;
        $annot{$acc}->{$pep_id}->{'gene_symbol'}=$gene_symbol if $gene_symbol;
        $annot{$acc}->{$pep_id}->{'product'}=$product if $product;
        #printf( "##%s:%s\t%s\n", $acc, $pep_id, $annot{$acc}->{$pep_id}->{'pep_start'});
        $pep_pro{$pep_id}=$acc;
        #last if $n==5;
        $n++;
    }
    close($IN1);
    
    #rank peptides
    foreach my $acc(keys %annot){
    	my %peps=%{$annot{$acc}};
    	my @sort_id=sort {$peps{$a}->{'pep_start'}<=>$peps{$b}->{'pep_start'}} keys(%peps);
    	#printf("%s\n", join(',', @sort_id) );
    	my $end_id=shift @sort_id;
    	$annot{$acc}->{$end_id}->{'pep_rank'}=@sort_id+1;
    	my $n=1;
    	for my $id(@sort_id){
    		$annot{$acc}->{$id}->{'pep_rank'}=$n;
    		#printf( "##%s\t%s\t\t##%s:%s##\n", $acc, $id, $annot{$acc}->{$id}->{'pep_rank'}, $annot{$acc}->{$id}->{'pep_start'});
    		$n++;
    	}
    	#printf( "##%s\t%s\t\t##%s:%s##\n\n", $acc, $end_id, $annot{$acc}->{$end_id}->{'pep_rank'}, $annot{$acc}->{$end_id}->{'pep_start'});
    }
    
    
    #export
    $attr{'annot_pointer'}=\%annot;
    $attr{'pep_pro_pointer'}=\%pep_pro;
    return(\%attr);
}

#############################
1;  # make sure the file returns true or require will not succeed!#
