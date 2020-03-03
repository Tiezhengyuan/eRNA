#! /usr/bin/perl -w
use strict;
use warnings;
use Archive::Extract;
use File::Find;
use File::Fetch;
use List::Util;
use List::MoreUtils;
use LWP::Simple;
#
use Bio::SeqIO;
use Bio::Seq::Quality;
use Bio::Perl;
use Bio::DB::GenBank;
use Bio::DB::Query::GenBank;
use Bio::DB::SwissProt;
use Bio::DB::EntrezGene;
use Bio::DB::Taxonomy;
use Bio::SeqFeatureI;
use Bio::Tools::Run::StandAloneBlast;
####################################################
#
#the file constains all subroutines required 
package sub_bioperl;

##############################


########################
#input taxonomy id and taxonomy level, and return taxonomy
sub taxon_id{
	my($taxon_id, $taxon_rank)=@_;
	
	my $db=Bio::DB::Taxonomy->new();
	my $node=$db->get_Taxonomy_Node(-taxonid=>$taxon_id);
	my $in_name=$node->scientific_name;
	my $in_rank=$node->rank;
	#print "$in_rank\n";
	if ($in_rank eq $taxon_rank){
		my $taxon_id=$node->id;
		my $taxon_name=$node->scientific_name;
		#printf("%s:%s:%s\n", $taxon_rank, $taxon_name, $taxon_id);
		return($taxon_name);
	}
	else{
		my $p=$node->parent_id;
		#printf ("Parent of %s is %s\n", $p, $taxon_id);
		taxon_id($p, $taxon_rank)
	}

}
########################
sub taxon_tree{
	my($taxon_id)=@_;
	my %h;
	
	my $db=Bio::DB::Taxonomy->new();
	my $node=$db->get_Taxonomy_Node(-taxonid=>$taxon_id);
	my $in_name=$node->scientific_name;
	my $in_rank=$node->rank;
	my $taxon_name=$node->scientific_name;
	print "$in_rank\t$taxon_name\n";
	$h{$in_rank}=$taxon_name;
	if($node->parent_id){
		my $parent_id=$node->parent_id;
		taxon_tree($parent_id);
	}else{
		return(\%h);
	}

}
##########
#input is accession number
sub query_EntrezGene{
	my ($transcript_id, $outdir)=@_;
	my %feature;
	#
	eval{#2
		my$query_obj = Bio::DB::EntrezGene->new();
		my $seq_obj=$query_obj->get_Seq_by_acc($transcript_id);
		#general annotations
		$feature{desc}=$seq_obj->desc() if $seq_obj->desc();
		$feature{display_id}=$seq_obj->display_id() if $seq_obj->display_id();
		$feature{'EntrezGene_acc'}=$seq_obj->accession_number if $seq_obj->accession_number;
		$feature{display_name}=$seq_obj->display_name if $seq_obj->display_name;
		$feature{seq_len}=$seq_obj->length if $seq_obj->length;
		#get annotations
		my $anno_collection=$seq_obj->annotation;
		foreach my $key ( $anno_collection->get_all_annotation_keys ) { #4 features circyling
			#print "$key\n";
			my @annotations=$anno_collection->get_Annotations($key);
			#print "@annotations\n";
			for my $value(@annotations){#5
				my $tagname=$value->tagname;
				my $tagtext=$value->display_text;
				#printf("\t%s###%s\n", $tagname, $tagtext);
				if ($tagname eq 'Function'){#6
					my $GO_term='GO:'.$tagtext;
					my $GO_cat='GO_MF';
					$feature{$GO_cat}= $feature{$GO_cat} ? $feature{$GO_cat}.','.$GO_term : $GO_term;
				}#6
				elsif ($tagname eq 'Process'){#6
					my $GO_term='GO:'.$tagtext;
					my $GO_cat='GO_BP';
					$feature{$GO_cat}= $feature{$GO_cat} ? $feature{$GO_cat}.','.$GO_term : $GO_term;
				}#6
				elsif ($tagname eq 'Component'){#6
					my $GO_term='GO:'.$tagtext;
					my $GO_cat='GO_CC';
					$feature{$GO_cat}= $feature{$GO_cat} ? $feature{$GO_cat}.','.$GO_term : $GO_term;
				}#6
				elsif (List::Util::first {$tagname eq $_} ('MIM', 'HPRD', 'HGNC', 'Exon count', 'Ensembl', 'cyto', 'Official Symbol', 'UniGene', 'Exon count') ){#6
					$feature{$tagname}=$tagtext;
				}#6
			}#5
		}#4
		$feature{'GO_Entrez'} = $feature{'GO_Entrez'} ? $feature{'GO_Entrez'}.','. $feature{'GO_BP'} : $feature{'GO_BP'} if $feature{'GO_BP'};
		$feature{'GO_Entrez'} = $feature{'GO_Entrez'} ? $feature{'GO_Entrez'}.','. $feature{'GO_CC'} : $feature{'GO_CC'} if $feature{'GO_CC'};
		$feature{'GO_Entrez'} = $feature{'GO_Entrez'} ? $feature{'GO_Entrez'}.','. $feature{'GO_MF'} : $feature{'GO_MF'} if $feature{'GO_MF'};
		printf("Get annotations of %s from EntrezGene\n", $transcript_id);
		#sub_data::print_hash(\%feature);
		sub_common::hash_to_file(\%feature, $outdir.$transcript_id, '=') if $outdir;
	};#2
	printf ("No %s found.\n", $transcript_id) if $@; #return error when query;
	#
	return(\%feature);
}

#######################
####################
sub query_GenBank{#1
	my ($accession, $outdir)=@_;
	my %info;
	
	#GenBank query;
	eval{
		my $db_obj = Bio::DB::GenBank->new;
		my $seq = $db_obj->get_Seq_by_id($accession);
		$info{GenBank_acc}=$seq->accession_number if $seq->accession_number;
		$info{description}=$seq->desc;
		$info{seq_len}=$seq->length;
		$info{seq}=$seq->seq;
		$info{taxon_specie}=$seq->species->taxon->scientific_name;
		$info{taxonid_specie}=$seq->species->taxon->id;
		eval{	$info{taxon_genus}= taxon_id($info{taxonid_specie}, 'genus');	};
		eval{	$info{taxon_subfamily}= taxon_id($info{taxonid_specie}, 'subfamily');	};
		eval{	$info{taxon_family}= taxon_id($info{taxonid_specie}, 'family');	};
		foreach my $feat ( $seq->get_SeqFeatures() ) { #3 features circyling
			my %feature;
			my $tag_name=$feat->primary_tag();
			$feature{name}=$tag_name;
			#get feature's info
			foreach my $tag ( $feat->get_all_tags() ) {  #4
				my $tag_value=join("", $feat->get_tag_values($tag));
				if($tag=~/db_xref/){
						my($a, $b)=split(":", $tag_value);
						$feature{$a}=$b;
				}
				else{	$feature{$tag}=$tag_value;	}
			}#4
			#print $tag_name, "\n";
			#print sub_data::print_hash(\%feature);
			#extract source_features
			if ( $tag_name=~/source/ ) { #4
				$info{chromosome}=$feature{chromosome} if exists $feature{chromosome};
				$info{map}=$feature{map} if exists $feature{map};
				$info{organism}=$feature{organism} if exists $feature{organism};
				#$info{taxon_id}=$feature{taxon} if exists $feature{taxon};
				$info{mol_type}=$feature{mol_type} if exists $feature{mol_type};
			}#4
			elsif ( $tag_name=~/gene/ ) { #4
				$info{gene_symbol}= $feature{gene} if $feature{gene};
				if ($feature{GeneID}){
					$info{GeneID}= $feature{GeneID} ;
					$info{GeneID}=~s/HGNC$//;
				}
				$info{gene_synonym}= $feature{gene_synonym} if $feature{gene_synonym};
				$info{note}= $feature{note} if $feature{note};
			}#4
			elsif($tag_name=~/Protein/){
				$info{pro_wt}= $feature{calculated_mol_wt} if $feature{calculated_mol_wt};
				$info{product}= $feature{product} if $feature{product};
			}
			elsif ( $tag_name=~/CDS/ ) { #4
				$info{GeneID}= $feature{GeneID} if $feature{GeneID};
				$info{GI}= $feature{GI} if $feature{GI};
				$info{product}= $feature{product} if $feature{product};
				$info{protein_id}= $feature{protein_id} if $feature{protein_id};
				$info{EC_number}= $feature{EC_number} if $feature{EC_number};
			}
		}#3
		if ($outdir){
			my $outfile=$outdir.$accession;
			sub_common::hash_to_file(\%info, $outfile, '=');
		}
	};
	printf ("No %s found.\n", $accession) if $@; #return error when query;
	#
	#sub_data::print_hash(\%info);

	return(\%info);
}

###################
#############
#$order could be omitted
sub query_SwissProt{
	my ($acc, $outdir)=@_;
	my @taxon_arr=qw/genus family superkingdom/;
	
	#
	my %pro_info;
	#eval{#3
		my$sp_obj = Bio::DB::SwissProt->new();
		my $seq_obj = $sp_obj->get_Seq_by_acc($acc);
		$pro_info{UniProt_acc}=$seq_obj->accession_number if $seq_obj->accession_number;
		$pro_info{display_id}=$seq_obj->display_id;
		$pro_info{display_name}=$seq_obj->display_name;
		$pro_info{pro_len}=$seq_obj->length;
		$pro_info{description}=$seq_obj->desc;
		#printf("%s\n",  $seq_obj->product);
	eval{#3
		#taxonomy name
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
	};#3
	
	#export
	if ($@){
		printf("\tNo %s in SwissProt.\n", $acc) ; #return error when query
		#sub_data::print_hash(\%pro_info);
	}else{
		printf("%s,%s\n", $pro_info{'taxon_superkingdom'}, $acc);
		#export
		#sub_data::print_hash(\%pro_info); #print 
		sub_common::hash_to_file(\%pro_info, $outdir.$acc, '=') if $outdir;
	}
	
	return(\%pro_info);
}

##############################################################
#extract information from rna.gbk
sub rna_gbk{#1
	my ($gbk_file)=@_;
  
	my %rna_info;
	my $in = Bio::SeqIO->new(-file => $gbk_file, -format => 'genbank');
	while ( my $seq = $in->next_seq() ) { #2 sequence circyling
		my $acc=$seq->accession();
		#$rna_info{$acc}->{display_id}=$seq->display_id();
		$rna_info{$acc}->{GI}=$seq->primary_id();
		$rna_info{$acc}->{description}=$seq->desc();
		#$rna_info{$acc}->{rna_seq}=$seq->seq();  
		$rna_info{$acc}->{rna_len}=$seq->length();
		#print "$acc\n";
		my $found=0;
		foreach my $feat ( $seq->get_SeqFeatures() ) { #3 features circyling
			my %feature;
			my $tag_name=$feat->primary_tag();
			$feature{name}=$tag_name;
			#get feature's info
			foreach my $tag ( $feat->get_all_tags() ) {  #4
				my $tag_value=join("", $feat->get_tag_values($tag));
				if($tag=~/db_xref/){
					my($a, $b)=split(":", $tag_value);
					$feature{$a}=$b;
				}
				else{
					$feature{$tag}=$tag_value;
				}
			}#4
         
			#extract source_features
			if ( $tag_name=~/source/ ) { #4
				$rna_info{$acc}->{chromosome}=exists $feature{chromosome} ? $feature{chromosome}: 'NA';
				$rna_info{$acc}->{map}=exists $feature{map} ? $feature{map} : 'NA';
				$rna_info{$acc}->{mol_type}=$feature{mol_type} if exists $feature{mol_type};
				$found++;
			}#4
			if ( $tag_name=~/gene/ ) { #4
				$rna_info{$acc}->{gene}= exists $feature{gene} ? $feature{gene} : 'NA';
				$rna_info{$acc}->{GeneID}= exists $feature{GeneID} ? $feature{GeneID} : 'NA';
				$rna_info{$acc}->{gene_synonym}=exists $feature{gene_synonym} ? $feature{gene_synonym}  : 'NA';
				$rna_info{$acc}->{note}=exists $feature{note} ? $feature{note}  : 'NA';
				$rna_info{$acc}->{HGNC}=exists $feature{HGNC} ? $feature{HGNC} : 'NA';
				$rna_info{$acc}->{HPRD}= exists $feature{HPRD} ? $feature{HPRD}  : 'NA';
				$rna_info{$acc}->{MIM}= exists $feature{MIM} ? $feature{MIM}  : 'NA';
				$found++;
			}#4
			last if $found==2;
		}#3 features circyling
		unless ($rna_info{$acc}->{mol_type} eq 'mRNA'){
				delete $rna_info{$acc};
		}
	} #2 sequence circyling
	return(\%rna_info);
}#1

##############################################################
#extract information from rna.gbk
#The table include 15columns
sub gbk_to_fa{#1
	my ($gbk_file, $result_dir, $specie)=@_;
	$result_dir=sub_common::format_directory($result_dir);
	my @annot_names=qw/display_id accession chromosome rna_map mol_type  rna_len description rna_seq/;
				  
	my %rna_feature;
	open my ($LOG), ">", $result_dir.$specie.'_RNA.log' or die;
	print $LOG "\nExtract RNA information from $gbk_file\n\n";
	open my ($RNA_txt), ">", $result_dir.$specie.'_RNA_info.txt' or die;
	print $RNA_txt join("\t", @annot_names), "\n";

	my $in = Bio::SeqIO->new(-file => $gbk_file, -format => 'genbank');
	while ( my $seq = $in->next_seq() ) { #2 sequence circyling
		my %rna_info=map {$_=>'NA'} @annot_names;
		$rna_info{display_id}=$seq->display_id();
		$rna_info{GI}=$seq->primary_id();
		$rna_info{accession}=$seq->accession();
		$rna_info{description}=$seq->desc();
		$rna_info{rna_seq}=$seq->seq();  
		$rna_info{rna_len}=$seq->length();
		#print "$rna_info{accession}\n";

		foreach my $feat ( $seq->get_SeqFeatures() ) { #3 features circyling
			my %feature;
			my $tag_name=$feat->primary_tag();
			$feature{name}=$tag_name;
			$feature{start}=$feat->start;
			$feature{end}=$feat->end;
			$feature{strand}=$feat->strand;
			$feature{seq}=substr($rna_info{rna_seq}, $feature{start}-1, $feature{end}-$feature{start}+1);
			$feature{gene}='NA';
			#get feature's info
			foreach my $tag ( $feat->get_all_tags() ) {  #4
				my $tag_value=join("", $feat->get_tag_values($tag));
				if($tag=~/db_xref/){
					my($a, $b)=split(":", $tag_value);
					$feature{$a}=$b;
				}else{
					$feature{$tag}=$tag_value;
				}
			}#4
         
			#get pseudo-genes
			if ($tag_name=~/misc_RNA/ and exists $feature{pseudo}){
				$feature{name}='pseudo';
			}
			#extract source_features
			if ( $tag_name=~/source/ ) { #4
				$rna_info{chromosome}=$feature{chromosome} if exists $feature{chromosome};
				$rna_info{rna_map}=$feature{map} if exists $feature{map};
				$rna_info{mol_type}=$feature{mol_type} if exists $feature{mol_type};
				#rna sequences information and source features 
				print $RNA_txt join("\t", map {$rna_info{$_}} @annot_names), "\n";
				#change feature name
				$feature{mol_type}=~s/\s/_/g;
				$feature{name}=$feature{mol_type};
			}elsif ($tag_name=~/ncRNA/) {#4
				$feature{name}=$feature{ncRNA_class};
			}#4
			my $feature_index=$rna_info{accession}.':'.$rna_info{rna_len}.':'.$feature{gene}.':'.$feature{start}.'_'.$feature{end};
			$rna_feature{$feature{name}}->{$feature_index}=$feature{seq};
        
			if ($tag_name=~/CDS/){#4
				#5'UTR
				if($feature{start}>10){
					my $UTR5_start=1; 
					my $UTR5_end=$feature{start}-1;
					my $UTR5_seq=substr($rna_info{rna_seq}, $UTR5_start-1, $UTR5_end-$UTR5_start+1);
					my $feature_index=$rna_info{accession}.':'.$rna_info{rna_len}.':'.$feature{gene}.':'.$UTR5_start.'_'.$UTR5_end; 
					$rna_feature{'5UTR'}->{$feature_index}=$UTR5_seq;
				}
				#3'UTR
				if ($rna_info{rna_len}-$feature{end}>10){
					my $UTR3_start=$feature{end}+1;
					my $UTR3_end=$rna_info{rna_len};
					my $UTR3_seq=substr($rna_info{rna_seq}, $UTR3_start-1, $UTR3_end-$UTR3_start+1);
					my $feature_index=$rna_info{accession}.':'.$rna_info{rna_len}.':'.$feature{gene}.':'.$UTR3_start.'_'.$UTR3_end; 
					$rna_feature{'3UTR'}->{$feature_index}=$UTR3_seq;
				}
			}#4
		}#3 features circyling
	} #2 sequence circyling
  
	print "export sequences of features\n";
	foreach my $feature_name(keys %rna_feature){
		print "$feature_name\n";
		my $pointer=$rna_feature{$feature_name};
		my %hash=%$pointer;
		open my $OUT, ">", $result_dir.$specie.'_RNA_'.$feature_name.'.fa' or die; 
		foreach my $feature_index(keys %hash){
			print $OUT ">$feature_name:$feature_index\n", "$hash{$feature_index}\n";
		}
		close($OUT);
		my $num=keys %hash;
		print $LOG "$feature_name=$num\n";
	}
	close($LOG);
}#1


############################################################
#extract information from protein.gbk
#save them into a table including 15columns
sub protein_gbk{#1
	my ($pro_gbk_file, $out_file)=@_;
	my @annot_names=qw/display_id GI accession secondary_acc description 
			pro_seq pro_len chromosome pro_map pro_start pro_end pro_note pro_product EC_number pro_weight 
			gene  gene_synonym   coded_by   geneid CDS_note      CDS_start      CDS_end/;
	
	my $num=0;
	my $CDS_num=0;
	my $protein_num=0;
	print "\nExtract protein information from $pro_gbk_file.\n\n";
	open my($OUT), ">", $out_file or die;
	print $OUT join("\t", @annot_names), "\n";
	my $in = Bio::SeqIO->new(-file => $pro_gbk_file, -format => 'genbank') or die;
	while ( my $seq = $in->next_seq() ) { #2
		my %pro_info=map {$_=>'NA'} @annot_names; #initiate

		$num++;    
		$pro_info{display_id}=$seq->display_id();
		$pro_info{GI}=$seq->primary_id();
		$pro_info{accession}=$seq->accession();
		$pro_info{description}=$seq->desc();
		$pro_info{pro_seq}=$seq->seq();  #protein sequence
		$pro_info{pro_len}=$seq->length();
		my @secondary_acc=$seq->get_secondary_accessions();   
		$pro_info{secondary_acc}=join(",", @secondary_acc) if @secondary_acc>0;
    
		my (@chromosome, @map, @product, @pro_note, @CDS_note, @MW, @gene, @gene_synonym, @coded_by, @geneid, @EC_number);
		foreach my $feat ( $seq->get_SeqFeatures() ) { #3
			my $tag_name=$feat->primary_tag();
			#extract source_features
			if ( $tag_name=~/Source/ ) { #4
				foreach my $tag ( $feat->get_all_tags() ) {  #5
					@chromosome = $feat->get_tag_values($tag) if $tag=~/chromosome/;
					@map = $feat->get_tag_values($tag) if $tag=~/map/;
				}#5
				$pro_info{chromosome}=join("", @chromosome) if @chromosome>0;
				$pro_info{pro_map}=join("", @map) if @map>0;
			}#4
           
			#extract protein_features
			if ( $tag_name=~/Protein/ ) { #4
				foreach my $tag ( $feat->get_all_tags() ) {  #5
				@product = $feat->get_tag_values($tag) if $tag=~/product/;
				@pro_note = $feat->get_tag_values($tag) if $tag=~/note/;
				@MW = $feat->get_tag_values($tag) if $tag=~/calculated_mol_wt/;
				@EC_number = $feat->get_tag_values($tag) if $tag=~/EC_number/;
				$pro_info{pro_start}=$feat->start;
				$pro_info{pro_end}=$feat->end;
				}#5
				$protein_num++;
				$pro_info{pro_product}=join("", @product) if @product>0; 
				$pro_info{pro_note}=join("", @pro_note) if @pro_note>0;
				$pro_info{pro_weight}=$MW[0] if @MW>0;
				$pro_info{EC_number}=join("", @EC_number) if @EC_number>0;      
			}#4
 
			#extract CDS_features
			if ( $tag_name=~/CDS/ ) { #4
				foreach my $tag ( $feat->get_all_tags() ) {  #5
					@gene = $feat->get_tag_values($tag) if $tag=~/gene/;
					@gene_synonym = $feat->get_tag_values($tag) if $tag=~/synonym/;
					@coded_by = $feat->get_tag_values($tag) if $tag=~/coded_by/;
					@geneid = $feat->get_tag_values($tag) if $tag=~/GeneID/;
					@CDS_note = $feat->get_tag_values($tag) if $tag=~/note/;
					$pro_info{CDS_start}=$feat->start;
					$pro_info{CDS_end}=$feat->end;
				}#5
				$CDS_num++; 
				$pro_info{gene}=join("", @gene) if @gene>0;
				$pro_info{gene_synonym}=join("", @gene_synonym) if @gene_synonym>0;  
				$pro_info{gene_synonym}=~s/;\s/,/g;  
				$pro_info{CDS_note}=join("", @CDS_note) if @CDS_note>0;
				$pro_info{coded_by}=join("", @coded_by) if @coded_by>0;
				$pro_info{geneid}=join("", @geneid) if @geneid>0;      
			}#4        
		}#3
		print $OUT join("\t", map {$pro_info{$_}} @annot_names), "\n";
		print "$num:\t", "protein=$protein_num\t", "CDS=$CDS_num\t", "$pro_info{display_id}\n";
		#summary features
		#print "$pro_info{accession}\t", "$pro_info{description}\t", "$pro_info{display_id}\t", "$pro_info{GI}\t", "$pro_info{secondary_acc}\t", "$pro_info{pro_len}\t", "$pro_info{pro_seq}\t"; 
		#source features 
		#print "$pro_info{chromosome}\t", "$pro_info{pro_map}\t"; 
		#protein features
		#print "$pro_info{EC_number}\t", "$pro_info{pro_start}\t", "$pro_info{pro_end}\t", "$pro_info{pro_note}\t", "$pro_info{pro_product}\t", "$pro_info{pro_weight}\t";
		#CDS features
		#print  "$pro_info{CDS_note}\t", "$pro_info{coded_by}\t", "$pro_info{gene}\t", "$pro_info{gene_synonym}\t", "$pro_info{geneid}\t";
	} #2 
	close($OUT);
	print "The items number: $num\n";   
	print "The number of proteins: $protein_num\n";
	print "The number of CDS regions: $CDS_num\n";
}#1


################
sub local_blast{
	my ($par_pointer, $outfile)=@_;
	my %par=%$par_pointer;
	
	#Alignment using Blast
	my $blast_obj=Bio::Tools::Run::StandAloneBlast->new(-program=>$par{program}, -F=>'F',
			-database=>$par{database},  -expect => $par{expect}, -outfile=>$par{outfile});
	#extract results
	my (%align, %align_percent);
	#input is fasta file
	my $n=1;
	my $seqio=Bio::SeqIO->new(-file=>$par{infile}, -format=>"Fasta");
	while(my $seq_obj=$seqio->next_seq()){#2
		my $query_name=$seq_obj->display_id;
		$align_percent{$query_name}={};
		print "$n: $query_name\n";
		my $report_obj=$blast_obj->blastall($seq_obj);
		while( my $result_obj = $report_obj->next_result ) {#3
			$query_name=$result_obj->query_name;
			my $query_acc=$result_obj->query_accession;
			my $query_desc=$result_obj->query_description;
			while( my $hit= $result_obj->next_hit ) {#4
				my $hit_name=$hit->name;
				my $hits_num=$hit->num_hsps;
				printf("Query name: %s\t Query accession:%s\t Description: %s\t",  $query_name, $query_acc, $query_desc);
				while( my $hsp = $hit->next_hsp ) {#5
					my $identity_percent= $hsp->percent_identity;
					my $hit_len=$hsp->length('total');
					my $hit_score=$hsp->score;
					if ( $identity_percent > $par{identity_percent}) {             
						printf("Hit name:%s\t Length:%s\t Identity_percent:%s\t Score:%s\n", $hit_name, $hit_len,  $identity_percent, $hit_score);
						$align{$query_name}->{$hit_name}->{query_acc}=$query_acc;
						$align{$query_name}->{$hit_name}->{query_desc}=$query_desc;
						$align{$query_name}->{$hit_name}->{identity_percent}=int($identity_percent+0.5);
						$align{$query_name}->{$hit_name}->{hit_len}=$hit_len;
						$align{$query_name}->{$hit_name}->{hit_score}=$hit_score;
						$align_percent{$query_name}->{$hit_name}=int($identity_percent*10+0.5)/10;
					}       
				}#5
			}#4
		}#3
		print "\n";
		#print $align_percent{$query_name}, "\n";
		$n++;
		#last if $n==100;
	}#2
	#
	sub_common::hash2_to_file(\%align_percent, $par{'outfile'}, "\t", 'query_names');
	return(\%align);
}

############################
#$out_dir default is the current directory
#$genome_type: DNA, RNA, protein in fasta or gbk format
#'DNA_fa', 'RNA_fa', 'protein_fa', 'RNA_gbk', 'protein_gbk', 'annot_gff'
#$specie: human
sub download_NCBI_genome{
	my($genome_type, $specie, $out_dir)=@_;
	my %species=('human'=>'H_sapiens', 'mouse'=>'Mus_musculus', 'rat'=>'Rattus_norvegicus', 
					'cattle'=>'Bos_taurus', 'dog'=>'Canis_familiaris', 'chicken'=>'Gallus_gallus', 'zebrafish'=>'Danio_rerio',
					'soyben'=>'Glycine_max', 'maize'=>'Zea_mays', 'rice'=>'Oryza_sativa_Japonica_Group');
	my $url_specie=($species{$specie}) ? $species{$specie} : $species{'human'};
	my $db_url='ftp://ftp.ncbi.nlm.nih.gov/genomes/'.$url_specie.'/';
	$out_dir=Cwd::getcwd() unless $out_dir;
	$out_dir=sub_common::format_directory($out_dir.$specie);
	
	#
	my @local_files;
	if($genome_type eq 'DNA_fa'){
		my $url = $db_url.'Assembled_chromosomes/seq/';
		my $web_file_names_pointer=sub_common::web_list_files($url, "file_names", "_ref_");
		my @web_file_names=@$web_file_names_pointer;
		@web_file_names=grep(/\.fa\./,  @web_file_names);
		@web_file_names=grep(/chr/,  @web_file_names);
		#print join("\n", @web_file_names);
		
		#download files
		foreach my $web_file_name(@web_file_names){
			my $web_file=$url.$web_file_name;
			my $saved_file=sub_common::web_download_file($web_file, $out_dir);
			push(@local_files, $saved_file);
		}
	}
	elsif($genome_type eq 'RNA_fa'){
		my $web_file = $db_url.'RNA/rna.fa.gz';
		my $saved_file=sub_common::web_download_file($web_file, $out_dir);
		push(@local_files, $saved_file);
	}
	elsif($genome_type eq 'RNA_gbk'){
		my $web_file = $db_url.'RNA/rna.gbk.gz';
		my $saved_file=sub_common::web_download_file($web_file, $out_dir);
		push(@local_files, $saved_file);
	}
	elsif($genome_type eq 'protein_fa'){
		my $web_file = $db_url.'protein/protein.fa.gz';
		my $saved_file=sub_common::web_download_file($web_file, $out_dir);
		push(@local_files, $saved_file);
	}
	elsif($genome_type eq 'protein_gbk'){
		my $web_file = $db_url.'protein/protein.gbk.gz';
		my $saved_file=sub_common::web_download_file($web_file, $out_dir);
		push(@local_files, $saved_file);
	}
	elsif($genome_type eq 'annot_gff'){
		#print $db_url.'GFF/', "\n";
		my $web_file=sub_common::web_list_files($db_url.'GFF/', 'file', 'top_level');
		print "$web_file\n";
		my $saved_file=sub_common::web_download_file($web_file, $out_dir);
		push(@local_files, $saved_file);
	}
	else{
		print "Error: No files downloaded\n";
	}
	#print join("\n", @local_files, "\n"); 
	return(\@local_files);
}

############################
#$out_dir default is the current directory
#$genome_type: DNA, RNA, protein in fasta or gbk format
#'DNA_fa', 'RNA_fa', 'protein_fa', 'RNA_gbk', 'protein_gbk', 'annot_gff'
sub download_Ensemble_genome{
	my($genome_type, $specie, $out_dir)=@_;
	my %species=('human'=>'homo_sapiens', 'mouse'=>'mus_musculus', 'rat'=>'rattus_norvegicus', 
					'cattle'=>'bos_taurus', 'dog'=>'canis_familiaris', 'chicken'=>'gallus_gallus');
	my $url_specie=($species{$specie}) ? $species{$specie} : $species{'human'};
	my $db_url='ftp://ftp.ensembl.org/pub/';
	$out_dir=Cwd::getcwd() unless $out_dir;
	$out_dir=sub_common::format_directory($out_dir.$specie);
	
	#
	my @local_files;
	if($genome_type eq 'DNA_fa'){
		my $url = $db_url.'current_fasta/'.$url_specie.'/dna/';
		my $web_files_pointer=sub_common::web_list_files($url, 'files', "dna.chromosome.");
		my @web_files=@$web_files_pointer;
		#print join("\n", @web_files);
		
		#download files
		foreach my $web_file(@web_files){
			my $saved_file=sub_common::web_download_file($web_file, $out_dir);
			push(@local_files, $saved_file);
		}
	}
	elsif($genome_type eq 'RNA_fa'){
		my $url=$db_url.'current_fasta/'.$url_specie.'/cdna/';
		#print "$url\n";
		my $web_file=sub_common::web_list_files($url, 'file', '.all.fa.');
		my $saved_file=sub_common::web_download_file($web_file, $out_dir);
		push(@local_files, $saved_file);
	}
	elsif($genome_type eq 'protein_fa'){
		my $url=$db_url.'current_fasta/'.$url_specie.'/pep/';
		#print "$url\n";
		my $web_file=sub_common::web_list_files($url, 'file', '.all.fa.');
		my $saved_file=sub_common::web_download_file($web_file, $out_dir);
		push(@local_files, $saved_file);
	}
	elsif($genome_type eq 'annot_gff'){
		my $url=$db_url.'current_gff3/'.$url_specie;
		#print "$url\n";
		my $web_files_pointer=sub_common::web_list_files($url, 'files');
		my @web_files=@$web_files_pointer;
		my @sub_list=grep {$_!~/\.abinitio|\.chr/} @web_files;
		@sub_list=grep {$_=~/\.gz$/} @sub_list;
		#print "###@sub_list\n";
		my $web_file=shift @sub_list;
		my $saved_file=sub_common::web_download_file($web_file, $out_dir);
		push(@local_files, $saved_file);
	}
	else{
		print "Error: No files downloaded\n";
	}
	#print join("\n", @local_files, "\n"); 
	return(\@local_files);
}

#########################
#download from PROSITE on protein motifs
################
sub PROSITE_protein_motifs{
	my( $motif_file)=@_;
	
	#download prosite.dat from PROSITE
	my $prosite_url='ftp://ftp.expasy.org/databases/prosite/prosite.dat';
	my $local_file=sub_common::web_download_file($prosite_url, '/home/yuan/Downloads/', 0);
	
	#open prosite.dat
	my $n=0;
	my @lines;
	my %prosite;
	open my($IN), "<", $local_file or die;
	while(<$IN>){#2
		chomp($_);
		if ($_=~/^\/\//){#3
			my ($accession,$description,$regular_pattern);
			my $motif="";
			foreach my $line(@lines){#4
				$line=~s/;$//;
				my ($tag, $value)=split(/\s+/, $line, 2);
				$accession=$value if $line=~/^AC/;
				$description=$value if $line=~/^DE/;
				if ($line=~/^PA/){
					$motif=$value;
					$regular_pattern="";
					#convert regular pattern
					my @a=split('-', $motif);
					foreach my $aa(@a){
						$aa=~s/\{(\D+)\}/\[\^$1\]/g; #exclude amino acid
						$aa=~s/\[(\D+)\>\]/$1\?\$/g; #exclude amino acid at C-end
						$aa=~s/x/\./g; # x amino acid
						$aa=~s/\((\d),(\d)\)/\{$1,$2\}/g; #repeat amino acid
						$aa=~s/\((\d)\)/\{$1\}/g;#repeat amino acid
						$aa=~s/>$/\$/; # C-end of peptide
						$aa=~s/^</\^/; # N-end of peptide
						$regular_pattern .=$aa;
					}
				}
			}#4
			if (length($motif) >2){
				$prosite{$accession}->{description}=$description;
				$prosite{$accession}->{motif}=$motif;
				$prosite{$accession}->{regular_pattern}=$regular_pattern;
				print join("\t\t", $accession, $motif, $regular_pattern), "\n";
				$n++;
			}
			undef @lines;
		}#3
		else{#
			push(@lines, $_);
		}
		#last if $n==20;
	}#2
	close($IN);
	
	#export
	printf ("Total number of motifs: %d", $n);
	sub_common::hash2_to_file(\%prosite, $motif_file, "\t");
	
}


############################################################
sub pre_references_info{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my %ref_info;
	
	#read old log file
	if( -f $variables{file_references_txt}){ 
		open my($IN), "<", $variables{file_references_txt}or die;
		my $firstline=<$IN>;
		while(<$IN>){
			chomp($_);
			my($index_name, $seq_num, $total_base, $ave_len)=split("\t", $_);
			$ref_info{$index_name}->{seq_num}=$seq_num;
			$ref_info{$index_name}->{total_base}=$total_base;
			$ref_info{$index_name}->{ave_len}=$ave_len;
		}
		close($IN);
	}
	
	#generate the new one
	my $files_pointer=sub_common::files_list($variables{dir_bowtie}, 'file');
	my @fasta_files=grep(/\.fasta$|\.fa$/, @$files_pointer);

	open my($OUT), ">", $variables{file_references_txt} or die;
	print $OUT join("\t", 'index_name', 'seq_num', 'total_base', 'ave_len'), "\n";
	foreach my $fasta_file(@fasta_files){#2
		print "\t$fasta_file\n";
		my $index_name=sub_common::file_operation($fasta_file, 'name_head');
		unless (exists $ref_info{$index_name}){#3
			my $seq_num=0;
			my $total_base=0;
			my $in_obj = Bio::SeqIO->new(-file => $fasta_file, -format => 'fasta');
			while (my $seq_obj = $in_obj->next_seq() ) {#3
				my $displayid=$seq_obj->display_id();
				my $seq_len=$seq_obj->length();
				$seq_num++;
				$total_base += $seq_len;
			}#3
			$ref_info{$index_name}->{seq_num}=$seq_num;
			$ref_info{$index_name}->{total_base}=$total_base;
			$ref_info{$index_name}->{ave_len}=int($total_base/$seq_num+0.5);
		}#3
		print $OUT join("\t", $index_name, $ref_info{$index_name}->{seq_num}, 
			$ref_info{$index_name}->{total_base}, $ref_info{$index_name}->{ave_len}), "\n";
	}#2
	close($OUT);

	return(\%ref_info);
}

#################
#end
1;