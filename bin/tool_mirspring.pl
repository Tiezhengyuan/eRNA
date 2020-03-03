#!/usr/bin/perl -w
use warnings;
use strict;
use Cwd;
use List::Util;
use Gtk2::Pango;
use threads 'exit' => 'threads_only';
use threads::shared;
use Glib qw/TRUE FALSE/;
use Gtk2 qw/-init -threads-init/;
use Bio::SeqIO;
use Bio::Perl;


############################################################
#get the directory of perl scripts involved in Pscore
our $perl_dir=Cwd::getcwd();

#get subroutines 
require $perl_dir."/eRNA_subroutines.pm";

#initiate variables
my $var_file=$ARGV[0];
#my $var_file='/home/yuan/mysql_pre/eRNA/result/variables.txt';
our $variables_pointer=E_RNA::process_info($var_file);
our %variables=%$variables_pointer;
$variables{mirspring_dir}=$perl_dir.'/mirspring';
$variables{mirspring_intermediate_script}=$perl_dir.'/mirspring/BAM_to_Intermediate.pl';
$variables{mirspring_script}=$perl_dir.'/mirspring/Intermediate_to_miRspring.pl';
#initiate envirinmental variables
$ENV{'PATH'} .= ":$variables{alignment_dir}:$variables{ref_seq_dir}:$variables{mirspring_dir}";
print $ENV{'PATH'}, "\n\n\n";

#setup shared hash to control thread
my %shash;
share(%shash); #will work for first level keys
$shash{'go'} = 0;
$shash{'data'} = '';
$shash{'work'} = '';
$shash{'die'} = 0;

our (@fasta_names, @selected_names, $selected_A, $selected_B);
#get index information
my $file_name_pointer=E_RNA::files_list($variables{alignment_dir}, 'file_names');
my @file_name=@$file_name_pointer;
my @index_names=grep(s/.rev.1.ebwt//, @file_name);;

#get fastq files
my $fastq_names_pointer=E_RNA::files_list($variables{rawdata_dir}, 'incrusive_file_names', 'fastq|fq');
my @fastq_names=@$fastq_names_pointer;

#species
our %species;
open my($IN), "<", $perl_dir.'/organisms.txt' or die;
while(<$IN>){
	chomp($_);
	unless($_=~/^#/){
		my($mir_organism, $mir_division, $mir_name, $mir_tree, $mir_NCBI_taxid)=split("\t", $_);
		$mir_name=~s/\s/_/g;
		$species{$mir_name}=$mir_organism;
	}
}
close($IN);

#################################################
#GUI interface
#create window
my $window = Gtk2::Window->new('toplevel');
$window->set_border_width(10);
$window->set_position('center');
$window->set_size_request('1000','700');
$window ->signal_connect( 'destroy' => \&exit );
$window->set_title('miRspring GUI');
	my $table=Gtk2::Table->new(12,6,TRUE);
	
		############frame: 
		my $frame=Gtk2::Frame->new();
			my $label=Gtk2::Label->new('Options of adaptor removal and sequence alignment');
		$frame->add($label);
	$table->attach_defaults($frame, 1,5,0,1);
		#frame: select fastq files
		$frame=Gtk2::Frame->new('Select fastq files');
		$frame->set_border_width(10);
			my $sub_table=items_selection(\@fastq_names);
		$frame->add($sub_table);
	$table->attach_defaults($frame, 0,6,1,5);
		#initiate sample name
		$frame=Gtk2::Frame->new('Enter sample name');
			my $entry_sample_name = Gtk2::Entry->new();
			$entry_sample_name->set_text('mirna');
		$frame->add($entry_sample_name);
	$table->attach_defaults($frame, 0,2,5,6);
		#initiate 3' adapter sequences
		$frame=Gtk2::Frame->new('Sequences of 3\' adapter');
			my $entry_3end_adapter = Gtk2::Entry->new();
			$entry_3end_adapter->set_text('AGATCGGAAGAGCACACGTCT');
		$frame->add($entry_3end_adapter);
	$table->attach_defaults($frame, 2,4,5,6);
		#frame: select genome sequences in fasta file
		$frame=Gtk2::Frame->new('Genome sequences in fasta format)');
			my $file_chooser_genome =Gtk2::FileChooserButton->new ('select a file' , 'open');
			$file_chooser_genome->set_filename($variables{ref_seq_dir});
		$frame->add($file_chooser_genome);
	$table->attach_defaults($frame, 4,6,5,6);
		#bowtie options
		$frame=Gtk2::Frame->new('Options of bowtie (default is bowtie1) for genome alignment');
			my $entry_bowtie_opt = Gtk2::Entry->new();
			$entry_bowtie_opt->set_text('-n 0 -m 5 -a --best --strata');
		$frame->add($entry_bowtie_opt);
	$table->attach_defaults($frame, 0,4,6,7);

	
		############frame: 
		$frame=Gtk2::Frame->new();
			$label=Gtk2::Label->new('Options of miRspring');
		$frame->add($label);
	$table->attach_defaults($frame, 1,5,7,8);
		#frame: select precursor miRNA file
		$frame=Gtk2::Frame->new('MirBase miRNA precursor file');
			my $file_chooser_precursor =Gtk2::FileChooserButton->new ('select a file' , 'open');
			$file_chooser_precursor->set_filename($variables{ref_seq_dir});
		$frame->add($file_chooser_precursor);
	$table->attach_defaults($frame, 0,2,8,9);
		#frame: select matued miRNA file
		$frame=Gtk2::Frame->new('MirBase matured miRNA file');
			my $file_chooser_mature =Gtk2::FileChooserButton->new ('select a file' , 'open');
			$file_chooser_mature->set_filename($variables{ref_seq_dir});
		$frame->add($file_chooser_mature);
	$table->attach_defaults($frame, 2,4,8,9);
		#frame: select gff file
		$frame=Gtk2::Frame->new('MirBase GFF file');
			my $file_chooser_gff =Gtk2::FileChooserButton->new ('select a file' , 'open');
			$file_chooser_gff->set_filename($variables{ref_seq_dir});
		$frame->add($file_chooser_gff);
	$table->attach_defaults($frame, 4,6,8,9);
		#miRNA species
		$frame=Gtk2::Frame->new('Species');
			my $cb_species = Gtk2::ComboBox->new_text;
			foreach (sort(keys %species)){
				$cb_species->append_text($_);
			}
			$cb_species->set_active(105);
		$frame->add($cb_species);
	$table->attach_defaults($frame, 0,2,9,10);
		#minimum length
		$frame=Gtk2::Frame->new('Minimum length of query reads');
			my $cb_min_len = Gtk2::ComboBox->new_text;
			foreach (16..30){
				$cb_min_len->append_text($_);
			}
			$cb_min_len->set_active(2);
		$frame->add($cb_min_len);
	$table->attach_defaults($frame, 2,4,9,10);
		#length of flanking sequences
		$frame=Gtk2::Frame->new('Length of flanking sequences');
			my $cb_flank_len = Gtk2::ComboBox->new_text;
			foreach (30..50){
				$cb_flank_len->append_text($_);
			}
			$cb_flank_len->set_active(5);
		$frame->add($cb_flank_len);
	$table->attach_defaults($frame, 4,6,9,10);
		#Mismatches
		$frame=Gtk2::Frame->new('Mismatches');
			my $cb_mismatches = Gtk2::ComboBox->new_text;
			foreach (0,1){
				$cb_mismatches->append_text($_);
			}
			$cb_mismatches->set_active(1);
		$frame->add($cb_mismatches);
	$table->attach_defaults($frame, 0,2,10,11);
	
		my $hbox = Gtk2::HBox->new( FALSE, 5 );
		$hbox->set_border_width(10);
			my $pbar = Gtk2::ProgressBar->new();
			$pbar->set_pulse_step(.1);
			$pbar->hide; #needs to be called after show_all
		$hbox->add($pbar);
			my $tbutton = Gtk2::Button->new_with_label('Run');
			my $sconnect;
			my $lconnect = $tbutton->signal_connect( clicked => \&launch);
		$hbox->add($tbutton);
			my $ebutton = Gtk2::Button->new_from_stock('gtk-quit');
			$ebutton->signal_connect( clicked =>\&exit );
		$hbox->add($ebutton);
	$table->attach_defaults($hbox, 0,6,11,12);
$window->add($table);

$window->show_all();
Gtk2->main;

#######################################
sub exit{
	$shash{'die'} = 1;
	foreach my $thread(threads->list()){
		$thread->join;
	}
	Gtk2->main_quit;
	return FALSE;
} 
#######################
sub launch{
	#print "$variables{mirspring_dir}\n", "$variables{mirspring_intermediate_script}\n", "$variables{mirspring_script}\n";
	$variables{mirspring_3end_adapter}=$entry_3end_adapter->get_text();
	$variables{mirspring_genome_file}=$file_chooser_genome->get_filename;
	$variables{mirspring_sample_name}=$entry_sample_name->get_text();
	$variables{mirspring_query_file}=$variables{mirspring_dir}.'/'.$variables{mirspring_sample_name}.'.fa';
	my @array=split("/", $variables{mirspring_genome_file});
	my $index_name=$array[-1];
	$index_name=~s/\.fa$|\.fasta$//;
	$variables{mirspring_index}=$variables{alignment_dir}.'/'.$index_name;
	$variables{mirspring_output_head}=$variables{mirspring_dir}.'/'.$variables{mirspring_sample_name}.'_'.$index_name;
	$variables{mirspring_sam_file}=$variables{mirspring_output_head}.'.sam';
	$variables{mirspring_bowtie_build_command}=join(" ", $variables{alignment_dir}.'/bowtie-build', 
											$variables{mirspring_genome_file}, $variables{mirspring_index}, 
										);
	$variables{mirspring_bowtie_command}=join(" ", $variables{alignment_dir}.'/bowtie', 
											$entry_bowtie_opt->get_text(), 
											$variables{mirspring_index}, 
											'-f', $variables{mirspring_query_file} ,
											'-S', $variables{mirspring_sam_file},  
										);
	$variables{mirspring_bam_file}=$variables{mirspring_output_head}.'.bam';
	$variables{mirspring_sorted_bam_file_head}=$variables{mirspring_output_head}.'.sorted';
	$variables{mirspring_sorted_bam_file}=$variables{mirspring_output_head}.'.sorted.bam';

	$variables{mirspring_min_query_len}=$cb_min_len->get_active_text();
	$variables{mirspring_mature_file}=$file_chooser_mature->get_filename;
	$variables{mirspring_precursor_file}=$file_chooser_precursor->get_filename;
	$variables{mirspring_species}=$species{ $cb_species->get_active_text() };
	$variables{mirspring_gff_file}=$file_chooser_gff->get_filename;
	$variables{mirspring_intermediate_result_file}=$variables{mirspring_output_head}.'_mirspring.txt';
	$variables{mirspring_result_file}=$variables{mirspring_output_head}.'_mirspring.html';
	
	$variables{mirspring_intermediate_command}=join(' ', 
				$variables{mirspring_intermediate_script},
				'-mm', $cb_mismatches->get_active_text(),
				'-ml', $variables{mirspring_min_query_len},
				'-flank', $cb_flank_len->get_active_text(),
				'-s', $variables{mirspring_species},
				'-bam', $variables{mirspring_sorted_bam_file},
				'-gff', $variables{mirspring_gff_file},
				'-mat', $variables{mirspring_mature_file},
				'-pre', $variables{mirspring_precursor_file},
				'-out', $variables{mirspring_intermediate_result_file},
			);
	$variables{mirspring_command}=join(' ', 
				$variables{mirspring_script}, 
				'-mm', $cb_mismatches->get_active_text(),
				'-ml', $cb_min_len->get_active_text(),
				'-flank', $cb_flank_len->get_active_text(),
				'-s', $variables{mirspring_species},
				'-in', $variables{mirspring_intermediate_result_file},
				'-gff', $variables{mirspring_gff_file},
				'-mat', $variables{mirspring_mature_file},
				'-pre', $variables{mirspring_precursor_file},
				'-out', $variables{mirspring_result_file},
			);
	my $thread_1;
	if ($shash{'go'}==0 ){
		#create 1 sleeping thread
		my $thread_1 = threads->new(\&work);
		
		$shash{'go'} = 1;
		$pbar->show;
		$tbutton->set_label('Running');
		$tbutton->signal_handler_block($lconnect);
		$sconnect = $tbutton->signal_connect( clicked => sub{ 	
			$shash{'go'} = 0;
			#$tbutton->set_label('Run');
			$tbutton->signal_handler_block ($sconnect);
			$tbutton->signal_handler_unblock ($lconnect);
		});

		Glib::Timeout->add (100, sub {
			if($shash{'go'} == 1){
				$pbar->set_text('Running!');
				$pbar->pulse;
				return TRUE;
			}
			else{	
				$thread_1->join();
				$pbar->set_text('OK! It is done!');
				$tbutton->set_label('Run');
				return FALSE;
			}
		});
	}


}

################################################## #######
sub work{
	$|++;

	while(1){
		return if $shash{'go'} == 0;
        
		if ( $shash{'go'} == 1 ){#2
			
			#current directory
			chdir $variables{mirspring_dir};
			
			#read fastq files and remove adapter sequences
			#get fastq files
			print "Read sequences from fastq files and remove adapter\n\n\n";
			my $fastq_files_pointer=E_RNA::files_list($variables{rawdata_dir}, 'incrusive_files', 'fastq|fq');
			my @fastq_files=@$fastq_files_pointer;
			my @sub_fastq_files;
			foreach my $f1(@selected_names){
				foreach my $f2(@fastq_files){
					push (@sub_fastq_files, $f2) if $f2=~/$f1/;
				}
			}
			read_fastq(\%variables, \@sub_fastq_files) unless -f $variables{mirspring_query_file}; #subroutine
			
			#alignment using bowtie
			print "Bowtie alignment\n\n\n";
			print " $variables{mirspring_bowtie_command}\n";
			system("$variables{mirspring_bowtie_build_command}") unless -f $variables{mirspring_index}.'.1.ebwt';
			system("$variables{mirspring_bowtie_command}");
			
			#bam file using samtools
			print "samtools transfering\n\n\n";
			system("samtools view -bS  $variables{mirspring_sam_file} -o $variables{mirspring_bam_file}");
			system("samtools sort $variables{mirspring_bam_file} $variables{mirspring_sorted_bam_file_head}");
			system("samtools index $variables{mirspring_sorted_bam_file}");
			
			#generate mirsrping files
			print "Run the first script of miRspring\n";
			E_RNA::refresh_log($var_file, 'mirspring_intermediate_command', "perl $variables{mirspring_intermediate_command}");
			system("perl $variables{mirspring_intermediate_command}");

			print "Run the second script of miRspring\n";
			E_RNA::refresh_log($var_file, 'mirspring_command', "perl $variables{mirspring_command}");
			system("perl $variables{mirspring_command}");

			$shash{'go'} = 0;   #turn off 
		}#2
		else{ sleep 1; }
	}
}

################
sub tree_store_str{
	my ($tree_store)=@_;
	my @items;
	$tree_store->foreach(sub{
			my($model, $path, $iter)=@_;
			my $value = $model->get ($iter, 0);
			push(@items, $value);
			return(FALSE);
		});
	my $str=join(',', @items);
	return($str);
}

#####################################
#used for items' selection
#return  @selected_names is fastq files
sub items_selection{
	my ($index_names_pointer)=@_;
	my @index_names=@$index_names_pointer;

	my $table=Gtk2::Table->new(4,8,TRUE);
	$table->set_sensitive (TRUE);
	
	#create a scrolled window that will host the treeview
	my $sw_1 = Gtk2::ScrolledWindow->new(undef, undef);
	$sw_1->set_shadow_type ('etched-out');
	$sw_1->set_policy ('automatic', 'automatic');
	$sw_1->set_border_width(5);
	
		#this is one of the provided base Gtk2::TreeModel classes.
		my $tree_store_1 = Gtk2::TreeStore->new(qw/Glib::String/);
		#fill it with arbitry data
		foreach (@index_names) {
			my $iter = $tree_store_1->append(undef); #the iter is a pointer in the treestore. We use to add data.
			$tree_store_1->set ($iter,0 => $_);
		}
	
		#this will create a treeview, specify $tree_store as its model
		my $tree_view_1 = Gtk2::TreeView->new($tree_store_1);
		#create a Gtk2::TreeViewColumn to add to $tree_view
		my $tree_column_1 = Gtk2::TreeViewColumn->new();
		$tree_column_1->set_title ("Sort candidates");
			#create a renderer that will be used to display info in the model
			my $renderer_1 = Gtk2::CellRendererText->new;
		#add this renderer to $tree_column. This works like a Gtk2::Hbox. so you can add more than one renderer to $tree_column			
		$tree_column_1->pack_start ($renderer_1, FALSE);
		$tree_column_1->add_attribute($renderer_1, text => 0);
		$tree_column_1->set_sort_column_id(0);# Allow sorting on the column
	
		#add $tree_column to the treeview
		$tree_view_1->append_column ($tree_column_1);
		$tree_view_1->set_search_column(0);# make it searchable
		$tree_view_1->set_reorderable(TRUE);# Allow drag and drop reordering of rows
	
	$sw_1->add($tree_view_1);
	$table->attach_defaults($sw_1, 0,3,0,6);

	#create a scrolled window that will host the treeview
	my $sw_2 = Gtk2::ScrolledWindow->new (undef, undef);
	$sw_2->set_shadow_type ('etched-out');
	$sw_2->set_policy ('automatic', 'automatic');
	$sw_2->set_border_width(5);
	
		#this is one of the provided base Gtk2::TreeModel classes.
		my $tree_store_2 = Gtk2::TreeStore->new(qw/Glib::String/);
		my $tree_view_2 = Gtk2::TreeView->new($tree_store_2);

		#create a Gtk2::TreeViewColumn to add to $tree_view
		my $tree_column_2 = Gtk2::TreeViewColumn->new();
		$tree_column_2->set_title ("Sort selected items");
			
			#create a renderer that will be used to display info in the model
			my $renderer_2 = Gtk2::CellRendererText->new;
		#add this renderer to $tree_column. 
		$tree_column_2->pack_start ($renderer_2, FALSE);
		
		 # set the cell "text" attribute to column 0 - retrieve text from that column in treestore 
		$tree_column_2->add_attribute($renderer_2, text => 0);
		$tree_column_2->set_sort_column_id(0);# Allow sorting on the column
	
		#add $tree_column to the treeview
		$tree_view_2->append_column ($tree_column_2);
		$tree_view_2->set_search_column(0);# make it searchable
		$tree_view_2->set_reorderable(TRUE);# Allow drag and drop reordering of rows

	$sw_2->add($tree_view_2);
	$table->attach_defaults($sw_2, 4,7,0,6);

	#Widget of table: button
	my $button_1=Gtk2::Button->new('>');        ######>> button: one transfered
	$button_1->signal_connect(clicked=>sub{
		my $tree_selection_1=$tree_view_1->get_selection;
		my($tree_store_1, $iter_1)=$tree_selection_1->get_selected;
		my $value=$tree_store_1->get($iter_1, 0);
		return unless $iter_1;
		$tree_store_1->remove($iter_1);
		my $iter=$tree_store_2->append(undef);  #add the value into the second tree 
		$tree_store_2->set($iter,0=> $value);
		push(@selected_names, $value);
	} );
	$table->attach_defaults($button_1, 3,4,1,2);
	my $button_2=Gtk2::Button->new('>>');        ###### >> button: all transfered
	$button_2->signal_connect(clicked=>sub{
		$tree_store_1->clear;
		foreach (@index_names) {
			my $iter = $tree_store_2->append(undef); #the iter is a pointer in the treestore. We use to add data.
			$tree_store_2->set ($iter,0 => $_);
		}
		@selected_names=@index_names;
	} );
	$table->attach_defaults($button_2, 3,4,2,3);
	my $button_3=Gtk2::Button->new('<');      ######<< button
	$button_3->signal_connect(clicked=>sub{
		my $tree_selection_2=$tree_view_2->get_selection;
		my($tree_store_2, $iter_2)=$tree_selection_2->get_selected;
		my $value=$tree_store_2->get($iter_2, 0);
		return unless $iter_2;
		$tree_store_2->remove($iter_2);
		my $iter=$tree_store_1->append(undef);  #add the value into the second tree 
		$tree_store_1->set($iter,0=> $value);
		foreach(my $i=0; $i<@selected_names; $i++){
			delete $selected_names[$i] if $selected_names[$i] eq $value;
		}
		@selected_names=grep($_, @selected_names);
	} );
	$table->attach_defaults($button_3, 3,4,3,4);
	my $button_4=Gtk2::Button->new('<<');      ######<< button
	$button_4->signal_connect(clicked=>sub{
		$tree_store_2->clear;
		foreach (@index_names) {
			my $iter = $tree_store_1->append(undef); #the iter is a pointer in the treestore. We use to add data.
			$tree_store_1->set ($iter,0 => $_);
		}
		undef @selected_names;
	} );
	$table->attach_defaults($button_4, 3,4,4,5);
	my $button_5=Gtk2::Button->new('Up');        ######Up button
	$button_5->signal_connect(clicked=>sub{
		my $tree_selection=$tree_view_2->get_selection;
		my($tree_store, $iter)=$tree_selection->get_selected;
		my $path = $tree_store->get_path($iter);
		$path->prev;
		my $prev_iter = $tree_store->get_iter($path);
		$tree_store->move_before($iter,$prev_iter );
		
		my $value=$tree_store->get($iter, 0);
		my $pre_value=$tree_store->get($prev_iter, 0);
		foreach(my $i=0; $i<@selected_names; $i++){
			if ( ($selected_names[$i] eq $value) and $i>0){
				$selected_names[$i-1]=$value;
				$selected_names[$i]=$pre_value;
			}
		}
		@selected_names=grep($_, @selected_names);
	} );
	$table->attach_defaults($button_5, 7,8,1,2);
	my $button_6=Gtk2::Button->new('Down');   ######Down button
	$button_6->signal_connect(clicked=>sub{
		my $tree_selection=$tree_view_2->get_selection;
		my($tree_store, $iter)=$tree_selection->get_selected;
		my $next_iter = $tree_store->iter_next($iter);
		$tree_store->move_before($next_iter,$iter );
				
		my $value=$tree_store->get($iter, 0);
		my $next_value=$tree_store->get($next_iter, 0);
		foreach(my $i=0; $i<@selected_names; $i++){
			if ( ($selected_names[$i] eq $value) and $i<@selected_names-1){
				$selected_names[$i+1]=$value;
				$selected_names[$i]=$next_value;
			}
		}
		@selected_names=grep($_, @selected_names);
	} );
	$table->attach_defaults($button_6, 7,8,4,5);

	return ($table);
}
#########################
sub read_fastq{
	my($variables_pointer, $fastq_pointer)=@_;
	my %variables=%$variables_pointer;
	my @fastq=@$fastq_pointer;
	
	my $num=1;
	open my($OUT), ">", $variables{mirspring_query_file} or die; 
	foreach my $fastq_file(@fastq){#2
		my $in_obj = Bio::SeqIO->new(-file => $fastq_file, -format => 'fastq');
		while (my $seq_obj = $in_obj->next_seq() ) {#3
			my $displayid=$variables{mirspring_sample_name}.'_'.$num;
			my $read_seq=$seq_obj->seq();
			my $read_seq_trunc=end3_seq_truncation($variables{mirspring_3end_adapter}, $read_seq);
			print $OUT ">$displayid\n", "$read_seq_trunc\n" if $read_seq_trunc cmp 'NA' and length($read_seq_trunc)>=$variables{mirspring_min_query_len};
			$num++;
		}#3
	}#2
	close($OUT);
	
}

##################################################################
#truncate sequence of adapter from read
sub end3_seq_truncation{#1
	my $adapter=$_[0];
	my $seq=$_[1];
	my $len_adapter=length($adapter);
	my $len_seq=length($seq);
	
	my $seq_trunc=$seq;
	if (length($seq_trunc)<1){#2  all adapter sequence
			$seq_trunc="NA"; #no insert
	} #2
	else{	#2
		#sequence is longer
		if ($len_seq>=$len_adapter){#3
			if ($seq_trunc=~/$adapter/){#4
				$seq_trunc=~s/$adapter.*//;
			}
			else{
				for (my $i=$len_adapter-1; $i>=4; $i--){
					my $adapter_index=substr($adapter, 0, $i);
					last if $seq_trunc=~s/$adapter_index$//;
				}
			}#4
		}#3
		#adapter sequence is longer
		else{#3
			for (my $i=$len_seq-1; $i>=4; $i--){#4
				my $adapter_index=substr($adapter, 0, $i);
				last if $seq_trunc=~s/$adapter_index$//;
			}#4
		}#3
	}#2

	return($seq_trunc);
}#1
