#!/usr/bin/perl -w
use warnings;
use strict;
use Cwd;
use List::Util;
use threads 'exit' => 'threads_only';
use threads::shared;
use Glib qw/TRUE FALSE/;
use Gtk2 qw/-init -threads-init/;
use Gtk2::Pango;
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

#initiate environmental variables
print "\n\nEnivironmental variables:\n";
$ENV{'PATH'} .= ":$variables{alignment_dir}:$variables{ref_seq_dir}:$variables{mirdeep_dir}:";
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
our %species=('Human'=>'hsa', 'Chimp'=>'ptr', 'Orangutan'=>'na', 'Rhesus'=>'na', 'Marmoset'=>'na', 'Mouse'=>'mmu', 'Rat'=>'rno', 
		'Guinea_Pig'=>'na', 'Cat'=>'lca', 'Dog'=>'cfa', 'Horse'=>'eca', 'Cow'=>'bta', 'Opossum'=>'na', 'Platypus'=>'na', 'Chicken'=>'gga',
		'Zebra_finch'=>'na', 'Lizard'=>'na', 'X_tropicalis'=>'xtr', 'Zebrafish'=>'tni', 'Tetraodon'=>'tni', 'Fugu'=>'fru', 'Stickleback'=>'na',
		'Medaka'=>'na', 'Lamprey'=>'na', 'Lancelet'=>'bfl', 'C_intestinalis'=>'cin', 'S_purpuratus'=>'spu', 'C_elegans'=>'cel',
		'C_brenneri'=>'na', 'C_briggsae'=>'cbr', 'C_remanei'=>'na', 'C_japonica'=>'sja', 'P_pacificus'=>'na', 'D_melanogaster'=>'dme',
		'D_simulans'=>'dsi', 'D_sechellia'=>'dse', 'D_yakuba'=>'dya', 'D_erecta'=>'der', 'D_ananassae'=>'dan', 'D_pseudoobscura'=>'dps',
		'D_persimilis'=>'dpe', 'D_virilis'=>'dvi', 'D_mojavensis'=>'dmo', 'D_grimshawi'=>'dgr', 'A_gambiae'=>'aga', 'A_mellifera'=>'ame',
		'S_cerevisiae'=>'na', 'worm'=>'cel', 'All'=>'All',
		);

#################################################
#GUI interface
#create window
my $window = Gtk2::Window->new('toplevel');
$window->set_border_width(10);
$window->set_position('center');
$window->set_size_request('1000','700');
$window ->signal_connect( 'destroy' => \&exit );
$window->set_title('miRDeep2 GUI');
	my $table=Gtk2::Table->new(13,6,TRUE);
	
		#frame: select fastq files
		my $frame=Gtk2::Frame->new('Select fastq files');
		$frame->set_border_width(10);
			my $sub_table=items_selection(\@fastq_names);
		$frame->add($sub_table);
	$table->attach_defaults($frame, 0,6,0,4);
		#initiate 3' adapter sequences
		$frame=Gtk2::Frame->new('Sequences of 3\' adapter');
			my $entry_3end_adapter = Gtk2::Entry->new();
			$entry_3end_adapter->set_text('AGATCGGAAGAGCACACGTCT');
		$frame->add($entry_3end_adapter);
	$table->attach_defaults($frame, 0,4,4,5);
		#frame: select genome sequences in fasta file
		$frame=Gtk2::Frame->new('Genome sequences in fasta format)');
			my $file_chooser_genome =Gtk2::FileChooserButton->new ('select a file' , 'open');
			$file_chooser_genome->set_filename($variables{ref_seq_dir});
		$frame->add($file_chooser_genome);
	$table->attach_defaults($frame, 4,6,4,5);
		#frame: select prematured miRNA fasta file
		$frame=Gtk2::Frame->new('miRNA precursor sequences from miRBase in fasta format');
			my $file_chooser_prematured =Gtk2::FileChooserButton->new ('select a file' , 'open');
			$file_chooser_prematured->set_filename($variables{ref_seq_dir});
		$frame->add($file_chooser_prematured);
	$table->attach_defaults($frame, 0,2,5,6);
		#frame: select matured miRNA fasta file
		$frame=Gtk2::Frame->new('Matured miRNA sequences from miRBase in fasta format');
			my $file_chooser_matured =Gtk2::FileChooserButton->new ('select a file' , 'open');
			$file_chooser_matured->set_filename($variables{ref_seq_dir});
		$frame->add($file_chooser_matured);
	$table->attach_defaults($frame, 2,4,5,6);
		#miRNA species
		$frame=Gtk2::Frame->new('miRNA species');
			my $cb_species = Gtk2::ComboBox->new_text;
			foreach (sort(keys %species)){
				$cb_species->append_text($_);
			}
			$cb_species->set_active(0);
		$frame->add($cb_species);
	$table->attach_defaults($frame, 4,6,5,6);
	
		###############frame: Other options
		$frame=Gtk2::Frame->new();
			my $label=Gtk2::Label->new('Other options');
		$frame->add($label);
	$table->attach_defaults($frame, 1,5,6,7);
		#frame: mirdeep directory
		$frame=Gtk2::Frame->new('MiRDeep directory');
			my $dir_chooser=Gtk2::Button->new ('Select another folder');
			$dir_chooser->signal_connect(clicked=>sub{
				my $file_chooser =Gtk2::FileChooserDialog->new ('Folder selection',  undef, 'select-folder',   'gtk-cancel' => 'cancel', 'gtk-ok' => 'ok');
					my $preview_widget = Gtk2::Label->new ('wheeee');
					$preview_widget->set_line_wrap (TRUE);
					$preview_widget->set_size_request (150, -1);
					$file_chooser->set (preview_widget => $preview_widget,  preview_widget_active => TRUE);
					$file_chooser->signal_connect (selection_changed => sub {
						my $filename = $file_chooser->get_preview_filename;
						my $active = defined $filename && not -d $filename;
						if ($active) {
							my $size = sprintf '%.1fK', (-s $filename) / 1024;
							my $desc = `file '$filename'`;
							$desc =~ s/^$filename:\s*//;
							$preview_widget->set_text ("$size\n$desc");
						}
						$file_chooser->set (preview_widget_active => $active);
					});
					if ('ok' eq $file_chooser->run) {	
						$variables{mirdeep_dir}=$file_chooser->get_filename;
					}
					$file_chooser->destroy;
					$dir_chooser->set_label($variables{mirdeep_dir});
				});
			$frame->add($dir_chooser);
		$table->attach_defaults($frame, 0,2,7,8);
		#initiate sample name
		$frame=Gtk2::Frame->new('Enter sample name');
			my $entry_fastq_name = Gtk2::Entry->new();
			$entry_fastq_name->set_text('mirna');
		$frame->add($entry_fastq_name);
	$table->attach_defaults($frame, 2,4,7,8);
		#initiate prefix
		$frame=Gtk2::Frame->new('Three-letter prefix for reads');
			my $entry_prefix = Gtk2::Entry->new();
			$entry_prefix->set_text('seq');
		$frame->add($entry_prefix);
	$table->attach_defaults($frame, 4,6,7,8);
		#-l discarded length
		$frame=Gtk2::Frame->new('Minimum discarded length of reads');
			my $cb_discard_len = Gtk2::ComboBox->new_text;
			foreach (18..20){
				$cb_discard_len->append_text($_);
			}
			$cb_discard_len->set_active(0);
		$frame->add($cb_discard_len);
	$table->attach_defaults($frame, 0,2,8,9);
		#-g 50000
		$frame=Gtk2::Frame->new('Maxium precursors in automatic excision gearing');
			my $cb_max_precursors = Gtk2::ComboBox->new_text;
			foreach (5e4,1e5,-1){
				$cb_max_precursors->append_text($_);
			}
			$cb_max_precursors->set_active(0);
		$frame->add($cb_max_precursors);
	$table->attach_defaults($frame, 2,4,8,9);
		#-r number of reads aligned to reference
		$frame=Gtk2::Frame->new('Minimum number of reads aligned to references');
			my $cb_min_aligned = Gtk2::ComboBox->new_text;
			foreach (1..10){
				$cb_min_aligned->append_text($_);
			}
			$cb_min_aligned->set_active(0);
		$frame->add($cb_min_aligned);
	$table->attach_defaults($frame, 4,6,8,9);
		#miRNA mismatch
		$frame=Gtk2::Frame->new('Mismatch allowed with miRDeep');
			my $cb_mismatch = Gtk2::ComboBox->new_text;
			foreach (1..3){
				$cb_mismatch->append_text($_);
			}
			$cb_mismatch->set_active(0);
		$frame->add($cb_mismatch);
	$table->attach_defaults($frame, 0,2,9,10);
		#number of nucleotides upstream
		$frame=Gtk2::Frame->new('Number of nucleotides upstream of matured miRNA');
			my $cb_upstream = Gtk2::ComboBox->new_text;
			foreach (1..6){
				$cb_upstream->append_text($_);
			}
			$cb_upstream->set_active(1);
		$frame->add($cb_upstream);
	$table->attach_defaults($frame, 2,4,9,10);
		#number of nucleotides downstream
		$frame=Gtk2::Frame->new('Number of nucleotides downstream of matured miRNA');
			my $cb_downstream = Gtk2::ComboBox->new_text;
			foreach (4..8){
				$cb_downstream->append_text($_);
			}
			$cb_downstream->set_active(1);
		$frame->add($cb_downstream);
	$table->attach_defaults($frame, 4,6,9,10);
		#-b 0
		$frame=Gtk2::Frame->new('Minimum score for prediction of novel miRNAs');
			my $cb_min_prediction = Gtk2::ComboBox->new_text;
			foreach (0..10){
				$cb_min_prediction->append_text($_);
			}
			$cb_min_prediction->set_active(0);
		$frame->add($cb_min_prediction);
	$table->attach_defaults($frame, 0,2,10,11);
		#
		$frame=Gtk2::Frame->new('Both precursor and matured miRNA mapping');
			my $check_button_2 = Gtk2::CheckButton->new ("miRNA precursor mapping allowed");
			$check_button_2->set_active (FALSE);
		$frame->add($check_button_2);
	$table->attach_defaults($frame, 2,4,10,11);
		#
		$frame=Gtk2::Frame->new('Multiple mapping');
			my $check_button_3 = Gtk2::CheckButton->new ("discard all multiple mapping");
			$check_button_3->set_active (FALSE);
		$frame->add($check_button_3);
	$table->attach_defaults($frame, 4,6,10,11);
		#mismatch allowed
		$frame=Gtk2::Frame->new('Mismatch in the seed in alignment with bowtie');
			my $check_button_1 = Gtk2::CheckButton->new ("One mismatch in the seed allowed");
			$check_button_1->set_active (FALSE);
		$frame->add($check_button_1);
	$table->attach_defaults($frame, 0,2,11,12);
	
		#rnafold analysis
		$frame=Gtk2::Frame->new('RNAfold analysis');
			my $check_button_RNAfold = Gtk2::CheckButton->new ("Enable RNAfold analysis");
			$check_button_RNAfold->set_active (TRUE);
		$frame->add($check_button_RNAfold);
	$table->attach_defaults($frame, 2,4,11,12);
	
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
	$table->attach_defaults($hbox, 0,6,12,13);
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
	
	$variables{mirdeep_sample_name}=$entry_fastq_name->get_text();
	#output of mapper
	$variables{mirdeep_processed_file}= $variables{mirdeep_dir}.'/'.$variables{mirdeep_sample_name}.'_collapsed.fa';
	$variables{mirdeep_arf_file}= $variables{mirdeep_dir}.'/'.$variables{mirdeep_sample_name}.'_collapsed_vs_genome.arf';
	$variables{mirdeep_prefix}=$entry_prefix->get_text();
	$variables{mirdeep_3end_adapter}=$entry_3end_adapter->get_text();
	$variables{mirdeep_discard_len}=$cb_discard_len->get_active_text();
	$variables{mirdeep_min_aligned}=$cb_min_aligned->get_active_text();
	$variables{mirdeep_prematured_fasta_file}=$file_chooser_prematured->get_filename ;
	$variables{mirdeep_matured_fasta_file}=$file_chooser_matured->get_filename ;
	#subroutine genome fasta file checking
	$variables{mirdeep_raw_genome_fasta_file}=$file_chooser_genome->get_filename;
	my @array=split("/", $variables{mirdeep_raw_genome_fasta_file});
	my $fasta_name=$array[-1];
	$variables{mirdeep_genome_fasta_file}=$variables{mirdeep_dir}.'/'.$fasta_name;
	$variables{mirdeep_index}=$fasta_name;
	$variables{mirdeep_index}=~s/\.fa$|\.fasta$//;
	$variables{mirdeep_index}=$variables{mirdeep_dir}.'/'.$variables{mirdeep_index};

	#subroutine: get miRNA of a certain specie for miRBase fasta files
	$variables{mirdeep_species}=$species{ $cb_species->get_active_text() };
	$variables{mirdeep_specie_fasta_files}=specie_miRNA($variables{mirdeep_prematured_fasta_file}, $variables{mirdeep_matured_fasta_file}, $variables{mirdeep_species}, $variables{mirdeep_dir});
	
	$variables{mapper_options}=join (" ", '-e', '-h', '-j', '-m', '-n',
									'-g', $variables{mirdeep_prefix},
									'-l', $variables{mirdeep_discard_len},
									'-p', $variables{mirdeep_index},
									'-r', $variables{mirdeep_min_aligned},
									'-s', $variables{mirdeep_processed_file},
									'-t', $variables{mirdeep_arf_file},
									);
	$variables{mapper_options}=join(" ", $variables{mapper_options}, '-k', $variables{mirdeep_3end_adapter}) if length($variables{mirdeep_3end_adapter})>3;
	$variables{mapper_options}=join(" ", $variables{mapper_options}, '-q') if $check_button_1->get_active;
	
	$variables{quantifier_options}=join(" ", '-p', $variables{mirdeep_prematured_fasta_file},
								'-m', $variables{mirdeep_matured_fasta_file},
								'-r',  $variables{mirdeep_processed_file},
								'-g', $cb_mismatch->get_active_text(),
								'-e', $cb_upstream->get_active_text(),
								'-f', $cb_downstream->get_active_text(),
							);
	$variables{quantifier_options}=join(" ", $variables{quantifier_options}, '-t', $variables{mirdeep_species}) unless $variables{mirdeep_species} eq 'All';
	$variables{quantifier_options}=join(" ", $variables{quantifier_options}, '-k') if $check_button_2->get_active;
	$variables{quantifier_options}=join(" ", $variables{quantifier_options}, '-U') if $check_button_3->get_active;
	
	$variables{mirdeep_options}=join(" ", $variables{mirdeep_processed_file}, 
								$variables{mirdeep_genome_fasta_file}, 
								$variables{mirdeep_arf_file},
								$variables{mirdeep_specie_fasta_files},
								'-b', $cb_min_prediction->get_active_text(),
								'-g', $cb_max_precursors->get_active_text(),
								 '-t', $variables{mirdeep_species},
							);
	$variables{mirdeep_options}=join(" ", $variables{mirdeep_options}, 'c') unless $check_button_RNAfold->get_active;
	print "##mirdeep output: $variables{quantifier_options}\n" ;

	
	
	my $thread_1;
	if ($shash{'go'}==0 and @selected_names>0){
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
			chdir $variables{mirdeep_dir};
			
			#generate genome fasta and bowtie index
			print "Generate genome fasta file and build bowtie index";
			genome_fasta_checking($variables{mirdeep_raw_genome_fasta_file}, $variables{mirdeep_genome_fasta_file}, $variables{alignment_dir}, $variables{mirdeep_index});
			
			#fastq files
			my $fastq_files_pointer=E_RNA::files_list($variables{rawdata_dir}, 'incrusive_files', 'fastq|fq');
			my @fastq_files=@$fastq_files_pointer;
			my @sub_fastq_files;
			foreach my $f1(@selected_names){
				foreach my $f2(@fastq_files){
					push (@sub_fastq_files, $f2) if $f2=~/$f1/;
				}
			}
			my $sub_fastq_files_str=join(' ', @sub_fastq_files);
			my $combined_fastq_name=$variables{mirdeep_sample_name}.'.fastq';
			my $combined_fastq_file=$variables{mirdeep_dir}.'/'.$variables{mirdeep_sample_name}.'.fastq';
			system(" cat $sub_fastq_files_str > $combined_fastq_file ");
			# run mapper
			E_RNA::refresh_log($var_file, "mirdeep_mapper", "mapper.pl $combined_fastq_name $variables{mapper_options}");
			system("mapper.pl $combined_fastq_name $variables{mapper_options} ");
			#run quantifier
			system("quantifier.pl $variables{quantifier_options} ");
			E_RNA::refresh_log($var_file, "mirdeep_quantifier", "quantifier.pl $variables{quantifier_options} ");
			#run miRDeep2.pl
			system("miRDeep2.pl $variables{mirdeep_options} ");
			E_RNA::refresh_log($var_file, "mirdeep_mirdeep2", "miRDeep2.pl $variables{mirdeep_options} 2>report.log");
			last if $shash{'go'} == 0;
			$shash{'go'} = 0;   #turn off 
		}#2
		else{ sleep 1; }
	}
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

###############################################################
#used for items selection
sub multiple_items_selection{
	my ($index_names_pointer)=@_;
	my @index_names=@$index_names_pointer;
	
	my $table=Gtk2::Table->new(4,8,TRUE);
	$table->set_sensitive (TRUE);
	
	#create a scrolled window that will host the treeview
	my $sw_1 = Gtk2::ScrolledWindow->new (undef, undef);
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
		$tree_column_1->set_title ("Sort sample names");
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
	$table->attach_defaults($sw_1, 4,7,0,6);

	  #group A
	#create a scrolled window that will host the treeview 
	my $sw_2 = Gtk2::ScrolledWindow->new (undef, undef);
	$sw_2->set_shadow_type ('etched-out');
	$sw_2->set_policy ('automatic', 'automatic');
	$sw_2->set_border_width(5);
		my $tree_store_2 = Gtk2::TreeStore->new(qw/Glib::String/);
		my $tree_view_2 = Gtk2::TreeView->new($tree_store_2);
		my $tree_column_2 = Gtk2::TreeViewColumn->new();
		$tree_column_2->set_title ("Sort selected group A");
			my $renderer_2 = Gtk2::CellRendererText->new;
		$tree_column_2->pack_start ($renderer_2, FALSE);
		$tree_column_2->add_attribute($renderer_2, text => 0);
		$tree_column_2->set_sort_column_id(0);# Allow sorting on the column
		$tree_view_2->append_column ($tree_column_2);
		$tree_view_2->set_search_column(0);# make it searchable
		$tree_view_2->set_reorderable(TRUE);# Allow drag and drop reordering of rows
	$sw_2->add($tree_view_2);
	$table->attach_defaults($sw_2, 0,3,0,6);

	#Widget of table: buttons of group A
	my $button_1=Gtk2::Button->new('>');        ######>> button: one transfered
	$button_1->signal_connect(clicked=>sub{
		my $tree_selection=$tree_view_2->get_selection;
		my($tree_store, $iter_1)=$tree_selection->get_selected;
		my $value=$tree_store->get($iter_1, 0);
		return unless $iter_1;
		$tree_store->remove($iter_1);
		my $iter=$tree_store_1->append(undef);  #add the value into the second tree 
		$tree_store_1->set($iter,0=> $value);
		$selected_A=tree_store_str($tree_store_2);#subroutine
	} );
	$table->attach_defaults($button_1, 3,4,1,2);

	my $button_3=Gtk2::Button->new('<');      ######<< button
	$button_3->signal_connect(clicked=>sub{
		my $tree_selection=$tree_view_1->get_selection;
		my($tree_store, $iter_1)=$tree_selection->get_selected;
		my $value=$tree_store->get($iter_1, 0);
		return unless $iter_1;
		$tree_store_1->remove($iter_1);
		my $iter=$tree_store_2->append(undef);  #add the value into the second tree 
		$tree_store_2->set($iter,0=> $value);
		$selected_A=tree_store_str($tree_store_2);#subroutine
	} );
	$table->attach_defaults($button_3, 3,4,4,5);
	
	  #group B
	#create a scrolled window that will host the treeview 
	my $sw_3 = Gtk2::ScrolledWindow->new (undef, undef);
	$sw_3->set_shadow_type ('etched-out');
	$sw_3->set_policy ('automatic', 'automatic');
	$sw_3->set_border_width(5);
		my $tree_store_3 = Gtk2::TreeStore->new(qw/Glib::String/);
		my $tree_view_3 = Gtk2::TreeView->new($tree_store_3);
		my $tree_column_3 = Gtk2::TreeViewColumn->new();
		$tree_column_3->set_title ("Sort selected group B");
			my $renderer_3 = Gtk2::CellRendererText->new;
		$tree_column_3->pack_start ($renderer_3, FALSE);
		$tree_column_3->add_attribute($renderer_3, text => 0);
		$tree_column_3->set_sort_column_id(0);# Allow sorting on the column
		$tree_view_3->append_column ($tree_column_3);
		$tree_view_3->set_search_column(0);# make it searchable
		$tree_view_3->set_reorderable(TRUE);# Allow drag and drop reordering of rows
	$sw_3->add($tree_view_3);
	$table->attach_defaults($sw_3, 8,11,0,6);
	
	#Widget of table: button
	my $button_4=Gtk2::Button->new('<');      ######< button
	$button_4->signal_connect(clicked=>sub{
		my $tree_selection=$tree_view_3->get_selection;
		my($tree_store, $iter_1)=$tree_selection->get_selected;
		my $value=$tree_store->get($iter_1, 0);
		return unless $iter_1;
		$tree_store_3->remove($iter_1);
		my $iter=$tree_store_1->append(undef);  #add the value into the second tree 
		$tree_store_1->set($iter,0=> $value);
		$selected_B=tree_store_str($tree_store_3); #subroutine
	} );
	$table->attach_defaults($button_4, 7,8,1,2);

	my $button_6=Gtk2::Button->new('>');        ######> button: one transfered
	$button_6->signal_connect(clicked=>sub{
		my $tree_selection=$tree_view_1->get_selection;
		my($model, $iter_1)=$tree_selection->get_selected;
		my $value=$model->get($iter_1, 0);
		return unless $iter_1;
		$tree_store_1->remove($iter_1);
		my $iter=$tree_store_3->append(undef);  #add the value into the second tree 
		$tree_store_3->set($iter,0=> $value);
		$selected_B=tree_store_str($tree_store_3); #subroutine
	} );
	$table->attach_defaults($button_6, 7,8,4,5);

	$table->show_all();
	return ($table);
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

###################################
#generate specie miRNA fasta file based on miRBase fasta files
sub specie_miRNA{
	my($miRNA_precursor_fasta_file, $matured_miRNA_fasta_file, $specie, $dir)=@_;
	
	my $specie_matured_fasta_file=$dir.'/specie_matured_miRNA.fa';
	my $other_matured_fasta_file=$dir.'/other_matured_miRNA.fa';
	my $in1_obj=Bio::SeqIO->new(-file=> $matured_miRNA_fasta_file, -format=>"fasta");
	open my($out1), ">", $specie_matured_fasta_file or die;
	open my($out2), ">", $other_matured_fasta_file or die;
	while (my $seq_obj=$in1_obj->next_seq){
		my $seq=$seq_obj->seq();
		my @array=split(' ', $seq_obj->display_id());
		my $display_id=$array[0];
		if ($display_id=~/$specie/){
			print $out1 ">$display_id\n", "$seq\n";
		}
		else{
			print $out2 ">$display_id\n", "$seq\n";
		}
	}
	close($out1);
	close($out2);
	
	my $specie_precursor_fasta_file=$dir.'/specie_precursor_miRNA.fa';
	my $in2_obj=Bio::SeqIO->new(-file=> $miRNA_precursor_fasta_file, -format=>"fasta");
	open my($out3), ">", $specie_precursor_fasta_file or die;
	while (my $seq_obj=$in2_obj->next_seq){
		my $seq=$seq_obj->seq();
		my @array=split(" ", $seq_obj->display_id());
		my $display_id=$array[0];
		if ($display_id=~/$specie/){
			#print "$display_id\n";
			print $out3 ">$display_id\n", "$seq\n";
		}
	}
	close($out3);
	
	my $out=join(" ", $specie_matured_fasta_file, $other_matured_fasta_file, $specie_precursor_fasta_file);
	return($out);
}

##########################################
#
sub genome_fasta_checking{
	my ($fasta_file, $checked_fasta_file, $alignment_dir, $index_name)=@_;
	
	unless (-f $checked_fasta_file){#2
		print "Copy genome fasta file into miRDeep output directory\n";
		open my($OUT), ">", $checked_fasta_file or die;
		my $in_obj=Bio::SeqIO->new(-file=> $fasta_file, -format=>"fasta");
		while (my $seq_obj=$in_obj->next_seq){
			my $display_id=$seq_obj->display_id();
			my $seq=$seq_obj->seq();
			$seq=~tr/a-z/A-Z/; #lower case 
			$seq=~s/[^A|T|G|C]/N/g; #transfer non A,T,G,C
			print $OUT ">$display_id\n", "$seq\n";
		}
		close($OUT);
	}#2
	
	unless(-f $index_name.'.1.ebwt' and -f $index_name.'.rev.1.ebwt'){
		print "Build bowtie index\n";
		system("$alignment_dir/bowtie-build $checked_fasta_file $index_name");
	}
	
}

