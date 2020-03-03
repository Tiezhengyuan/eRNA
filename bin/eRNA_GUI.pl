#! /usr/bin/perl -w     
use warnings;
use strict;
use Cwd;
use Gtk2 qw/-init -threads-init/;
#use Glib qw/TRUE FALSE/; 
use Gtk2::Pango;
use threads; 
use threads::shared;
use List::Util;
use List::MoreUtils;




#get subroutines
require "func_gui.pm";
require "func_rna.pm";
require "func_basic.pm";
require "func_common.pm";
require "func_ncrna.pm";

#############################3
#global variables
#species
our %species=('Human'=>'hsa', 'Mouse'=>'mmu', 'Rat'=>'rno', 'Dog'=>'cfa', 'Horse'=>'eca', 'Cow'=>'bta', 'Chicken'=>'gga',
							'Zebrafish'=>'tni', 'C_elegans'=>'cel', 'D_melanogaster'=>'dme', 'worm'=>'cel', 'maize'=>'zma', 'All'=>'All', );
		
#get E_RNA variables 
our %variables=( 'software_aligner'=>'bowtie2', 'software_spliced_aligner'=>'tophat2', 'software_assembler'=>'cufflinks',
				'mRNA_genome_mapping'=>'no', 'mRNA_sequencing_control'=>'no',
				'mRNA_transcripts_assembling'=>'no', 'mRNA_cuffmerge'=>'no', 'mRNA_cuffdiff'=>'no',
				'raw_file_format'=>'FASTQ', 'quality_filter'=>13, 'QC_compression'=>1000, 'region_len'=>10000, 
				#'adapter_3'=>'AGATCGGAAGAGCACACGTCT', 'adapter_5'=>'GTTCAGAGTTCTACAGTCCGACGATC', #NEB
				'adapter_5'=>'GTTCAGAGTTCTACAGTCCGACGATC', 'adapter_3'=>'TGGAATTCTCGGGTGCCAAGG', #Illumina
				'match_len'=>8, 'mismatch_allowed'=>'yes', 'mismatch_len'=>12,
				'ncRNA_read_rawdata'=>'no', 'ncRNA_adapter_removal'=>'no',
				'ncRNA_seperate_alignment'=>'no', 'ncRNA_iterative_alignment'=>'no', 'ncRNA_counting'=>'no',
				'index_seperate'=>'no', 'index_iterative'=>'no',
				'UN_query'=>'yes', 'query_len'=>16, 'alignment_mismatch'=>1, 'reporting_mode'=>'-a --best --strata',
				'suppress_num'=>3, 'statistical_items'=>'NA', 'unaligned_read_counts_filter'=>3, 
				'RC_background'=>1, 'threads_num'=>1,
			);
			
#get the directory of exe scripts involved in
our $perl_dir=Cwd::abs_path(Cwd::getcwd().'/..' );

#initiate directory and file
my $variables_pointer=sub_gui::initiate_eRNA_variables(\%variables, $perl_dir);
%variables=%$variables_pointer;

#Illumina TruSeq small RNA prep kits
#adapter 5: GTTCAGAGTTCTACAGTCCGACGATC
#adapter 3: TGGAATTCTCGGGTGCCAAGG

#setup shared hash to control thread
my %shash;
share(%shash); #will work for first level keys
$shash{'go'} = 0;
$shash{'predict'}=0;
$shash{'work'} = '';
$shash{'die'} = 0;


################################


##############################################################################
#standard window creation, placement, and signal connecting
my $window= sub_gui::default_window("eRNA: A tool for data analysis from high-throughput RNA sequencing", '900', '600');
 	my $table = Gtk2::Table->new(12,6,1); #the last '1' is true
	$table->set_border_width(10);
		#create a menu bar with two menu items, one will have $menu_edit as a submenu
		my $menu_bar = Gtk2::MenuBar->new;
			#add first menu item
			my $menu_item_outer= Gtk2::MenuItem->new('Outer GUIs');
			my $menu_outer=&menu_outer_GUIs();
			$menu_item_outer->set_submenu ($menu_outer);#set its sub menu
		$menu_bar->append($menu_item_outer);
			#add second menu item
			my $menu_item_view= Gtk2::MenuItem->new('Quality viewers');
			my $menu_view=&menu_view();
			$menu_item_view->set_submenu ($menu_view);#set its sub menu
		$menu_bar->append($menu_item_view);
			#add second menu item
			my $menu_item_tools = Gtk2::MenuItem->new('Mini tools');
			my $menu_tools=&menu_tools();
			$menu_item_tools->set_submenu ($menu_tools);#set its sub menu
		$menu_bar->append($menu_item_tools);
			#add third menu item
			my $menu_item_help = Gtk2::MenuItem->new('Help');
			$menu_item_help->set_sensitive(1); #1 is true
			$menu_item_help->set_right_justified(1);#1 is true push it off to the right
			my $menu_help=&menu_help();
			$menu_item_help->set_submenu ($menu_help);#set its sub menu
		$menu_bar->append($menu_item_help);	
	$table->attach_defaults($menu_bar,0,6,0,1);

	#add the second widget vpaned of vbox
	my $vpaned = Gtk2::VPaned->new;
	$vpaned->set_position (300);

	#add widget of notebook under vpaned
	my $notebook = Gtk2::Notebook->new;
		############ notebok_1
		my $nb_sw=&nb_self();
		my $nb_label = Gtk2::Label->new_with_mnemonic ('Self-testing');
		$notebook->append_page ($nb_sw, $nb_label);
		#### notebook: identification of known ncRNA
		$nb_sw=&nb_ncRNA();
		$nb_label = Gtk2::Label->new_with_mnemonic ('miRNA identification');
		$notebook->append_page ($nb_sw, $nb_label);
		#### notebook: identification of known mRNA
		$nb_sw=&nb_mRNA();
		$nb_label = Gtk2::Label->new_with_mnemonic ('mRNA identification');
		$notebook->append_page ($nb_sw, $nb_label);
		#notebook:expression profiling
		#$nb_sw=&nb_M4_targets();
		#$nb_label = Gtk2::Label->new_with_mnemonic ('Target screening');
		#$notebook->append_page ($nb_sw, $nb_label);
	$notebook->show_all();
	$vpaned->add1 ($notebook);

	#add widget of frame-scrollwindow-textview under vpaned
	my $frame=Gtk2::Frame->new();
	$frame->set_label('Output viewing');
	$frame->set_shadow_type('out');
	$frame->set_border_width(5);
		my $sw = Gtk2::ScrolledWindow->new ();
		$sw->set_shadow_type ('etched-out');
		$sw->set_policy ('automatic', 'automatic');
		$sw->set_border_width(10);
			my $text_view=Gtk2::TextView->new();
			$text_view->can_focus(1);
			$text_view->set_editable(0);
			$text_view->set_left_margin(5);
			$text_view->set_right_margin(5);
			$text_view->set_wrap_mode('GTK_WRAP_WORD_CHAR');
			#our $text_buffer = Gtk2::TextBuffer->new();
			#our $iter=$text_buffer->get_iter_at_offset(0);
			#$text_buffer->insert($iter, "Self-testing begins before running:\n\n\n");
			#$text_view->set_buffer($text_buffer);
			 #refresh screen per 1 second
		$text_view->show_all;
		$sw->add($text_view);
	$frame->add($sw);
	$vpaned->add2 ($frame);
	$vpaned->show_all();
	$table->attach_defaults ($vpaned, 0,6,1,11);

	#progress 
	$frame=Gtk2::Frame->new();
	$frame->set_label('Progress bar');
	$frame->set_shadow_type('out');
	$frame->set_border_width(5);
		my $pbar = Gtk2::ProgressBar->new();
		$pbar->set_pulse_step(.1);
		$pbar->hide; #needs to be called after show_all
	$frame->add($pbar);
	$table->attach_defaults ($frame, 0,3,11,12);
	
	$frame=Gtk2::Frame->new();
	$frame->set_shadow_type('out');
	$frame->set_border_width(5);
		my $button_1=Gtk2::Button->new('Refresh output');
		$button_1->signal_connect(clicked=>sub{
				my $content=`cat -e $variables{file_log}`; #refresh view
				my @lines=split("\n", $content);
				my $lines_num=(@lines>1000) ? 1000: @lines;
				my $last_1000=join("\n",@lines[-$lines_num..-1]);
				print "$variables{file_log}";
				my $text_buffer = $text_view->get_buffer();
				$text_buffer->set_text($last_1000);
		});
	$frame->add($button_1);
	$table->attach_defaults ($frame, 3,4,11,12);
	
	$frame=Gtk2::Frame->new();
	$frame->set_shadow_type('out');
	$frame->set_border_width(5);
		my $button_2=Gtk2::Button->new('Stop and quit');
		$button_2->signal_connect(clicked=>sub{
			$shash{'go'} =0;
			$shash{'die'} = 1;
			#foreach my $thread(threads->list()){
				#print "$thread\n";
				#$thread->join();
			#}
			$window->destroy;
			threads->exit;
			return 1; # 1 is true
		});
	$frame->add($button_2);
	$table->attach_defaults ($frame, 4,6,11,12);

$window->add($table);
$window->show_all;
Gtk2->main(); #our main event-loop



####################################################################
#used in widget of menu
sub menu_outer_GUIs{
	#Start of with a Gtk2::Menu widget
	my $menu_edit = Gtk2::Menu->new();
	
	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add  a tearoff item
		my $menu_item_cut = Gtk2::ImageMenuItem->new_from_stock ('Bowtie1 GUI', undef);#add an Image Menu item using stock
		$menu_item_cut->signal_connect('activate' => sub {
			my $perl_script=$perl_dir."/tool_bowtie1.pl";
			system("perl $perl_script $variables{file_var}");
		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item_cut);#add second: a item named cut

	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add  a tearoff item
		$menu_item_cut = Gtk2::ImageMenuItem->new_from_stock ('Bowtie2 GUI', undef);#add an Image Menu item using stock
		$menu_item_cut->signal_connect('activate' => sub {
			my $perl_script=$perl_dir."/tool_bowtie2.pl";
			system("perl $perl_script $variables{file_var}");
		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item_cut);#add second: a item named cut
	
	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add  a tearoff item
		$menu_item_cut = Gtk2::ImageMenuItem->new_from_stock ('miRDeep2 GUI', undef);#add an Image Menu item using stock
		$menu_item_cut->signal_connect('activate' => sub {
			my $perl_script=$perl_dir."/tool_mirdeep2.pl";
			system("perl $perl_script $variables{file_var}");
		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item_cut);#add second: a item named cut
	
		$menu_edit->append(Gtk2::TearoffMenuItem->new);#add  a tearoff item
		$menu_item_cut = Gtk2::ImageMenuItem->new_from_stock ('miRspring GUI', undef);#add an Image Menu item using stock
		$menu_item_cut->signal_connect('activate' => sub {
			my $perl_script=$perl_dir."/tool_mirspring.pl";
			system("perl $perl_script $variables{file_var}");
		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item_cut);#add second: a item named cut
	
	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add a tearoff item
	
	return($menu_edit);
}


###############################################################################
sub menu_view{
	#Start of with a Gtk2::Menu widget
	my $menu_edit = Gtk2::Menu->new();
	
	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add a tearoff item
		my $menu_item = Gtk2::ImageMenuItem->new_from_stock ('QS Viewer: quality scores', undef);#add an Image Menu item using stock
		$menu_item->signal_connect('activate' => sub { 
				my $perl_script=$perl_dir."/viewer_QS.pl";
				system("perl $perl_script $variables{file_var}");
		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item);#add first a item named cut
	
	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add a tearoff item
		$menu_item = Gtk2::ImageMenuItem->new_from_stock ('IL Viewer: insert length distribution', undef);#add an Image Menu item using stock
		$menu_item->signal_connect('activate' => sub { 
				my $perl_script=$perl_dir."/viewer_IL.pl";
				system("perl $perl_script $variables{file_var}");
		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item);#add first a item named cut
	
	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add a tearoff item
		$menu_item = Gtk2::ImageMenuItem->new_from_stock ('SD Viewer: sequencing depth', undef);#add an Image Menu item using stock
		$menu_item->signal_connect('activate' => sub { 
				my $perl_script=$perl_dir."/viewer_SD.pl";
				system("perl $perl_script $variables{file_var}");
		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item);#add first a item named cut
	
	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add a tearoff item
	return($menu_edit);
}
################################
#used in widget of menu
sub menu_tools{
	#Start of with a Gtk2::Menu widget
	my $menu_edit = Gtk2::Menu->new();
	
	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add a tearoff item
		my $menu_item_cut = Gtk2::ImageMenuItem->new_from_stock ('GBK analyzer', undef);#add an Image Menu item using stock
		$menu_item_cut->signal_connect('activate' => sub { 
			my $perl_script=$perl_dir."/tool_GBK_analyzer.pl";
			system("perl $perl_script");
		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item_cut);#add first a item named cut

	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add a tearoff item
		$menu_item_cut = Gtk2::ImageMenuItem->new_from_stock ('FASTQ pre-processor', undef);#add an Image Menu item using stock
		$menu_item_cut->signal_connect('activate' => sub { 
			my $perl_script=$perl_dir."/tool_FASTQ_pre-processor.pl";
			system("perl $perl_script");
		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item_cut);#add first a item named cut
	
	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add a tearoff item
	
	return($menu_edit);
}
#################################
#used in widget of menu
sub menu_help{
	#Start of with a Gtk2::Menu widget
	my $menu_edit = Gtk2::Menu->new();
	
	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add a tearoff item
		my $menu_item_cut = Gtk2::ImageMenuItem->new_from_stock ('Help document', undef);#add an Image Menu item using stock
		$menu_item_cut->signal_connect(activate => sub { 

		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item_cut);#add first a item named cut

	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add  a tearoff item
		$menu_item_cut = Gtk2::ImageMenuItem->new_from_stock ('About', undef);#add an Image Menu item using stock
		$menu_item_cut->signal_connect('activate' => sub {
			#standard window
			my $window = sub_gui::default_window('About E_RNA', '400', '300');
			$window->signal_connect('delete_event' => sub { Gtk2->main_quit; });
				$frame = Gtk2::Frame->new();
				$frame->set_border_width(5);
					my $table=Gtk2::Table->new(6,6, 1); # 1 is true
						my $label = Gtk2::Label->new("Version: 1.3");
					$table->attach_defaults($label, 0,6,0,1);
						$label = Gtk2::Label->new('Developer: Tiezheng Yuan, tyuan10@jhmi.edu');
					$table->attach_defaults($label, 0,6,1,2);
						$label = Gtk2::Label->new("Copyright: 2012.09-2016.02 TiezhengYuan\n   eRNA can be freely applied for non-commercial use!");
					$table->attach_defaults($label, 0,6,2,5);
						my $button=Gtk2::Button->new('Close');
					$table->attach_defaults($button, 2,4,5,6);
				$frame->add($table);
			$window->add($frame);
			$window->show_all();
			Gtk2->main();
		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item_cut);#add second: a item named cut

	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add a tearoff item
	
	return($menu_edit);
}


#########################################
sub nb_self_S1_basic_setup{
	use constant TRUE=>1;
	use constant FALSE=>1;

	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Basic setup', '900','500');
	my $table = Gtk2::Table->new(6,3,TRUE);
	print "##########3\n";
	#raw data
		my($frame_1, $entry_1, $button_1)=sub_gui::entry_control('Directory of raw data storage', $variables{dir_raw_data}, 'add');
	$table->attach_defaults($frame_1,0,3,0,1);
	#result
		my($frame_2, $entry_2, $button_2)=sub_gui::entry_control('Directory of result storage', $variables{dir_result_array}, 'add');
	$table->attach_defaults($frame_2,0,3,1,2);
	
	#frame: threads num
		my @entry_text=(1,2,4,8,16,24,32);
		my($frame_3, $cb_entry)=sub_gui::combobox_frame('Number of multiple threads', \@entry_text, 1);
	$table->attach_defaults($frame_3,0,1,2,3);
	
	#textview
		my($frame, $textview)=sub_gui::textview_control('Checking information');
	$table->attach_defaults($frame,0,3,3,6);
			
	$frame = Gtk2::Frame->new();
	$frame->set_border_width(5);
		my $hbox = Gtk2::HBox->new(FALSE,5);
			my $button_3=Gtk2::Button->new('Check');
			$button_3->set_border_width(10);
			$button_3->signal_connect(pressed=>sub{
				my $textbuffer=Gtk2::TextBuffer->new();
				my $iter=$textbuffer->get_iter_at_offset(0);
				my @rawdata_dirs=split(',', $entry_1->get_text());
				$textbuffer->insert($iter, "Storage of raw data\n");
				foreach my $dir(@rawdata_dirs){
					my $dir_judging= ( -d $dir) ? "\tYes\t$dir\n" : "\tNo\t$dir\n";
					$textbuffer->insert($iter, $dir_judging);
				}
				my @result_dirs=split(',', $entry_2->get_text());
				$textbuffer->insert($iter, "\nStorage of restuls\n");
				foreach my $dir(@result_dirs){
					my $dir_judging= ( -d $dir) ? "\tYes\t$dir\n" : "\tNo\t$dir\n";
					$textbuffer->insert($iter, $dir_judging);
				}
				$textview->set_buffer($textbuffer);
				$textview->show_all();
			});
		$hbox->pack_start($button_3,TRUE,TRUE,10);
			my $button_4=Gtk2::Button->new('Save and close');
			$button_4->set_border_width(10);
			$button_4->signal_connect('clicked' => sub { 
				$variables{dir_raw_data}=$entry_1->get_text();
				$variables{dir_result_array}=$entry_2->get_text();
				my @result_dirs=split(',', $entry_2->get_text());
				$variables{dir_result}=$result_dirs[0];
				#multi-threading
				$variables{threads_num}=$cb_entry->get_active_text();
				
				#refresh eRNA directory and file
				my $variables_pointer=sub_gui::initiate_eRNA_variables(\%variables);
				%variables=%$variables_pointer;
				#
				$window->destroy;
			} );
		$hbox->pack_start($button_4,TRUE,TRUE,10);
	$frame->add($hbox);
	$table->attach_defaults($frame,0,3,6,7);
	
	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}

###################
sub nb_self_S3_outer_tools{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Outer tools checking', '600','700');
	my $table = Gtk2::Table->new(8,3,TRUE);
	
	#Bowtie1
	my($frame_1, $entry_1, $button_1)=sub_gui::entry_control('Directory of the aligner: bowtie1', $variables{dir_bowtie1});
	my $row=0;
	$table->attach_defaults($frame_1,0,3,$row,$row+1);
	
	#Bowtie2
	my($frame_2, $entry_2, $button_2)=sub_gui::entry_control('Directory of the aligner: bowtie2', $variables{dir_bowtie2});
	$row++;
	$table->attach_defaults($frame_2,0,3,$row,$row+1);
	
	#tophat
	my($frame_3, $entry_3, $button_3)=sub_gui::entry_control('Directory of the mapper: tophat2', $variables{dir_mapper});
	$row++;
	$table->attach_defaults($frame_3,0,3,$row,$row+1);
	
	#cufflinks
	my($frame_4, $entry_4, $button_4)=sub_gui::entry_control('Directory of the assembler: cufflinks', $variables{dir_assembler});
	$row++;
	$table->attach_defaults($frame_4,0,3,$row,$row+1);
	
	#miRDeep
	my($frame_5, $entry_5, $button_5)=sub_gui::entry_control('Directory of the miRNA tool: miRDeep', $variables{dir_mirdeep});
	$row++;
	$table->attach_defaults($frame_5,0,3,$row,$row+1);
	
	#miRspring
	my($frame_6, $entry_6, $button_6)=sub_gui::entry_control('Directory of the miRNA tool: mirspring', $variables{dir_mirspring});
	$row++;
	$table->attach_defaults($frame_6,0,3,$row,$row+1);
	
	#checking information
	my($frame_textview, $textview)=sub_gui::textview_control('Checking information');
	$row++;
	$table->attach_defaults($frame_textview,0,2,$row,$row+3);
			
	my $frame = Gtk2::Frame->new();
	$frame->set_border_width(5);
	my $hbox = Gtk2::VBox->new(FALSE,5);
	my $button_7=Gtk2::Button->new('Check');
	$button_7->set_border_width(10);
	$button_7->signal_connect(pressed=>sub{
		my $textbuffer=Gtk2::TextBuffer->new();
		my $iter=$textbuffer->get_iter_at_offset(0);
		#bowtie checking
		my $dir=$entry_1->get_text();
		my $judging= ( -f $dir.'/bowtie') ? "Bowtie v.1 \t OK \n" : "Bowtie v.1 \t NO \n";
 		$textbuffer->insert($iter, $judging);
		#
		$dir=$entry_2->get_text(); 		
		$judging= ( -f $dir.'/bowtie2') ? "Bowtie v.2 \t OK \n" : "Bowtie v.2 \t NO \n";
		$textbuffer->insert($iter, $judging);
		#TopHat checking
		$dir=$entry_3->get_text();
		$judging= ( -f $dir.'/tophat2') ?  "TopHat v.2 \t OK \n" : "TopHat v.2 \t NO \n";
		$textbuffer->insert($iter, $judging);
		#cufflinks checking
		$dir=$entry_4->get_text();
		$judging= ( -f $dir.'/cufflinks') ? "cufflinks \t OK \n" : "cufflinks \t NO \n";
		$textbuffer->insert($iter, $judging);
		$judging= ( -f $dir.'/cuffmerge') ? "cuffmerge \t OK \n" : "cuffmerge \t NO \n";
		$textbuffer->insert($iter, $judging);
		$judging= ( -f $dir.'/cuffdiff') ? "cuffdiff \t OK \n" : "cuffdiff \t NO \n";
		$textbuffer->insert($iter, $judging);
		#miRDeep2 checking
		$dir=$entry_5->get_text();
		$judging= ( -f $dir.'/miRDeep2.pl') ? "miRDeep2.pl \t OK \n" : "miRDeep2.pl \t NO \n";
		$textbuffer->insert($iter, $judging);
		#miRspring checking
		$dir=$entry_6->get_text();
		$judging= ( -f $dir.'/BAM_to_Intermediate.pl' and -f $dir.'/Intermediate_to_miRspring.pl' and -f $dir.'/javascriptTEMPLATE.txt') ? "miRspring \t OK \n" : "miRspring \t NO \n";
		$textbuffer->insert($iter, $judging);
		
		$textview->set_buffer($textbuffer);
		$textview->show_all();
	},);
	$hbox->pack_start($button_7,TRUE,TRUE,10);
	my $button_8=Gtk2::Button->new('Save and Close');
	$button_8->set_border_width(10);
	$button_8->signal_connect('clicked' => sub {
		$variables{dir_bowtie1}=$entry_1->get_text();
		$variables{dir_bowtie2}=$entry_2->get_text();
		$variables{dir_assembler}=$entry_4->get_text();
		$variables{dir_mapper}=$entry_3->get_text();
		$variables{dir_mirdeep}=$entry_5->get_text();
		$variables{dir_mirspring}=$entry_6->get_text();
		#initiate variables in variables.txt
		sub_common::hash_to_file(\%variables, $variables{file_var}, '=');
		#
		$window->destroy;
	} );
	$hbox->pack_start($button_8,TRUE,TRUE,10);
	$frame->add($hbox);
	$table->attach_defaults($frame,2,3,$row,$row+3);
	
	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}

##################################
sub nb_self_S2_sample_management{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Sample management', 700, 400);
	my $table = Gtk2::Table->new(7,4,TRUE);
	
	#get sample names and fastq files
	my (%sample_hash, %fq_sample);
	my $sample_names_pointer=sub_basic::auto_sample_names($variables{dir_raw_data}, $variables{raw_file_format} );
	my @sample_names=@$sample_names_pointer;
	my $files_pointer=sub_common::files_list($variables{dir_raw_data}, 'incrusive_file');
	my @files=@$files_pointer;
	my @raw_files=grep(/\.fastq$|\.fq$/, @files);
	
		#the list of sample names and fastq names
		my $vbox = Gtk2::VBox->new;
			my $scroll = Gtk2::ScrolledWindow->new;
			$scroll->set_policy ('never', 'automatic');

			# create and load the model
			my $liststore = Gtk2::ListStore->new ('Glib::String', 'Glib::String');
			foreach (@sample_names) {
				my $iter = $liststore->append;
				$liststore->set ($iter, 0, $_, 1, $_);
			}

			# now a view
			my $treeview_1 = Gtk2::TreeView->new ($liststore);
			$treeview_1->set_rules_hint (TRUE);
			$treeview_1->set_reorderable (TRUE);
				# regular editable text column for column 0
				my $renderer = Gtk2::CellRendererText->new;
				$renderer->set (editable => FALSE);
				my $column = Gtk2::TreeViewColumn->new_with_attributes ('Fastq names', $renderer, text => 0);
				$renderer->signal_connect (edited => sub {
							my ($cell, $text_path, $new_text, $model) = @_;
							my $path = Gtk2::TreePath->new_from_string ($text_path);
							my $iter = $model->get_iter ($path);
							$model->set ($iter, 0, $new_text);
						}, $liststore);
			$treeview_1->append_column ($column);
				#renderer column 1
				$renderer = Gtk2::CellRendererText->new;
				$renderer->set (editable => TRUE);
				$column = Gtk2::TreeViewColumn->new_with_attributes ('Sample names', $renderer, text => 1);
				$renderer->signal_connect (edited => sub {
							my ($cell, $text_path, $new_text, $model) = @_;
							my $path = Gtk2::TreePath->new_from_string ($text_path);
							my $iter = $model->get_iter ($path);
							$model->set ($iter, 1, $new_text);
						}, $liststore);
			$treeview_1->append_column ($column);

			$scroll->add ($treeview_1);
		$vbox->pack_start ($scroll, TRUE, TRUE, 0);
	$table->attach_defaults($vbox,0,3,0,3);
	
		#frame: check sample_info
		$frame= Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Checking information');
			my $sw = Gtk2::ScrolledWindow->new ();
			$sw->set_shadow_type ('etched-out');
			$sw->set_policy ('automatic', 'automatic');
			$sw->set_border_width(10);
				my $textview_2 = Gtk2::TextView->new();
				$textview_2->set_wrap_mode('GTK_WRAP_WORD_CHAR');
			$sw->add($textview_2);
		$frame->add($sw);
	$table->attach_defaults($frame,0,4,3,7);
	
		#button
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button_1=Gtk2::Button->new('Check');
			$button_1->set_border_width(5);
			$button_1->signal_connect(clicked=>sub{
						%sample_hash=(); #empty out hash
						#$liststore=$textview_2->get_selection;
						$liststore->foreach(sub{
							my($model, $path, $iter)=@_;
							my $fastq=$model->get ($iter, 0);
							my $sample=$model->get ($iter, 1);
							$sample_hash{$sample}= $sample_hash{$sample} ? $sample_hash{$sample}.','.$fastq : $fastq;
							$fq_sample{$fastq}=$sample;
							#print %sample_hash;
							FALSE;
						});
						print %sample_hash;
						#connection between sample name and fastq files
						my $textbuffer=Gtk2::TextBuffer->new();
						my $iter=$textbuffer->get_iter_at_offset(0);
						my $n=1;
						foreach my $sample_name(keys %sample_hash){#4
							my @raw_names=split(",", $sample_hash{$sample_name});
							my @sub_raw_files;
							foreach my $raw_name(@raw_names){
								foreach my $raw_file(@raw_files){
									push(@sub_raw_files, $raw_file) if $raw_file=~/\/$raw_name/;
								}
							}
							my $files_num=@sub_raw_files;
							$textbuffer->insert($iter, $n.": $sample_name (number of fastq files=$files_num):\n");
							foreach(@sub_raw_files){
								$textbuffer->insert($iter, "\t$_\n");
							}
							$n++;
						}#4
					$textview_2->set_buffer($textbuffer);
					$textview_2->show_all();
				});
		$frame->add($button_1);
	$table->attach_defaults($frame,3,4,0,1);
	
		#button
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button_2=Gtk2::Button->new('Save');
			$button_2->set_border_width(5);
			$button_2->signal_connect(pressed=>sub{
							my $sample_info_file=$variables{dir_result}.'/sample_info.csv';
							sub_common::hash_to_file(\%fq_sample, $sample_info_file, ",");
						});
		$frame->add($button_2);
	$table->attach_defaults($frame,3,4,1,2);
		
		# Close button
		my $frame_close = sub_gui::close_button_control($window);
	$table->attach_defaults($frame_close,3,4,2,3);
	
	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}


#########################
sub nb_self_S4_package{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Package dependency checking', '400','400');
	my $table = Gtk2::Table->new(6,3,TRUE);
	
		my $frame_1 = Gtk2::Frame->new();
		$frame_1->set_border_width(5);
		$frame_1->set_label('Checking information');
			my $textview = Gtk2::TextView->new();
			$textview->set_left_margin(5);
			$textview->set_wrap_mode('GTK_WRAP_WORD_CHAR');
		$frame_1->add($textview);
	$table->attach_defaults($frame_1,0,3,0,5);
			
		my $frame_2 = Gtk2::Frame->new();
		$frame_2->set_border_width(5);
			my $hbox = Gtk2::HBox->new(FALSE,5);
				my $button_1=Gtk2::Button->new('Check');
				$button_1->set_border_width(10);
				$button_1->signal_connect(pressed=>sub{
					my $textbuffer=Gtk2::TextBuffer->new();
					my $iter=$textbuffer->get_iter_at_offset(0);
					$textbuffer->insert($iter, "Perl package dependency checking:\n");
					my @ISA=qw(Cwd File::Find Glib Gtk2 List::Util List::MoreUtils threads threads::shared 
									Bio::Perl Bio::SeqIO Bio::Seq::Quality Statistics::R );
					foreach (sort @ISA){
						$textbuffer->insert($iter, (eval "require $_") ? "\tOK\t"."$_\n" : "\tNO\t"."$_\n");
					}
					$textview->set_buffer($textbuffer);
					$textview->show_all();
				});
			$hbox->pack_start($button_1,TRUE,TRUE,10);
				my $button_2=Gtk2::Button->new('Close');
				$button_2->set_border_width(10);
				$button_2->signal_connect(clicked => sub { $window->destroy;}  );
			$hbox->pack_start($button_2,TRUE,TRUE,10);
		$frame_2->add($hbox);
	$table->attach_defaults($frame_2,0,3,5,6);
	
	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}
#############################################
sub nb_self{
	my $scroller = Gtk2::ScrolledWindow->new();
	$scroller->set_shadow_type ('etched-out');
	$scroller->set_policy ('automatic', 'automatic');
	$scroller->set_size_request ('200', '200');

	my $table_stock=Gtk2::Table->new(3, 5, TRUE);
	$table_stock->set_border_width(30);
	$table_stock->set_row_spacings(30);
	$table_stock->set_size_request (100, 100);
	$table_stock->attach_defaults(&sub_gui::nb_button('Basic setup', \&nb_self_S1_basic_setup), 0,2,0,1);
	$table_stock->attach_defaults(Gtk2::Arrow->new('right', 'etched-out'), 2,3,0,1);
	$table_stock->attach_defaults(&sub_gui::nb_button('Sample management', \&nb_self_S2_sample_management), 3,5,0,1);
	$table_stock->attach_defaults(Gtk2::Arrow->new('down', 'etched-out'), 4,5,1,2);
	$table_stock->attach_defaults(&sub_gui::nb_button('Outer tools checking', \&nb_self_S3_outer_tools), 3,5,2,3);
	$table_stock->attach_defaults(Gtk2::Arrow->new('left', 'etched-out'), 2,3,2,3);
	$table_stock->attach_defaults(&sub_gui::nb_button('Package dependency checking', \&nb_self_S4_package), 0,2,2,3);
	$scroller->add_with_viewport($table_stock);

	return($scroller);
}
######################################################################
sub nb_ncRNA{
	my $scroller = Gtk2::ScrolledWindow->new();
	$scroller->set_shadow_type ('etched-out');
	$scroller->set_policy ('automatic', 'automatic');
	$scroller->set_size_request (200, 200);
	
	my $table_stock=Gtk2::Table->new(3, 8, TRUE);
	$table_stock->set_border_width(30);
	$table_stock->set_row_spacings(30);
	$table_stock->set_size_request (100, 100);
	$table_stock->attach_defaults(&sub_gui::nb_button('Reference miRNAs', \&nb_ncRNA_S1_reference_miRNAs), 0,2,0,1);
	$table_stock->attach_defaults(Gtk2::Arrow->new('right', 'etched-out'), 2,3,0,1);
	
	$table_stock->attach_defaults(&sub_gui::nb_button('Sequence trimming', \&nb_ncRNA_S2_sequence_trimming), 3,5,0,1);
	$table_stock->attach_defaults(Gtk2::Arrow->new('right', 'etched-out'), 5,6,0,1);
	
	$table_stock->attach_defaults(&sub_gui::nb_button('Sequence alignment', \&nb_ncRNA_S3_sequence_alignment), 6,8,0,1);
	$table_stock->attach_defaults(Gtk2::Arrow->new('down', 'etched-out'), 6,8,1,2);
	
	$table_stock->attach_defaults(&sub_gui::nb_button('Reads counting', \&nb_ncRNA_S4_reads_counting), 6,8,2,3);
	$table_stock->attach_defaults(Gtk2::Arrow->new('left', 'etched-out'), 5,6,2,3);
	
	my $tbutton=Gtk2::Button->new('Apply and Run');
	$tbutton->signal_connect(clicked=>sub{
		my $thread_monitor=threads->new( sub{
			while(1){
				return if $shash{'die'} == 1;
				if ( $shash{'go'} == 1 ){
					my $perl_script=$perl_dir.'/time_monitor.pl';
					system("perl $perl_script $variables{file_var}");
				}
				else{	sleep 1;	}
			}
		});

	my $thread_main=threads->new( sub{
		#local $SIG{KILL} = sub { threads->exit };
		while(1){
			return if $shash{'die'} == 1;

			if ( $shash{'go'} == 1 ){
				my $perl_script=$perl_dir."/eRNA_ncRNA.pl";
				system("perl $perl_script $variables{file_var} > $variables{file_log}");
				$shash{'go'} = 0; #turn off main running 
			}
			else{	sleep 1;	}
		}
	});

		if($shash{'go'} == 0){
			$shash{'go'} = 1;
			$pbar->show;
			Glib::Timeout->add (100, sub {
				if($shash{'go'} == 1){
					$pbar->set_text('Running!');
					$pbar->pulse;
					return TRUE;
				}
				else{	
					$pbar->set_text('OK! It is done!');
					return FALSE;
				}
			});
		}
	});  
	$table_stock->attach_defaults($tbutton, 2,5,2,3);
	$scroller->add_with_viewport($table_stock);
	return($scroller);
}


###########################
sub nb_ncRNA_S1_reference_miRNAs{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Options for extraction of reference miRNAs', '450','300');
	my $table = Gtk2::Table->new(4,6,TRUE);
	print "##########3\n";
		#miRNA species
		my @entry_text=sort(keys %species);
		my($frame, $cb_species)=sub_gui::combobox_frame('Species of microRNA', \@entry_text, 28);
	$table->attach_defaults($frame, 0,6,0,1);

		my $button=Gtk2::Button->new('Fetch sequences from miRbase');
		$button->set_border_width(10);
		$button->signal_connect(pressed=>sub{
			$variables{ncRNA_species}=$species{ $cb_species->get_active_text() };
			sub_basic::refresh_log($variables{file_var}, 'ncRNA_species', $variables{ncRNA_species} );
			
			#get fasta files of precursor and mature miRNA from miRbase
			sub_ncrna::fetch_mirbase(\%variables);
			$variables{ncRNA_mirbase_precursor_file}=$variables{dir_bowtie1}.'hairpin.fa';
			$variables{ncRNA_mirbase_mature_file}=$variables{dir_bowtie1}.'mature.fa';
			
			#subroutine: get miRNA of a certain specie for miRBase fasta files
			sub_ncrna::specie_miRNA(\%variables);
			sub_gui::popup_window("Note", "get precursor and mature miRNA of ".$variables{ncRNA_species}." from miRBase fasta files!");
		} );
	$table->attach_defaults($button,0,6,1,2);
		$button=Gtk2::Button->new('Close');
		$button->set_border_width(10);
		$button->signal_connect(pressed=>sub{		$window->destroy;		} );
	$table->attach_defaults($button,4,6,3,4);
	
	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}

#####################################
sub nb_ncRNA_S2_sequence_trimming{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Options for low quality and adapter sequences removal', '800','350');
	my $table = Gtk2::Table->new(4,5,TRUE);
	
		#frame
		my @entry_text=(10,13,20,30,40);
		my($frame_1, $cb_entry_Q_filter)=sub_gui::combobox_frame('Quality control filter (Q value)', \@entry_text, 1);
	$table->attach_defaults($frame_1,0,3,0,1);

		#frame
		@entry_text=map {10**$_} (0..10);
		my($frame_2, $cb_entry_Q_compression)=sub_gui::combobox_frame('Quality control filter (Q value)', \@entry_text, 3);
	$table->attach_defaults($frame_2,3,6,0,1);
	
		#frame: 3'adapter
		my($frame_3, $entry_3end_adapter)=sub_gui::entry_frame("Sequences of 3' adapter (5'-)", $variables{adapter_3});
	$table->attach_defaults($frame_3, 0,3,1,2);
	
		#frame: 5'adapter
		my($frame_4, $entry_5end_adapter)=sub_gui::entry_frame("Sequences of 5' adapter (5'-)", $variables{adapter_5});
	$table->attach_defaults($frame_4, 3,6,1,2);

		#frame: match len
		@entry_text=(6..16);
		my($frame_5, $cb_match_len)=sub_gui::combobox_frame('Length of adapter sequence with exact match', \@entry_text, 2);
	$table->attach_defaults($frame_5, 0,3,2,3);
	
		#frame:mismatching of adapte sequences
		my $check_button=Gtk2::CheckButton->new('One mismatch is allowed in matching seqences.');
		$check_button->set_active(TRUE);
	$table->attach_defaults($check_button,3,6,2,3);
	
		#frame:query len
		@entry_text=(10..40);
		my($frame_6, $cb_query_len)=sub_gui::combobox_frame('Minimum query length', \@entry_text, 8);
	$table->attach_defaults($frame_6, 0,3,3,4);
	
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button=Gtk2::Button->new('Save and close');
			$button->signal_connect(clicked=>sub{
				#read raw data
				sub_basic::refresh_log($variables{file_var}, 'quality_filter', $cb_entry_Q_filter->get_active_text() );
				sub_basic::refresh_log($variables{file_var}, 'QC_compression', $cb_entry_Q_compression->get_active_text() );
				#print "$variables{quality_control}\n$variables{read_rawdata}\n";
				#adapter removal
				$variables{adapter_3}=sub_basic::format_DNA( $entry_3end_adapter->get_text() );
				$variables{adapter_5}=sub_basic::format_DNA( $entry_5end_adapter->get_text() );
				sub_basic::refresh_log($variables{file_var}, 'adapter_3', $variables{adapter_3} );
				sub_basic::refresh_log($variables{file_var}, 'adapter_5', $variables{adapter_5} );
				sub_basic::refresh_log($variables{file_var}, 'match_len', $cb_match_len->get_active_text );
				sub_basic::refresh_log($variables{file_var}, 'query_len', $cb_query_len->get_active_text() );
				sub_basic::refresh_log($variables{file_var}, 'ncRNA_adapter_removal', 'yes');
				
				if ($check_button->get_active ){
					sub_basic::refresh_log($variables{file_var}, 'mismatch_allowed', 'yes' );
					sub_basic::refresh_log($variables{file_var}, 'mismatch_len', $variables{match_len}+4);
				}
				else{	sub_basic::refresh_log($variables{file_var}, 'mismatch_allowed', 'no');		}
				
				sub_basic::refresh_log($variables{file_var}, 'ncRNA_read_rawdata', 'yes');
				$window->destroy;
				});
		$frame->add($button);
	$table->attach_defaults($frame,1,5,4,5);
		
	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}

#######################
sub nb_ncRNA_S3_sequence_alignment{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Options of sequence alignment', '800','700');
	my $table = Gtk2::Table->new(11,6,TRUE);
	my $file_names_pointer=sub_common::files_list($variables{dir_bowtie1}, 'file_name');
	my @file_names=@$file_names_pointer;
	my @index_names=grep(s/\.fa$|\.fasta$//, @file_names); #get index information
	print $variables{dir_bowtie1};
	print @index_names;
	
	##################bowtie1
	#Frame:-m, default -m 1 as unique alignment
	my @entry_text=(1..5);
	my($frame_1, $cb_supressed_multi)=sub_gui::combobox_frame("Cutoff of surpressed multiple alignment (-m)", \@entry_text, 1);
	$table->attach_defaults($frame_1, 0,3,0,1);
	
	#Frame: -l 
	@entry_text=(20..40);
	my($frame_2, $cb_seed_length)=sub_gui::combobox_frame("Seed length (-l)", \@entry_text, 4);
	$table->attach_defaults($frame_2, 3,6,0,1);
	
	#Frame:-v/-n
	@entry_text=('-v 0', '-v 1', '-v 2', '-v 3', '-n 0', '-n 1', '-n 2');
	my($frame_3, $cb_mismatches)=sub_gui::combobox_frame("Maximum mismatches (-n/-v)", \@entry_text, 1);
	$table->attach_defaults($frame_3, 0,3,1,2);
	
	#Frame: --nofw/--norc
	@entry_text=('Forward/Reverse', 'Forward only','Reverse only');
	my($frame_4, $cb_strand_alignment)=sub_gui::combobox_frame("Strand specific alignmen", \@entry_text, 1);
	$table->attach_defaults($frame_4, 3,6,1,2);

	#Frame:  items selection
	my $frame = Gtk2::Frame->new();
	$frame->set_border_width(5);
		my $check_button_1 = Gtk2::CheckButton->new ("Separate alignment: references selection");
		$check_button_1->set_active (TRUE);
			my $sub_table_3=items_selection(\@index_names, $variables{file_var}, 'index_seperate');
			$check_button_1->signal_connect (clicked => sub {
				my $self = shift;
				if ($self->get_active){
					$sub_table_3->set_sensitive (TRUE);
					$sub_table_3->show_all();
				}
				else{	$sub_table_3->set_sensitive (FALSE);	}
			});
	$frame->add($sub_table_3);
	$frame->set_label_widget($check_button_1);
	$table->attach_defaults($frame,0,6,2,6);
	
	#Frame 4 items selection
	$frame = Gtk2::Frame->new();
	$frame->set_border_width(5);
		my $check_button_2 = Gtk2::CheckButton->new ("Iterative alignment: references selection");
		$check_button_2->set_active (TRUE);
		my $sub_table_4=items_selection(\@index_names, $variables{file_var}, 'index_iterative' );
		$check_button_2->signal_connect (clicked => sub {
			my $self = shift;
			if ($self->get_active){
				$sub_table_4->set_sensitive (TRUE);
				$sub_table_4->show_all();
			}
			else{	$sub_table_4->set_sensitive (FALSE);	}
		});
	$frame->add($sub_table_4);
	$frame->set_label_widget($check_button_2);
	$table->attach_defaults($frame,0,6,6,10);
	
	#Frame of save and close
	$frame = Gtk2::Frame->new();
	$frame->set_border_width(5);
	my $button=Gtk2::Button->new('Save and close');
	$button->signal_connect(clicked=>sub{
		#print "$variables{seperate_alignment} \t$variables{iterative_alignment} \n";
		#bowtie options
		$variables{bowtie_options}=join(' ',$variables{dir_bowtie1}.'bowtie', ' --best --strata',
									'-l', $cb_seed_length->get_active_text(), 
									'-m', $cb_supressed_multi->get_active_text(),
									$cb_mismatches->get_active_text(),	);
		if ($cb_strand_alignment->get_active_text() eq 'Forward only'){
			$variables{bowtie_options} .= ' --norc';
		}
		elsif ($cb_strand_alignment->get_active_text() eq 'Reverse only'){
			$variables{bowtie_options} .= ' --nofw';
		}
		#print "$variables{bowtie_options}\n";
		$variables{ncRNA_seperate_alignment}= 'yes' unless sub_basic::read_log($variables{file_var}, 'index_seperate') eq 'no';
		$variables{ncRNA_iterative_alignment}= 'yes' unless sub_basic::read_log($variables{file_var}, 'index_iterative') eq 'no';
		#
		sub_basic::refresh_log($variables{file_var}, 'bowtie_options', $variables{bowtie_options} );
		sub_basic::refresh_log($variables{file_var}, 'ncRNA_seperate_alignment', $variables{ncRNA_seperate_alignment} );
		sub_basic::refresh_log($variables{file_var}, 'ncRNA_iterative_alignment', $variables{ncRNA_iterative_alignment} );
		$window->destroy;
	});
	$frame->add($button);
	$table->attach_defaults($frame,1,5,10,11);

	######################################
	$table->show_all();
	$window->add($table);
	$window->show_all();
	Gtk2->main();
}
###############################################################
#used for items' selection
sub items_selection{
	my ($index_names_pointer, $variable_file, $variable_name)=@_;
	my @index_names=@$index_names_pointer;
	my @selected_names;
	
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
		sub_basic::refresh_log($variable_file, $variable_name, join(",", @selected_names) );
	} );
	$table->attach_defaults($button_1, 3,4,1,2);
	my $button_2=Gtk2::Button->new('>>');        ###### >> button: all transfered
	$button_2->signal_connect(clicked=>sub{
		$tree_store_1->clear;
		foreach (@index_names) {
			my $iter = $tree_store_2->append(undef); #the iter is a pointer in the treestore. We use to add data.
			$tree_store_2->set ($iter,0 => $_);
		}
		sub_basic::refresh_log($variable_file, $variable_name, join(",", @index_names) );
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
		sub_basic::refresh_log($variable_file, $variable_name, join(",", @selected_names) );
	} );
	$table->attach_defaults($button_3, 3,4,3,4);
	my $button_4=Gtk2::Button->new('<<');      ######<< button
	$button_4->signal_connect(clicked=>sub{
		$tree_store_2->clear;
		foreach (@index_names) {
			my $iter = $tree_store_1->append(undef); #the iter is a pointer in the treestore. We use to add data.
			$tree_store_1->set ($iter,0 => $_);
		}
		sub_basic::refresh_log($variable_file, $variable_name, 'NA' );
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
		for(my $i=0; $i<@selected_names; $i++){
			if ( ($selected_names[$i] eq $value) and $i>0){
				$selected_names[$i-1]=$value;
				$selected_names[$i]=$pre_value;
			}
		}
		@selected_names=grep($_, @selected_names);
		sub_basic::refresh_log($variable_file, $variable_name, join(",", @selected_names) );
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
		for(my $i=0; $i<@selected_names; $i++){
			if ( ($selected_names[$i] eq $value) and $i<@selected_names-1){
				$selected_names[$i+1]=$value;
				$selected_names[$i]=$next_value;
			}
		}
		@selected_names=grep($_, @selected_names);
		sub_basic::refresh_log($variable_file, $variable_name, join(",", @selected_names) );
	} );
	$table->attach_defaults($button_6, 7,8,4,5);

	return ($table);
}

############################
sub nb_ncRNA_S4_reads_counting{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Options of reads counting', '300','250');
	my $table = Gtk2::Table->new(3,4,TRUE);
	
		#frame:read counts noise
		my @entry_text=(1..20);
		my($frame_2, $cb_entry_2)=sub_gui::combobox_frame('Reads backgroud ', \@entry_text, 0);
	$table->attach_defaults($frame_2,0,4,0,1);
	
		#frame:normalization
		@entry_text=(1e3, 1e6,1e9);
		my($frame_3, $cb_entry_3)=sub_gui::combobox_frame('Scaling of normalization', \@entry_text, 1);
	$table->attach_defaults($frame_3,0,4,1,2);
	
		#Frame 5
		my $frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button=Gtk2::Button->new('Save and close');
			$button->set_border_width(5);
			$button->signal_connect(clicked=>sub{
				$variables{RC_background}=$cb_entry_2->get_active_text();
				sub_basic::refresh_log($variables{file_var}, 'RC_background', $variables{RC_background});
				$variables{normalization_scaling}=$cb_entry_3->get_active_text();
				sub_basic::refresh_log($variables{file_var}, 'normalization_scaling', $variables{normalization_scaling} );
				sub_basic::refresh_log($variables{file_var}, 'ncRNA_counting', 'yes' );
				$window->destroy;
			});
		$frame->add($button);
	$table->attach_defaults($frame,0,4,2,3);
		
	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}

############################
#used for items selection
sub multiple_items_selection{
	my ($index_names_pointer, $variable_file, $variable_name_A, $variable_name_B)=@_;
	my @index_names=@$index_names_pointer;
	my (@selected_A, @selected_B);
	
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
	my $button_3=Gtk2::Button->new('<');      ######<< button
	$button_3->signal_connect(clicked=>sub{
		my $tree_selection=$tree_view_1->get_selection;
		my($tree_store, $iter_1)=$tree_selection->get_selected;
		my $value=$tree_store->get($iter_1, 0);
		return unless $iter_1;
		$tree_store_1->remove($iter_1);
		my $iter=$tree_store_2->append(undef);  #add the value into the second tree 
		$tree_store_2->set($iter,0=> $value);
		sub_basic::refresh_log($variable_file, $variable_name_A, sub_gui::tree_store_str($tree_store_2) );
	} );
	$table->attach_defaults($button_3, 3,4,1,2);
	my $button_1=Gtk2::Button->new('>');        ######>> button: one transfered
	$button_1->signal_connect(clicked=>sub{
		my $tree_selection=$tree_view_2->get_selection;
		my($tree_store, $iter_1)=$tree_selection->get_selected;
		my $value=$tree_store->get($iter_1, 0);
		return unless $iter_1;
		$tree_store->remove($iter_1);
		my $iter=$tree_store_1->append(undef);  #add the value into the second tree 
		$tree_store_1->set($iter,0=> $value);
		sub_basic::refresh_log($variable_file, $variable_name_A, sub_gui::tree_store_str($tree_store_2) );
	} );
	$table->attach_defaults($button_1, 3,4,3,4);
	my $button_2=Gtk2::Button->new('>>');        ###### >> button: all transfered
	$button_2->signal_connect(clicked=>sub{
		my @names=split(',', sub_gui::tree_store_str($tree_store_2) );
		foreach(@names){
			unless($_ eq 'NA'){
				my $iter=$tree_store_1->append(undef);  #add the value into the second tree 
				$tree_store_1->set($iter,0=> $_);
			}
		}
		$tree_store_2->clear;
		sub_basic::refresh_log($variable_file, $variable_name_A, 'NA' );
	} );
	$table->attach_defaults($button_2, 3,4,4,5);

	
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
		my $button_6=Gtk2::Button->new('>');        ######>> button: one transfered
	$button_6->signal_connect(clicked=>sub{
		my $tree_selection=$tree_view_1->get_selection;
		my($model, $iter_1)=$tree_selection->get_selected;
		my $value=$model->get($iter_1, 0);
		return unless $iter_1;
		$tree_store_1->remove($iter_1);
		my $iter=$tree_store_3->append(undef);  #add the value into the second tree 
		$tree_store_3->set($iter,0=> $value);
		sub_basic::refresh_log($variable_file, $variable_name_B, sub_gui::tree_store_str($tree_store_3) );
	} );
	$table->attach_defaults($button_6, 7,8,1,2);
	my $button_4=Gtk2::Button->new('<');      ######<< button
	$button_4->signal_connect(clicked=>sub{
		my $tree_selection=$tree_view_3->get_selection;
		my($tree_store, $iter_1)=$tree_selection->get_selected;
		my $value=$tree_store->get($iter_1, 0);
		return unless $iter_1;
		$tree_store_3->remove($iter_1);
		my $iter=$tree_store_1->append(undef);  #add the value into the second tree 
		$tree_store_1->set($iter,0=> $value);
		sub_basic::refresh_log($variable_file, $variable_name_B, tree_store($tree_store_3) );
	} );
	$table->attach_defaults($button_4, 7,8,3,4);

	my $button_5=Gtk2::Button->new('<<');        ###### >> button: all transfered
	$button_5->signal_connect(clicked=>sub{
		my @names=split(',', sub_gui::tree_store_str($tree_store_3));
		foreach (@names) {
			unless($_ eq 'NA'){
				my $iter = $tree_store_1->append(undef); #the iter is a pointer in the treestore. We use to add data.
				$tree_store_1->set ($iter,0 => $_);
			}
		}
		$tree_store_3->clear;
		sub_basic::refresh_log($variable_file, $variable_name_B, 'NA' );
	} );
	$table->attach_defaults($button_5, 7,8,4,5);


	$table->show_all();
	return ($table);
}

######################################################################
sub nb_mRNA{
	my $scroller = Gtk2::ScrolledWindow->new();
	$scroller->set_shadow_type ('etched-out');
	$scroller->set_policy ('automatic', 'automatic');
	$scroller->set_size_request (200, 200);

	my $table_stock=Gtk2::Table->new(3, 8, TRUE);
	$table_stock->set_border_width(30);
	$table_stock->set_row_spacings(30);
	$table_stock->set_size_request (100, 100);
	$table_stock->attach_defaults(&sub_gui::nb_button('Reference selection', \&nb_mRNA_S1_reference_sequences), 0,2,0,1);
	$table_stock->attach_defaults(Gtk2::Arrow->new('right', 'etched-out'), 2,3,0,1);
	$table_stock->attach_defaults(&sub_gui::nb_button('Spliced alignment', \&nb_mRNA_S2_spliced_alignment), 3,5,0,1);
	$table_stock->attach_defaults(Gtk2::Arrow->new('right', 'etched-out'), 5,6,0,1);
	$table_stock->attach_defaults(&sub_gui::nb_button('Transcript assembling', \&nb_mRNA_S3_transcript_assembling), 6,8,0,1);
	$table_stock->attach_defaults(Gtk2::Arrow->new('down', 'etched-out'), 6,8,1,2);
	$table_stock->attach_defaults(&sub_gui::nb_button('Differential expression', \&nb_mRNA_S4_differential_transcripts), 6,8,2,3);
	$table_stock->attach_defaults(Gtk2::Arrow->new('left', 'etched-out'), 5,6,2,3);
	my $tbutton=Gtk2::Button->new('Apply and Run');
	$tbutton->signal_connect(clicked=>sub{
		my $thread_monitor=threads->new( sub{
			while(1){
				return if $shash{'die'} == 1;
				if ( $shash{'go'} == 1 ){
					my $perl_script=$perl_dir.'/time_monitor.pl';
					system("perl $perl_script $variables{file_var}");
				}
				else{	sleep 1;	}
			}
		});

	my $thread_main=threads->new( sub{
		#local $SIG{KILL} = sub { threads->exit };
		while(1){
			return if $shash{'die'} == 1;

			if ( $shash{'go'} == 1 ){
				my $perl_script=$perl_dir."/eRNA_mRNA.pl";
				system("perl $perl_script $variables{file_var} > $variables{file_log}");
				$shash{'go'} = 0; #turn off main running 
			}
			else{	sleep 1;	}
		}
	});

		if($shash{'go'} == 0){
			$shash{'go'} = 1;
			$pbar->show;
			Glib::Timeout->add (100, sub {
				if($shash{'go'} == 1){
					$pbar->set_text('Running!');
					$pbar->pulse;
					return TRUE;
				}
				else{	
					$pbar->set_text('OK! It is done!');
					return FALSE;
				}
			});
		}
	});  
	$table_stock->attach_defaults($tbutton, 2,5,2,3);
	$scroller->add_with_viewport($table_stock);
	return($scroller);
}

########################
sub nb_mRNA_S1_reference_sequences{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Reference selection', '400','250');
	my $table = Gtk2::Table->new(3,4,TRUE);
	
	#get index information
	my $file_names_pointer=sub_common::files_list($variables{dir_bowtie2}, 'file_name');
	my @file_names=@$file_names_pointer;
	
	#frame: select bowtie index
	my @fa_names=grep(/\.fa$|\.fasta$/, @file_names);
	my($frame_1, $cb_fasta_file)=sub_gui::combobox_frame('Select genome sequences (.fa file)', \@fa_names, 0);
	$table->attach_defaults($frame_1, 0,4,0,1);
	
	#frame: gtf file
	my @gtf_file_names=grep(/\.gtf$|\.gff$/, @file_names);
	my($frame_2, $cb_gtf_file)=sub_gui::combobox_frame('Select annotation file (.gtf file)', \@gtf_file_names, 0);
	$table->attach_defaults($frame_2, 0,4,1,2);
	
	$frame = Gtk2::Frame->new();
	$frame->set_border_width(5);
		my $button=Gtk2::Button->new('Save and close');
		$button->set_border_width(10);
		$button->signal_connect('clicked' => sub { 
				my $genome_fasta_file=$cb_fasta_file->get_active_text();
				$variables{genome_fasta_file}=$variables{dir_bowtie2}.$genome_fasta_file;
				$variables{genome_index_name}=$genome_fasta_file;
				$variables{genome_index_name}=~s/\.fa$|\.fasta$//;
				$variables{genome_index}=$variables{dir_bowtie2}.$variables{genome_index_name};
				$variables{genome_annotation_file}=$variables{dir_bowtie2}.$cb_gtf_file->get_active_text();
				#
				sub_basic::refresh_log($variables{file_var}, 'genome_index', $variables{genome_index} );
				sub_basic::refresh_log($variables{file_var}, 'genome_index_name', $variables{genome_index_name});
				sub_basic::refresh_log($variables{file_var}, 'genome_annotation_file', $variables{genome_annotation_file});
				sub_basic::refresh_log($variables{file_var}, 'genome_fasta_file', $variables{genome_fasta_file});
				
				$window->destroy;
			} );
	$frame->add($button);
	$table->attach_defaults($frame,0,4,2,3);
	
	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}


#######################
sub nb_mRNA_S2_spliced_alignment{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Options of spliced aligner (TopHat v.2)', '1200','800');
	my $table=Gtk2::Table->new(11,6,TRUE);
	
		#frame: select sample names files
		$frame=Gtk2::Frame->new('Select sample names');
		$frame->set_border_width(5);
			#get sample names
			my $pointer=sub_basic::SM_sample_info(\%variables);
			%variables=%$pointer;
			my @sample_names=split(',', $variables{sample_names});
			#print @sample_names;
			my $sub_table=items_selection(\@sample_names, $variables{file_var}, 'tophat_sample_names');
		$frame->add($sub_table);
	$table->attach_defaults($frame, 0,6,0,6);
		
	#frame: mismatch
	my @entry_text=(0..5);
	my($frame_1, $cb_read_mismatches)=sub_gui::combobox_frame('Mismatches of read alignments', \@entry_text, 2);
	$table->attach_defaults($frame_1, 0,2,6,7);
	
	#frame: read gap
	@entry_text=(0..4);
	my($frame_2, $cb_read_gap_length)=sub_gui::combobox_frame('Length of discarded gaps of read alignments', \@entry_text, 2);
	$table->attach_defaults($frame_2, 2,4,6,7);
	
	#frame: read gap
	@entry_text=(0..4);
	my($frame_3, $cb_read_edit_dist)=sub_gui::combobox_frame('Edit distance of read alignments', \@entry_text, 2);
	$table->attach_defaults($frame_3, 4,6,6,7);
	
	#frame: mate dist
	@entry_text=(1..8);
	my($frame_4, $cb_mate_inner_dist)=sub_gui::combobox_frame('Inner distance between mate pairs', \@entry_text, 1);
	$table->attach_defaults($frame_4, 0,2,7,8);
	
	#frame: mate
	@entry_text=(10, 20, 30, 40);
	my($frame_5, $cb_mate_std_dev)=sub_gui::combobox_frame('Standard deviation for the distribution on inner distances', \@entry_text, 1);
	$table->attach_defaults($frame_5, 2,4,7,8);
	
	#frame: anchor length
	@entry_text=(3..12);
	my($frame_6, $cb_min_anchor_length)=sub_gui::combobox_frame('Anchor length', \@entry_text, 5);
	$table->attach_defaults($frame_6, 4,6,7,8);
	
	#frame:
	@entry_text=(0..2);
	my($frame_7, $cb_splice_mismatches)=sub_gui::combobox_frame('Splice mismatches', \@entry_text, 0);
	$table->attach_defaults($frame_7, 0,2,8,9);
	
	#frame:
	@entry_text=(1..10);
	my($frame_8, $cb_min_intron_length)=sub_gui::combobox_frame('Minimum of intron length', \@entry_text, 6);
	$table->attach_defaults($frame_8, 2,4,8,9);
	
	#frame:
	@entry_text=map {$_*50000} (1..20);
	my($frame_9, $cb_max_intron_length)=sub_gui::combobox_frame('Maximum of intron length', \@entry_text, 9);
	$table->attach_defaults($frame_9, 4,6,8,9);
	
	#frame:
	@entry_text=(1..6);
	my($frame_10, $cb_max_insertion_length)=sub_gui::combobox_frame('Maximum of insertion length', \@entry_text, 2);
	$table->attach_defaults($frame_10, 0,2,9,10);
	
	#frame:
	@entry_text=(1..6);
	my($frame_11, $cb_max_deletion_length)=sub_gui::combobox_frame('Maximum of deletion length', \@entry_text, 2);
	$table->attach_defaults($frame_11, 2,4,9,10);

	#frame:
	@entry_text=(10..20);
	my($frame_12, $cb_integer_quals)=sub_gui::combobox_frame('Quality values', \@entry_text, 3);
	$table->attach_defaults($frame_12, 4,6,9,10);
	
	#frame:
	@entry_text=('fr-unstranded', 'fr-firststrand ', 'fr-secondstrand');
	my($frame_13, $cb_strand_specific)=sub_gui::combobox_frame('Strand specific of libraries', \@entry_text, 0);
	$table->attach_defaults($frame_13, 0,2,10,11);
	
	#Frame
	my $read_realign_edit_dist=$cb_read_edit_dist->get_active_text()-1;
	$read_realign_edit_dist=0 if $read_realign_edit_dist<0;
	$frame = Gtk2::Frame->new();
	$frame->set_border_width(5);
	my $button=Gtk2::Button->new('Save and close');
	$button->signal_connect(clicked=>sub{
		$variables{tophat_annotation_dir}=$variables{dir_result}.'genome_annotation';
		$variables{tophat_options}=join (" ", '-p 1', '-N', $cb_read_mismatches->get_active_text(), 
								'--read-gap-length', $cb_read_gap_length->get_active_text(),
								'--read-edit-dist', $cb_read_edit_dist->get_active_text(),
								'--read-realign-edit-dist', $read_realign_edit_dist,
								'-r', $cb_mate_inner_dist->get_active_text(),
								'--mate-std-dev', $cb_mate_std_dev->get_active_text(),
								'-a', $cb_min_anchor_length->get_active_text(),
								'-m', $cb_splice_mismatches->get_active_text(),
								'-i', $cb_min_intron_length->get_active_text(),
								'-I', $cb_max_intron_length->get_active_text(),
								'--max-insertion-length', $cb_max_insertion_length->get_active_text(),
								'--max-deletion-length', $cb_max_deletion_length->get_active_text(),
								'--library-type', $cb_strand_specific->get_active_text(),
								#'--integer-quals', $cb_integer_quals->get_active_text(),
							);
		#default directory and files
		$variables{library_type}=$cb_strand_specific->get_active_text();
		$variables{mRNA_genome_mapping}='yes';
		$variables{mRNA_sequencing_control}='yes';
		sub_basic::refresh_log($variables{file_var}, 'tophat_options', $variables{tophat_options} );
		sub_basic::refresh_log($variables{file_var}, 'mRNA_genome_mapping', $variables{mRNA_genome_mapping});
		sub_basic::refresh_log($variables{file_var}, 'mRNA_sequencing_control', $variables{mRNA_sequencing_control});
		#print "$variables{tophat_output_dir}\n", "$variables{tophat_options}\n";
		$window->destroy;
	});
	$button->set_border_width(5);
	$frame->add($button);
	$table->attach_defaults($frame,1,5,11,12);
	
	$table->show_all();
	$window->add($table);
	$window->show_all;
	Gtk2->main();
}

#######################
sub nb_mRNA_S3_transcript_assembling{
	#get index information
	my $file_names_pointer=sub_common::files_list($variables{dir_bowtie2}, 'file_name');
	my @file_names=@$file_names_pointer;
	my @fa_names=grep(/\.fa$/, @file_names);
	my @gtf_file_names=grep(/\.gtf$/, @file_names);
		
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Options of transcripts assembling (Cufflinks)', '1200','600');
	my $table=Gtk2::Table->new(8,6,TRUE);
	
	#frame: --max-bundle-frags
	my @entry_text=(1e4,5e4,1e5,5e5,1e6,5e6);
	my($frame_1, $cb_max_bundle_frags)=sub_gui::combobox_frame('Number of fragments of a locus before being skipped', \@entry_text, 4);
	$table->attach_defaults($frame_1, 0,2,1,2);

	#frame
	@entry_text=('--total-hits-norm', '--compatible-hits-norm');
	my($frame_2, $cb_hits_counting)=sub_gui::combobox_frame('Alignment hit counting', \@entry_text, 1);
	$table->attach_defaults($frame_2, 2,4,1,2);
	
	#frame
	@entry_text=('length-correction', '--no-effective-length-correction', '--no-length-correction');
	my($frame_3, $cb_length_correction)=sub_gui::combobox_frame('Length correction for transcript FPKM', \@entry_text, 0);
	$table->attach_defaults($frame_3, 4,6,1,2);
	
	#-F/--min-isoform-fraction
	@entry_text=map {$_*0.05} (0..20);
	my($frame_4, $cb_isoform_fraction)=sub_gui::combobox_frame('Suppress low abundant isoforms for a gene', \@entry_text, 1);
	$table->attach_defaults($frame_4, 0,2,2,3);
	
	#-j/--pre-mrna-fraction
	@entry_text=map {$_*0.05} (0..20);
	my($frame_5, $cb_intra_intronic_fraction)=sub_gui::combobox_frame('Suppress low abundant intra-intronic transcripts', \@entry_text, 1);
	$table->attach_defaults($frame_5, 2,4,2,3);
	
	#-I/--max-intron-length
	@entry_text=map {$_*1e5} (1..10);
	my($frame_6, $cb_intron_length)=sub_gui::combobox_frame('Maximum intron length', \@entry_text, 2);
	$table->attach_defaults($frame_6, 4,6,2,3);
	
	#-a/--junc-alpha
	@entry_text=(0.0001,0.001,0.01,0.1,0.2,0.5,1);
	my($frame_7, $cb_alpha_value)=sub_gui::combobox_frame('The alpha value of binomial test', \@entry_text, 1);
	$table->attach_defaults($frame_7, 0,2,3,4);
	
	#-A/--small-anchor-fraction
	@entry_text=(0.01,0.05,0.09,0.2,0.5,1);
	my($frame_8, $cb_anchor_fraction)=sub_gui::combobox_frame("Percent read overhang taken as 'suspiciously small' ", \@entry_text, 2);
	$table->attach_defaults($frame_8, 2,4,3,4);
	
	#--min-frags-per-transfrage
	@entry_text=(1,3,5,10,20);
	my($frame_9, $cb_min_frags)=sub_gui::combobox_frame("Minimum number of fragments needed for new transfrags", \@entry_text, 3);
	$table->attach_defaults($frame_9, 4,6,3,4);
	
	#--overhang-tolerance
	@entry_text=(3..20);
	my($frame_10, $cb_overhang_tolerance)=sub_gui::combobox_frame('Maximum basepairs of a gap size to fill between transfrags', \@entry_text, 5);
	$table->attach_defaults($frame_10, 0,2,4,5);
	
	#--max-bundle-length
	@entry_text=map {$_*5e5} (1..20);
	my($frame_11, $cb_max_bundle_length)=sub_gui::combobox_frame('Maximum genomic length allowed for a given bundle', \@entry_text, 6);
	$table->attach_defaults($frame_11, 2,4,4,5);
	
	#--min-intron-length
	@entry_text=map {$_*10} (1..10);
	my($frame_12, $cb_min_intron_length)=sub_gui::combobox_frame('Minimum intron size allowed in genome', \@entry_text, 4);
	$table->attach_defaults($frame_12, 4,6,4,5);
	
	#--max-multiread-fraction
	@entry_text=map {$_*0.05} (10..20);
	my($frame_13, $cb_multiread_fraction)=sub_gui::combobox_frame('Maximum fraction of allowed multireads per transcript', \@entry_text, 5);
	$table->attach_defaults($frame_13, 0,2,5,6);
	
	#--overlap-radius
	@entry_text=map {$_*10} (1..10);
	my($frame_14, $cb_overlap_radius)=sub_gui::combobox_frame('Maximum gap size to fill between transfrags (in bp)', \@entry_text, 4);
	$table->attach_defaults($frame_14, 2,4,5,6);
	
	#frame: --max-mle-iterations
	@entry_text=map {$_*1e3} (5..10);
	my($frame_15, $cb_max_iterations)=sub_gui::combobox_frame('Number of iterations during maximum likelihood estimation', \@entry_text, 0);
	$table->attach_defaults($frame_15, 4,6,5,6);
	
	#-u/--multi-read-correct
	my $check_button_1 = Gtk2::CheckButton->new ("Make an initial procedure to more accurately weight reads mapping to multiple locations");
	$check_button_1->set_active (FALSE);
	$table->attach_defaults($check_button_1, 0,6,6,7);

	
	#Frame
	$frame = Gtk2::Frame->new();
	$frame->set_border_width(5);
	my $button=Gtk2::Button->new('Save and close');
	$button->signal_connect(clicked=>sub{
		$variables{cufflinks_options}=join (" ", '-p 1', '--library-type', $variables{library_type},
									'--max-mle-iterations', $cb_max_iterations->get_active_text(), 
									'--max-bundle-frags', $cb_max_bundle_frags->get_active_text(), 
									'-F', $cb_isoform_fraction->get_active_text(), 
									'-j', $cb_intra_intronic_fraction->get_active_text(), 
									'-I', $cb_intron_length->get_active_text(), 
									'-a', $cb_alpha_value->get_active_text(), 
									'-A', $cb_anchor_fraction->get_active_text(), 
									'--min-frags-per-transfrag', $cb_min_frags->get_active_text(), 
									'--overhang-tolerance', $cb_overhang_tolerance->get_active_text(),
									'--max-bundle-length', $cb_max_bundle_length->get_active_text(),
									'--min-intron-length', $cb_min_intron_length->get_active_text(),
									'--max-multiread-fraction', $cb_multiread_fraction->get_active_text(),
									'--overlap-radius', $cb_overlap_radius->get_active_text(),
									$cb_hits_counting->get_active_text(),
		);
		$variables{cufflinks_options}=join(' ', $variables{cufflinks_options}, $cb_length_correction->get_active_text() ) unless $cb_length_correction->get_active_text() eq 'length-correction';
		$variables{cufflinks_options}=join(' ', $variables{cufflinks_options}, '--multi-read-correct' ) if $check_button_1->get_active;
		sub_basic::refresh_log($variables{file_var}, 'mRNA_transcripts_assembling', 'yes' );
		sub_basic::refresh_log($variables{file_var}, 'cufflinks_options', $variables{cufflinks_options} );
		
		#cuffmerge
		sub_basic::refresh_log($variables{file_var}, 'mRNA_cuffmerge', 'yes' );
		
		$window->destroy;
	});
	$frame->add($button);
	$table->attach_defaults($frame,1,5,7,8);
	
	$table->show_all();
	$window->add($table);
	$window->show_all;
	Gtk2->main();
}
#######################
sub nb_mRNA_S4_differential_transcripts{
	
	#get index information
	my $file_names_pointer=sub_common::files_list($variables{dir_bowtie2}, 'file_name');
	my @file_names=@$file_names_pointer;
	my @fa_names=grep(/\.fa$/, @file_names);
	my @gtf_file_names=grep(/\.gtf$/, @file_names);
	
	print  "################\n";
	my @DEG_candidates;
	#get sample names
	my $pointer=sub_basic::SM_sample_info(\%variables);
	%variables=%$pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	@DEG_candidates=sort (keys %sample_info);
	#print "@DEG_candidates\n";
	
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Options of differential expression analysis (Cuffdiff)', '1200','700');
	my $table=Gtk2::Table->new(10,6,TRUE);
	
	#items selection
		my $frame=Gtk2::Frame->new('Select mRNA samples for differential expression analysis ');
			my $sub_table=multiple_items_selection(\@DEG_candidates, $variables{file_var}, 'diff_A', 'diff_B' );
		$frame->add($sub_table);
	$table->attach_defaults($frame, 0,6,0,7);
	
		#frame: fasta files
		my($frame_1, $entry_1)=sub_gui::entry_frame('Name of differential expression analysis', 'cuffdiff');
	$table->attach_defaults($frame_1, 0,2,7,8);
	
		#frame:--num-frag-count-draws
		my @entry_text=(50,100,200,300,400,500);
		my($frame_2, $cb_num_frag_count_draws)=sub_gui::combobox_frame('Number of fragment generation samples', \@entry_text, 1);
	$table->attach_defaults($frame_2, 2,4,7,8);
	
		#frame:-c/--min-alignment
		@entry_text=(5,10,15,20,30);
		my($frame_3, $cb_min_alignments)=sub_gui::combobox_frame('Minimum number of alignments', \@entry_text, 1);
	$table->attach_defaults($frame_3, 4,6,7,8);
	
		#frame:--num-frag-assign-draws
		@entry_text=(20,50,100,200,300,400,500);
		my($frame_4, $cb_num_frag_assign_draws)=sub_gui::combobox_frame('Number of fragment assignment samples per generation', \@entry_text, 1);
	$table->attach_defaults($frame_4, 0,2,8,9);
	
		#frame:--min-reps-for-js-test
		@entry_text=(1,2,3,4,5,6);
		my($frame_5, $cb_min_reps)=sub_gui::combobox_frame('Replicates needed for relative isoform shift testing', \@entry_text, 2);
	$table->attach_defaults($frame_5, 2,4,8,9);
	
		#frame:--FDR
		@entry_text=(0.005,0.01,0.05,0.1);
		my($frame_6, $cb_FDR)=sub_gui::combobox_frame('False discovery rate allowed', \@entry_text, 2);
	$table->attach_defaults($frame_6, 4,6,8,9);
	
		#frame:--max-mle-iterations
		@entry_text=(1e3,5e3,1e4,5e4);
		my($frame_7, $cb_max_mle_iterations)=sub_gui::combobox_frame('Maximum iterations allowed for MLE calculation', \@entry_text, 1);
	$table->attach_defaults($frame_7, 0,2,9,10);
	
		#frame:--max-bundle-frages
		@entry_text=(1e4,1e5,1e6,1e7,1e8);
		my($frame_8, $cb_max_bundle_frags)=sub_gui::combobox_frame('Maximum fragments allowed in a bundle before skipping', \@entry_text, 2);
	$table->attach_defaults($frame_8, 2,4,9,10);
	#############
	#Frame
	$frame = Gtk2::Frame->new();
	my $button=Gtk2::Button->new('Save and close');
	$button->signal_connect(clicked=>sub{
		$variables{cuffmerge_output_dir}=$variables{dir_result}.'cuffmerge';
		$variables{cuffdiff_output_dir}=$variables{dir_result}.$entry_1->get_text().'_'.time;
		$variables{cuffdiff_options}=join (" ", '-p', $variables{threads_num}, 
									'-u', $variables{cuffmerge_output_dir}.'/merged.gtf',
									'-c', $cb_min_alignments->get_active_text(), 
									'--FDR', $cb_FDR->get_active_text(), 
									'--max-mle-iterations', $cb_max_mle_iterations->get_active_text(),
									'--max-bundle-frags', $cb_max_bundle_frags->get_active_text(), 
									'--num-frag-count-draws', $cb_num_frag_count_draws->get_active_text(), 
									'--num-frag-assign-draws', $cb_num_frag_assign_draws->get_active_text(), 
									'--min-reps-for-js-test', $cb_min_reps->get_active_text(),
									'-o', $variables{cuffdiff_output_dir},
									'-b', $variables{genome_fasta_file}, 		);
		sub_basic::refresh_log($variables{file_var}, 'diff_candidates', join(',', @DEG_candidates) );
		sub_basic::refresh_log($variables{file_var}, 'cuffdiff_options', $variables{cuffdiff_options} );
		sub_basic::refresh_log($variables{file_var}, 'cuffdiff_output_dir', $variables{cuffdiff_output_dir} );
		sub_basic::refresh_log($variables{file_var}, 'mRNA_cuffdiff', 'yes' );
		$window->destroy;
	});
	$frame->add($button);
	$table->attach_defaults($frame,1,5,10,11);
	
	$table->show_all();
	$window->add($table);
	$window->show_all;
	Gtk2->main();
}

######################################################################
sub nb_targets{
	my $scroller = Gtk2::ScrolledWindow->new();
	$scroller->set_shadow_type ('etched-out');
	$scroller->set_policy ('automatic', 'automatic');
	$scroller->set_size_request (200, 200);

	my $table_stock=Gtk2::Table->new(6, 10, TRUE); #row_number, and column number 
	$table_stock->set_border_width(20);
	$table_stock->set_row_spacings(20);
	#
	$table_stock->attach_defaults(&sub_gui::nb_button('Data pre-processing', sub{ 
			my $perl_script=$perl_dir."/eRNA_expression_data_processing.pl";
			system("perl $perl_script $variables{file_var}");
		}), 0,1,0,2);#column, row
		my $arrow = Gtk2::Arrow->new('right', 'etched-out');
	$table_stock->attach_defaults($arrow, 1,2,1,2);
	
	$table_stock->attach_defaults(&sub_gui::nb_button('Differential expression profiles analysis', sub{ 
			my $perl_script=$perl_dir."/eRNA_differential_expression_analysis.pl";
			system("perl $perl_script $variables{file_var}");
		}), 2,4,0,1);
	$table_stock->attach_defaults(&sub_gui::nb_button('Gene significance through incrusive analysis', sub{ 
			my $perl_script=$perl_dir."/eRNA_gene_significance_analysis.pl";
			system("perl $perl_script $variables{file_var}");
		}), 2,4,1,2);
	$table_stock->attach_defaults(&sub_gui::nb_button('Weighted Gene Co-Expression Network Analysis (WGCNA)', sub{ 
			my $perl_script=$perl_dir."/eRNA_WGCNA_analysis.pl";
			system("perl $perl_script $variables{file_var}");
		}), 2,4,2,3);
		$arrow = Gtk2::Arrow->new('right', 'etched-out');
	$table_stock->attach_defaults($arrow, 4,5,1,2);
	
	$scroller->add_with_viewport($table_stock);
	return($scroller);
}

#################################
sub nb_targets_setup{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Setup of target screening', '600','300');
	my $table = Gtk2::Table->new(4,4,TRUE);
		
		#frame: methods
		my @entry_text=('1: Identification of differential expressed genes using DEseq', '2: Prediction of crucial genes using RandomForest');
		my($frame_1, $cb_targets_method)=sub_gui::combobox_frame('Methods of targets screening', \@entry_text, 0);
	$table->attach_defaults($frame_1,0,4,0,1);
	
		#frame: data frame file
		$frame=Gtk2::Frame->new('Transcriptional level file');
			my $input_chooser =Gtk2::FileChooserButton->new ('select a file' , 'open');
			$input_chooser->set_filename($variables{dir_result});
		$frame->add($input_chooser);
	$table->attach_defaults($frame, 0,2,1,2);
		#frame: sample_info file
		$frame=Gtk2::Frame->new('Sample traits file');
			my $sample_info_chooser =Gtk2::FileChooserButton->new ('select a file' , 'open');
			$sample_info_chooser->set_filename($variables{dir_result});
		$frame->add($sample_info_chooser);
	$table->attach_defaults($frame, 2,4,1,2);

		#Frame 5
		$frame = Gtk2::Frame->new();
		my $button_2=Gtk2::Button->new('Save and close');
		$button_2->set_border_width(5);
		$button_2->signal_connect(clicked=>sub{	
			my $targets_method=$cb_targets_method->get_active_text();
			if($targets_method=~/^1/){	$variables{targets_method}='DEseq';	}
			#elsif($targets_method=~/^2/){	$variables{targets_method}='cummeRbund';	}
			elsif($targets_method=~/^2/){	$variables{targets_method}='RF';	}
			sub_basic::refresh_log($variables{file_var}, 'R_targets_method', $variables{targets_method});
			#dat frame file
			if( -f $input_chooser->get_filename){
				$variables{targets_df_file}=$input_chooser->get_filename; 
			}
			sub_rna::refresh_R_script($variables{R_DEseq_file}, 'R_df_file', $variables{targets_df_file});
			sub_rna::refresh_R_script($variables{R_RF_file}, 'R_df_file', $variables{targets_df_file});
			
			#cover the value set in the window named nb_self_sample
			if (-f $sample_info_chooser->get_filename){
				$variables{file_sample_info}=$sample_info_chooser->get_filename; 
				sub_basic::refresh_log($variables{file_var}, 'sample_info_file', $variables{file_sample_info});
			}
			$window->destroy;
		});
		$frame->add($button_2) ;
		$table->attach_defaults($frame,1,3,3,4);

	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}
#################################
sub nb_targets_DEG{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Differential expression profiling analysis', '1000','700');
	my $table = Gtk2::Table->new(11,4,TRUE);
		
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Name of group A'); #group A name
			my $entry_A=Gtk2::Entry->new();
			$entry_A->set_text('A');
		$frame->add($entry_A);
	$table->attach_defaults($frame,0,2,0,1);
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Name of group B'); # name of group B
			my $entry_B=Gtk2::Entry->new();
			$entry_B->set_text('B');
		$frame->add($entry_B);
	$table->attach_defaults($frame,2,4,0,1);
		#frame: p value cutoff
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Cutoff of p value');
			my $cb_p_value=Gtk2::ComboBox->new_text();
			foreach(0.1,0.05,0.01,0.001){
				$cb_p_value->append_text($_);
			}
			$cb_p_value->set_active(1);
		$frame->add($cb_p_value);
	$table->attach_defaults($frame,0,2,1,2);
			#frame: transcriptional level cutoff
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Cutoff of transcriptional level on average');
			my $cb_read_cutoff=Gtk2::ComboBox->new_text();
			foreach(1..10){
				$cb_read_cutoff->append_text($_);
			}
			$cb_read_cutoff->set_active(4);
		$frame->add($cb_read_cutoff);
	$table->attach_defaults($frame,2,4,1,2);
		#Frame: list all sample traits
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Sample Traits');
			my $matrix_pointer=sub_rna::read_sample_info($variables{file_sample_info}, 'attr');
			my $box_sample_traits=CellRenderer_window($matrix_pointer);
		$frame->add($box_sample_traits);
	$table->attach_defaults($frame,0,6,2,5);
	
		#display selected samples
		$frame = Gtk2::Frame->new();
		$frame->set_label('Checking information');
			my $textview = Gtk2::TextView->new();
			$textview->set_wrap_mode('GTK_WRAP_WORD_CHAR');
		$frame->add($textview);
	$table->attach_defaults($frame,0,6,5,8);
		#Frame: list all sample traits
		$frame = Gtk2::Frame->new();
		$frame->set_label('Formula of group A for differential comparison');
			my $entry_group_A=Gtk2::Entry->new();
		$frame->add($entry_group_A);
	$table->attach_defaults($frame,0,5,8,9);
		$frame = Gtk2::Frame->new();
		$frame->set_label('Formula of group B for differential comparison');
			my $entry_group_B=Gtk2::Entry->new();
		$frame->add($entry_group_B);
	$table->attach_defaults($frame,0,5,9,10);
		#display button
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		my $button_display=Gtk2::Button->new('Display');
		$button_display->signal_connect(clicked=>sub{
			my @cond_A=split(";", $entry_group_A->get_text());
			my @cond_B=split(";", $entry_group_B->get_text());
			my $sample_info_pointer=sub_rna::read_sample_info($variables{file_sample_info}, 'hash');
			my %sample_info=%$sample_info_pointer;
			my %group_A= %sample_info;
			my %group_B= %sample_info;
			foreach my $sample_name(keys %sample_info){
				foreach(@cond_A){
					my($col_name, $cond_value)=split(":", $_);
					my $sample_value=$sample_info{$sample_name}->{$col_name};
					if($cond_value=~/-/){ #numeric value
						my ($value_min, $value_max)=split("-", $cond_value);
						delete $group_A{$sample_name} if $sample_value<$value_min or $sample_value>$value_max;
					}
					else{
						delete $group_A{$sample_name} unless $sample_value eq $cond_value;
					}
				}
				foreach(@cond_B){
					my($col_name, $cond_value)=split(":", $_);
					my $sample_value=$sample_info{$sample_name}->{$col_name};
					if($cond_value=~/\-/){ #numeric value
						my ($value_min, $value_max)=split(/\-/, $cond_value);
						delete $group_B{$sample_name} if $sample_value<$value_min or $sample_value>$value_max;
					}
					else{
						delete $group_B{$sample_name} unless $sample_value eq $cond_value;
					}
				}
			}
			$variables{R_group_A}=join(",", keys %group_A);
			$variables{R_group_B}=join(",", keys %group_B);
			my $A_num=split(",", $variables{R_group_A});
			my $B_num=split(",", $variables{R_group_B});
			#display
			my $textbuffer=Gtk2::TextBuffer->new();
			my $iter=$textbuffer->get_iter_at_offset(0);
			$textbuffer->insert($iter, "Selected samples in group A ($A_num): $variables{R_group_A}\n");
			$textbuffer->insert($iter, "\nSelected samples in group B ($B_num): $variables{R_group_B}\n");
			$textview->set_buffer($textbuffer);
			$textview->show_all();
		});
		$frame->add($button_display);
	$table->attach_defaults($frame,5,6,8,10);
	
		#Frame 5
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		my $hbox=Gtk2::HBox->new();
			my $button_1=Gtk2::Button->new('Clear');
			$button_1->set_border_width(5);
			$button_1->signal_connect(clicked=>sub{
				delete $variables{R_group_A};
				delete $variables{R_group_B};
			});
		$hbox->add($button_1);
		my $button_2=Gtk2::Button->new('Save and close');
		$button_2->set_border_width(5);
		$button_2->signal_connect(clicked=>sub{
			#refresh names of A and B 
			my $A_name=$entry_A->get_text();
			my $B_name=$entry_B->get_text();
			sub_rna::refresh_R_script($variables{R_DEseq_file}, 'R_A_name', $A_name);
			sub_rna::refresh_R_script($variables{R_DEseq_file}, 'R_B_name', $B_name);
			#refresh samples A and B
			sub_rna::refresh_R_script($variables{R_DEseq_file}, 'R_group_A', $variables{R_group_A});
			sub_rna::refresh_R_script($variables{R_DEseq_file}, 'R_group_B', $variables{R_group_B});
			#refresh files
			
			sub_rna::refresh_R_script($variables{R_DEseq_file}, 'R_diff_file', $variables{dir_result}.'DEG_diff_'.$A_name.'_'.$B_name.'.txt');
			sub_rna::refresh_R_script($variables{R_DEseq_file}, 'R_scatterplot_file', $variables{dir_result}.'DEG_scatterplot_'.$A_name.'_'.$B_name.'.jpg');
			sub_rna::refresh_R_script($variables{R_DEseq_file}, 'R_boxplot_file', $variables{dir_result}.'DEG_boxplot_'.$A_name.'_'.$B_name.'.jpg');
			sub_rna::refresh_R_script($variables{R_DEseq_file}, 'R_p_value', $cb_p_value->get_active_text());
			sub_rna::refresh_R_script($variables{R_DEseq_file}, 'R_read_cutoff', $cb_read_cutoff->get_active_text());
			
			if(exists $variables{R_group_A} and exists $variables{R_group_B}){
				sub_basic::refresh_log($variables{file_var}, 'R_DEG_analysis', 'yes');
				$window->destroy;
			}
			else{	
				sub_gui::popup_window('Error', 'Samples used for comparison should be selected before running!');
			}
		});
		$hbox->add($button_2);
		$frame->add($hbox) ;
		$table->attach_defaults($frame,2,4,10,11);

	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}
#################################
sub nb_targets_RF{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Recursive partitioning analysis of expression', '1000','700');
	my $table = Gtk2::Table->new(11,4,TRUE);
		
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Name of RandomForest classification'); 
			my $entry_RF_name=Gtk2::Entry->new();
			$entry_RF_name->set_text('A');
		$frame->add($entry_RF_name);
	$table->attach_defaults($frame,0,2,0,1);
			#frame: transcriptional level cutoff
		$frame = Gtk2::Frame->new();
		$frame->set_label('Cutoff of transcriptional level on average');
			my $cb_read_cutoff=Gtk2::ComboBox->new_text();
			foreach(1..10){
				$cb_read_cutoff->append_text($_);
			}
			$cb_read_cutoff->set_active(4);
		$frame->add($cb_read_cutoff);
	$table->attach_defaults($frame,2,4,0,1);
		#frame: tree number
		$frame = Gtk2::Frame->new();
		$frame->set_label('Overall number of trees in the forest');
			my $cb_tree_num=Gtk2::ComboBox->new_text();
			foreach(500,2000,5000,10000){
				$cb_tree_num->append_text($_);
			}
			$cb_tree_num->set_active(2);
		$frame->add($cb_tree_num);
	$table->attach_defaults($frame,4,6,0,1);
		#frame: number of input variables
		$frame = Gtk2::Frame->new();
		$frame->set_label('Number of randomly preselected variables for each split');
			my $cb_variables_num=Gtk2::ComboBox->new_text();
			foreach(3,5,7,9,10){
				$cb_variables_num->append_text($_);
			}
			$cb_variables_num->set_active(1);
		$frame->add($cb_variables_num);
	$table->attach_defaults($frame,0,2,1,2);
		#frame: Symbolic description variables to be fit in the RF modeling
		$frame = Gtk2::Frame->new();
		$frame->set_label('Dependent variables fited in RF modeling');
			my @sample_traits=split(";", sub_rna::file_info($variables{file_sample_info}, 'first_line') );
			my $sub_table=multiple_items_selection(\@sample_traits, $variables{file_var}, 'targets_RF_dependent', 'targets_RF_independent');
		$frame->add($sub_table);
		$table->attach_defaults($frame,0,6,2,4);
		
		#Frame: list all sample traits
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Sample Traits');
			my $attr_matrix_pointer=sub_rna::read_sample_info($variables{file_sample_info}, 'attr');
			my $box_sample_traits=CellRenderer_window($attr_matrix_pointer);
		$frame->add($box_sample_traits);
	$table->attach_defaults($frame,0,6,4,7);
	
		#display selected samples
		$frame = Gtk2::Frame->new();
		$frame->set_label('Checking information');
			my $textview = Gtk2::TextView->new();
			$textview->set_wrap_mode('GTK_WRAP_WORD_CHAR');
		$frame->add($textview);
	$table->attach_defaults($frame,0,6,7,9);
		#Frame: list all sample traits
		$frame = Gtk2::Frame->new();
		$frame->set_label('Formular of sample selection');
			my $entry_RF_samples=Gtk2::Entry->new();
		$frame->add($entry_RF_samples);
	$table->attach_defaults($frame,0,4,9,10);
		#display button
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		my $button_display=Gtk2::Button->new('Display');
		$button_display->signal_connect(clicked=>sub{
			my @cond_samples=split(";", $entry_RF_samples->get_text());
			my $sample_info_pointer=sub_rna::read_sample_info($variables{file_sample_info}, 'hash');
			my %sample_info=%$sample_info_pointer;
			my %RF_samples= %sample_info;
			foreach my $sample_name(keys %sample_info){
				foreach(@cond_samples){
					my($col_name, $cond_value)=split(":", $_);
					my $sample_value=$sample_info{$sample_name}->{$col_name};
					if($cond_value=~/-/){ #numeric value
						my ($value_min, $value_max)=split("-", $cond_value);
						delete $RF_samples{$sample_name} if $sample_value<$value_min or $sample_value>$value_max;
					}
					else{
						delete $RF_samples{$sample_name} unless $sample_value eq $cond_value;
					}
				}
			}
			$variables{R_RF_samples}=join(",", keys %RF_samples);
			my $RF_num=split(',', $variables{R_RF_samples});
			#display
			my $textbuffer=Gtk2::TextBuffer->new();
			my $iter=$textbuffer->get_iter_at_offset(0);
			$textbuffer->insert($iter, "Selected samples for classification ($RF_num): $variables{R_RF_samples}\n");
			$textbuffer->insert($iter, "\nSelected input variables for classification: $variables{R_RF_variables}\n");
			$textview->set_buffer($textbuffer);
			$textview->show_all();
		});
		$frame->add($button_display);
	$table->attach_defaults($frame,4,6,9,10);
	
		#Frame 5
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		my $hbox=Gtk2::HBox->new();
			my $button_1=Gtk2::Button->new('Clear');
			$button_1->set_border_width(5);
			$button_1->signal_connect(clicked=>sub{
				delete $variables{R_RF_samples};
				delete $variables{R_RF_variables};
			});
		$hbox->add($button_1);
		my $button_2=Gtk2::Button->new('Save and close');
		$button_2->set_border_width(5);
		$button_2->signal_connect(clicked=>sub{
			#refresh names of RF
			my $RF_name=$entry_RF_name->get_text();
			sub_rna::refresh_R_script($variables{R_RF_file}, 'R_RF_name', $RF_name);
			#refresh selected samples
			sub_rna::refresh_R_script($variables{R_RF_file}, 'R_RF_samples', $variables{R_RF_samples});
			
			#refresh files
			sub_rna::refresh_R_script($variables{R_RF_file}, 'R_sample_traits_file', $variables{file_sample_info});
			sub_rna::refresh_R_script($variables{R_RF_file}, 'R_importance_file', $variables{dir_result}.'Prediction_importance'.$RF_name.'.txt');
			sub_rna::refresh_R_script($variables{R_RF_file}, 'R_mining_plot_file', $variables{dir_result}.'Prediction_mining_'.$RF_name.'.jpg');
			#refresh other variables
			sub_rna::refresh_R_script($variables{R_RF_file}, 'R_read_cutoff', $cb_read_cutoff->get_active_text());
			sub_rna::refresh_R_script($variables{R_RF_file}, 'R_tree_num', $cb_tree_num->get_active_text());
			sub_rna::refresh_R_script($variables{R_RF_file}, 'R_mtry_num', $cb_variables_num->get_active_text());
			#subroutine: refresh factors
			my $dependent_var=sub_basic::read_log($variables{file_var}, 'targets_RF_dependent');
			$dependent_var='.' if $dependent_var eq 'NA';
			$dependent_var=~s/,/+/g;
			my $independent_var=sub_basic::read_log($variables{file_var}, 'targets_RF_independent');
			$independent_var='.' if $independent_var eq 'NA';
			$independent_var=~s/,/+/g;
			my $new_line='data.cforest<-cforest('.$dependent_var.'~'.$independent_var.', data=countTable_Sig, controls=data.controls)';
			sub_rna::refresh_cforest_line($variables{R_RF_file}, '\<\-cforest\(', $new_line);
			
			if(exists $variables{R_RF_samples} and $dependent_var cmp 'NA' ){
				sub_basic::refresh_log($variables{file_var}, 'R_RF_analysis', 'yes');
				$window->destroy;
			}
			else{	
				sub_gui::popup_window('Error', 'Samples used for comparison should be selected before running!');
			}
		});
		$hbox->add($button_2);
		$frame->add($hbox) ;
		$table->attach_defaults($frame,2,4,10,11);

	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}

###########################
sub CellRenderer_window{
	my ($matrix_pointer)=@_;
	my @matrix=@$matrix_pointer;
	my $row_num=@matrix;
	my $first_line_pointer=shift @matrix;
	my @first_line=@$first_line_pointer;
	my $col_num=@first_line;
	
	my @model_attr;
	for(my $i=0; $i<$col_num; $i++){
		$model_attr[$i]='Glib::String';
	}
	my $model = Gtk2::ListStore->new (@model_attr);
	#display content
	for(my $m=0; $m<$row_num; $m++){
		if($col_num==1){	$model->set ($model->append, 0 => $matrix[$m][0]);	}
		elsif($col_num==2){	$model->set ($model->append, 0 => $matrix[$m][0], 1=>$matrix[$m][1]);	}
		else{	$model->set ($model->append, 0 => $matrix[$m][0], 1=>$matrix[$m][1], 2=>$matrix[$m][2]);	}
	}
	
	#add column
	my $view = Gtk2::TreeView->new ($model);
	for(my $i=0; $i<$col_num; $i++){
		my $renderer = Gtk2::CellRendererText->new;
		my $column = Gtk2::TreeViewColumn->new_with_attributes ($first_line[$i], $renderer, text => $i);
		$view->append_column ($column);
	}
	
	#scroller window
	my $scroller = Gtk2::ScrolledWindow->new;
	$scroller->set_policy (qw(automatic automatic));
	$scroller->add ($view);
	
	#checkbutton
	my $check = Gtk2::CheckButton->new ('resizable columns');
	$check->set_active (FALSE);
	$check->signal_connect (toggled => sub {
		map { $_->set_resizable ($check->get_active); } $view->get_columns;
	});

	my $box = Gtk2::VBox->new;
	$box->add ($scroller);
	$box->pack_start ($check, FALSE, FALSE, 0);

	return($box);
}

