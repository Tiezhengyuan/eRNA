#! /usr/bin/perl -w
use warnings;
use strict;
use Cwd;
use Gtk2 qw/-init -threads-init/;
#use Glib qw/TRUE FALSE/; 
use Gtk2::Pango;
use Gtk2::SimpleList;
use threads; 
use List::Util;
use List::MoreUtils;
use File::Find;
#use constant {TRUE=>1, FALSE=>0};

#get the directory of perl scripts
our $perl_dir=Cwd::getcwd();

#get subroutines
use func_common; #sub_common::
use func_basic; #sub_basic::
use func_rna;  #sub_rna::


#
package sub_gui;

#
###################################################
#used in GUI
sub initiate_eRNA_variables{
	my ($variables_pointer, $perl_dir)=@_;
	my %variables=%$variables_pointer;
	
	#
	if ($perl_dir){
		$variables{dir_home}=$perl_dir;
		$variables{dir_raw_data}=$perl_dir.'/raw_data';
		$variables{dir_result}=$perl_dir.'/result';
	}
	#print $perl_dir;
	#
	#$variables{dir_raw_data}='/home/yuan/data_1/rawdata/a,/home/yuan/data_1/rawdata/b,/home/yuan/data_1/rawdata/c';
	#$variables{dir_result_array}='/home/yuan/eRNA/result,/home/yuan/eRNA/result/a,/home/yuan/eRNA/result/b';
	#initiate variables in variables.txt
	
	#
	$variables{dir_home}=sub_common::format_directory($variables{dir_home});
	$variables{dir_raw_data}=sub_common::format_directory($variables{dir_raw_data});
	$variables{dir_result}=sub_common::format_directory($variables{dir_result});
	$variables{dir_result_array}=$variables{dir_result};
	$variables{dir_log}=sub_common::format_directory($variables{dir_result}.'sample_log');
	$variables{file_total_log}=$variables{dir_log}.'Total.log';
	
	#default files
	$variables{file_var}=$variables{dir_result}.'variables.txt';
	$variables{file_log}=$variables{dir_result}.'output.log';
	$variables{file_time_monitor_log}=$variables{dir_result}.'time_monitor.log';
	$variables{file_system_monitor_log}=$variables{dir_result}.'system_monitor.log';
	$variables{file_references_txt}=$variables{dir_result}.'references.txt';
	$variables{file_targets_df}=$variables{dir_result}.'expression_profilings.txt';
	$variables{file_sample_info}=$variables{dir_result}.'sample_info.csv';
	#rule the directories used for saving results determined by eRNA when large data analysis
	$variables{file_storage_log}=$variables{dir_result}.'storage.log';
	$variables{R_DEseq_file}=$variables{dir_home}.'R_DEseq.R';
	$variables{R_RF_file}=$variables{dir_home}.'R_RF.R';
	
	#directories
	$variables{dir_bowtie1}=$variables{dir_home}.'bowtie1/';
	$variables{dir_bowtie2}=$variables{dir_home}.'bowtie2/';
	$variables{dir_mapper}=$variables{dir_home}.'tophat/';
	$variables{dir_assembler}=$variables{dir_home}.'cufflinks/';
	$variables{dir_mirdeep}=$variables{dir_home}.'mirdeep/';
	$variables{dir_mirspring}=$variables{dir_home}.'mirspring/';
	$variables{dir_statistics}=sub_common::format_directory($variables{dir_result}.'statistics');
	$variables{file_statistics}=$variables{dir_statistics}.'statistics.txt';
	
	#parallel processing
	$variables{'threads_num'}=1;
	#
	#print $variables{file_var};
	sub_common::hash_to_file(\%variables, $variables{file_var}, '=');
		
	return(\%variables);
}

############################################
#Pops up a standard file chooser. # Specify a header to be displayed Specify a type depending on your needs 
#Optionally add a filter to show only certain files-will return a path, if valid

sub show_chooser {
	my($heading, $type, $filter) =@_;
	#$type can be:* 'open', * 'save', * 'select-folder' or * 'create-folder' 
	my $file_chooser =  Gtk2::FileChooserDialog->new ( 
				$heading,
				undef,
				$type,
				'gtk-cancel' => 'cancel',
				'gtk-ok' => 'ok'
		);
	(defined $filter)&&($file_chooser->add_filter($filter));
    
	#if action = 'save' suggest a filename
	($type eq 'save')&&($file_chooser->set_current_name("suggeste_this_file.name"));

	my $filename;
	if ('ok' eq $file_chooser->run){    
		$filename = $file_chooser->get_filename;
	}

	$file_chooser->destroy;

	if (defined $filename){
		if ((-f $filename)&&($type eq 'save')) {
			my $overwrite =show_message_dialog( my $window,
							'question'
							,'Overwrite existing file:'."<b>\n$filename</b>"
							,'yes-no'
				);
			return  if ($overwrite eq 'no');
		}
		return $filename;
	}
	return;
}

#######################################################
#used in subroutine of notebook
sub nb_button{
	my ($label, $window)=@_;
	my $button=Gtk2::Button->new($label);
	$button->signal_connect(clicked => $window);
	return($button);
}
#####################
sub close_button_control{
	my($window)=@_;
	
	my $frame = Gtk2::Frame->new();
	$frame->set_border_width(5);
		my $button=Gtk2::Button->new('Close');
		$button->set_border_width(5);
		$button->signal_connect(pressed=>sub{		$window->destroy;			});
	$frame->add($button);
	
	return($frame);
}
#################
sub default_window{
    my($title, $width,  $height)=@_;
    
	my $window = Gtk2::Window->new('toplevel');
	$window->set_title($title);
	$window->signal_connect('delete_event' => sub { Gtk2->main_quit; });
	$window->set_border_width(5);
	$window->set_position('center');
	$window->set_size_request($width,$height);
	
	return($window);
}

################
sub progress_bar{
	my($interval)=@_;
	
	my $frame=Gtk2::Frame->new();
	$frame->set_label('Progress bar');
	$frame->set_shadow_type('out');
	$frame->set_border_width(5);
		my $pbar = Gtk2::ProgressBar->new();
		$pbar->set_pulse_step($interval);
		$pbar->hide; #needs to be called after show_all
	$frame->add($pbar);
	
	return($frame)
}

############################
sub popup_window{
	my ($title, $str)=@_;
	
	my $window = Gtk2::Window->new();
	$window->set_title($title);
	$window->signal_connect('delete_event' => sub { Gtk2->main_quit; });
	$window->set_border_width(5);
	$window->set_position('center');
	$window->set_size_request('400','100');
		my $vbox=Gtk2::VBox->new(0, 5); #here the first 0 indicates FALSE
			my $label=Gtk2::Label->new($str);
		$vbox->add($label);
			my $button=Gtk2::Button->new_with_label('Close');
			$button->signal_connect(clicked=>sub{		$window->destroy;	 });
		$vbox->add($button);
	$window->add($vbox);
	$window->show_all;
	Gtk2->main();
}

##########
sub entry_control{
	my($entry_label, $entry_text, $type)=@_;
	
	#Bowtie1
	my $frame = Gtk2::Frame->new();
	$frame->set_border_width(5);
	$frame->set_label($entry_label);
		my $subtable=Gtk2::Table->new(1,4, 1);#the last 1 is TRUE
			my $entry= Gtk2::Entry->new();
			$entry->set_text($entry_text);
		$subtable->attach_defaults($entry, 0,4,0,1);
			my $button=Gtk2::Button->new('Select a folder');
			$button->signal_connect('clicked' => sub{ 
				my $text=show_chooser('File Chooser type select-folder','select-folder');
				$text=$entry->get_text().','.$text if $type;
				$entry->set_text($text);
			 });
		$subtable->attach_defaults($button, 4,5,0,1);
	$frame->add($subtable);
	
	return($frame, $entry, $button);
}
################
sub entry_frame{
	my($label, $entry_text)=@_;
	
	my $frame = Gtk2::Frame->new();
	$frame->set_border_width(5);
	$frame->set_label($label);
		my $entry=Gtk2::Entry->new();
		$entry->set_text($entry_text);
	$frame->add($entry);
	
	return($frame, $entry);
}
######################
sub combobox_frame{
	my($label, $entry_text_pointer, $active)=@_;
	my @entry_text=@$entry_text_pointer;
	
	my $frame = Gtk2::Frame->new();
	$frame->set_border_width(5);
	$frame->set_label($label);
		my $entry = Gtk2::ComboBox->new_text;
		foreach (@entry_text){
			$entry->append_text($_);
		}
		$entry->set_active($active);
	$frame->add($entry);
	
	return($frame, $entry);
}
#####################
sub textview_control{
	my($label)=@_;
	
	my $frame = Gtk2::Frame->new();
	$frame->set_border_width(5);
	$frame->set_label($label);
		my $sw = Gtk2::ScrolledWindow->new ();
		$sw->set_shadow_type ('etched-out');
		$sw->set_policy ('automatic', 'automatic');
		$sw->set_border_width(10);
			my $textview = Gtk2::TextView->new();
			$textview->set_wrap_mode('GTK_WRAP_WORD_CHAR');
		$sw->add($textview);
	$frame->add($sw);
	
	return($frame, $textview);
}


################
sub tree_store_str{
	my ($tree_store)=@_;
	
	my @items;
	$tree_store->foreach(sub{
			my($model, $path, $iter)=@_;
			my $value = $model->get ($iter, 0);
			push(@items, $value);
			return(0);
		});
	my $str=join(',', @items);
	return($str);
}

##############################
sub column_treeview{
	my ($liststore, $site_info_items_pointer)=@_;
	my @site_info_items=@$site_info_items_pointer;
	
	my $treeview = Gtk2::TreeView->new ($liststore);
	$treeview->set_rules_hint (1); #TRUE
	$treeview->set_reorderable (1);#TRUE
	my $num=0;
	foreach (@site_info_items){
		# renderer column 0:site_name
		my $renderer = Gtk2::CellRendererText->new;
		$renderer->set (editable => 1);#TRUE
		my $column = Gtk2::TreeViewColumn->new_with_attributes ($_, $renderer, text => $num);
		$renderer->signal_connect (edited => sub {
						my ($cell, $text_path, $new_text, $model) = @_;
						my $path = Gtk2::TreePath->new_from_string ($text_path);
						my $iter = $model->get_iter ($path);
						$model->set ($iter, $num, $new_text);
				}, $liststore);
		$treeview->append_column ($column);
		$num++;
	}
	
	return($treeview);
}


##########################################################
1; 
# make sure the file returns true or require will not succeed!#

