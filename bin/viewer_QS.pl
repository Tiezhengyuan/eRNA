#!/usr/bin/perl -w
use strict;
use Cwd;
use Glib qw(TRUE FALSE);
use Gtk2 -init;
use threads;
use threads::shared;
use List::Util;
use Statistics::R;
use Gtk2::Pango;

############################################################
#get the directory of perl scripts involved in Pscore
our $perl_dir=Cwd::getcwd();

#get subroutines 
require $perl_dir."/eRNA_subroutines.pm";

#setup shared hash to control thread
my %shash;
share(%shash); #will work for first level keys
$shash{'go'} = 0;
$shash{'data'} = '';
$shash{'work'} = '';
$shash{'die'} = 0;

#get variables
my $var_file=$ARGV[0];
our $variables_pointer=E_RNA::process_info($var_file);
our %variables=%$variables_pointer;
our @sample_names=split( ",", E_RNA::read_log($var_file, 'sample_names') );

our (@candidates, @selected_names);
foreach my $sample_name( @sample_names ){
	push(@candidates, $sample_name) if -f $variables{result_dir}.'/'.$sample_name.'.R1_QC' or -f $variables{result_dir}.'/'.$sample_name.'.R2_QC';
}
@candidates=List::MoreUtils::uniq(@candidates);
#print "@candidates\n";

######################################
#GUI of monitor
my $window = Gtk2::Window->new;
$window->set_size_request('1000', '700');
$window->signal_connect (delete_event => sub {Gtk2->main_quit});
$window->set_title('QS Viewer for viewing sequencing quality scores');
$window->set_position('center');
	my $table=Gtk2::Table->new(10,6,TRUE);
		#setup input
		my $frame=Gtk2::Frame->new('Select samples');
		$frame->set_border_width(10);
			my $sub_table=items_selection();
		$frame->add($sub_table);
	$table->attach_defaults($frame, 0,4,0,3);
		my $tbutton = Gtk2::Button->new_with_label('Apply and run');
		my $sconnect;
		my $lconnect = $tbutton->signal_connect( clicked => \&launch);
	$table->attach_defaults($tbutton, 4,6,1,2);
	
		#picture viewer
		my $sw=Gtk2::ScrolledWindow->new();
		$sw->set_policy('automatic', 'automatic');
		$sw->set_border_width(10);
	$table->attach_defaults($sw, 0,6,3,9);

		#running mode
		my $hbox = Gtk2::HBox->new( FALSE, 5 );
		$hbox->set_border_width(20);
			my $pbar = Gtk2::ProgressBar->new();
			$pbar->set_pulse_step(.1);
			$pbar->hide; #needs to be called after show_all
		$hbox->add($pbar);
			my $button_view = Gtk2::Button->new_with_label('View graphs');
			$button_view->signal_connect(clicked=>sub{
					my $view = create_view();  #subroutine
					$sw->add($view);
					$sw->show_all;
			});
		$hbox->add($button_view);
			my $ebutton = Gtk2::Button->new_from_stock('gtk-quit');
			$ebutton->signal_connect( clicked =>\&exit );
		$hbox->add($ebutton);
	$table->attach_defaults($hbox, 0,6,9,10);
	$table->show_all;
$window->add ($table);
$window->show_all;
Gtk2->main;

#################################################################

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
		#create 1 sleeping thread
		my $thread_1 = threads->new(\&work);
			
		$shash{'go'} = 1;
		$pbar->show;
		$tbutton->set_label('Stop');
		$tbutton->signal_handler_block($lconnect);
		$sconnect = $tbutton->signal_connect( clicked => sub{ 	
			$shash{'go'} = 0;
			$tbutton->set_label('Run');
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
				$pbar->set_text('OK! It is done!');
				$tbutton->set_label('Run');
				return FALSE;
			}
		});

}
################################################## #######
sub work{
	$|++;
	while(1){
		return if $shash{'die'} == 1;

		if ( $shash{'go'} == 1 ){#2
			foreach my $sample_name(@selected_names){
				my $plot_file=QC_plot($variables{result_dir}, $sample_name);
				last if $shash{'go'} == 0;
				return if $shash{'die'} == 1;
			} #3
			$shash{'go'} = 0;   #turn off 
		}#2
		else{ sleep 1; }
	}
}
#################################
sub create_view {
	my $view = Gtk2::TextView->new;
	$view->can_focus(1);
	$view->set_editable(0);
	$view->set_left_margin(10);
	$view->set_right_margin(10);
	$view->set_wrap_mode('GTK_WRAP_CHAR');
	my $buffer = Gtk2::TextBuffer->new();
	my $iter=$buffer->get_iter_at_offset(0);
	$buffer->create_tag("bold", weight=>PANGO_WEIGHT_BOLD);
	$buffer->insert_with_tags_by_name($iter, "Selected samples: @selected_names\n\n\n", "bold");

	foreach my $sample_name(@selected_names){
		my $plot_file=$variables{result_dir}.'/QC_'.$sample_name.'.jpg';
		if(-f $plot_file){
			$buffer->insert($iter, "QC plot of $sample_name ($plot_file):\n");
			my $pix_buffer=Gtk2::Gdk::Pixbuf->new_from_file($plot_file);
			$buffer->insert_pixbuf($iter, $pix_buffer);
			$buffer->insert($iter, "\n");
		}
	}
	$view->set_buffer($buffer);
	$view->show_all;

	return $view;
}

######################################
sub QC_plot{
	my ($dir, $sample_name)=@_;
	$dir .='/' unless $dir=~/\/$/;
	my $plot_file=$dir.'QC_'.$sample_name.'.jpg';
	
	# Create a communication bridge with R and start R
	my $R=Statistics::R->new();
	$R->set('out_file', $plot_file);

	my $R1_file=$dir.$sample_name.'.R1_QC';
	my $R2_file=$dir.$sample_name.'.R2_QC';
	$R->set('R1_file', $R1_file);
	$R->set('R2_file', $R2_file);
	$R->set('sample_name', $sample_name);
	if(-f $R1_file and -f $R2_file){
		$R->run(q'df_1<-read.delim(R1_file, header=F, sep="\t")',
				qq'jpeg(out_file, units="mm", height=120, width=5*ncol(df_1),res=200, pointsize=8)',
				q'par(mfrow=c(2,1), mar=c(4,4,2,2))',
				q'names(df_1)<-seq(1,ncol(df_1))',
				q'boxplot(df_1, cex=0.5, xlab=paste("Sequencing cycles, ", sample_name), ylab="Quality score", main="Quality score distribution per sequencing circyle (R1)")',
				q'df_2<-read.delim(R2_file, header=F, sep="\t")',
				q'names(df_2)<-seq(1,ncol(df_2))',
				q'boxplot(df_2, cex=0.5, xlab=paste("Sequencing cycles, ", sample_name), ylab="Quality score", main="Quality score distribution per sequencing circyle (R2)")',
			);
	}
	elsif(-f $R1_file){
				$R->run(q'df_1<-read.delim(R1_file, header=F, sep="\t")',
				qq'jpeg(out_file, units="mm", height=60, width=5*ncol(df_1),res=200, pointsize=8)',
				q'names(df_1)<-seq(1,ncol(df_1))',
				q'boxplot(df_1, cex=0.5, xlab=paste("Sequencing cycles, ", sample_name), ylab="Quality score", main="Quality score distribution per sequencing circyle (R1)")',
			);
	}
	$R->run(q'dev.off()');
	
	return($plot_file);
}

##########################################################
#used for items' selection
sub items_selection{

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
		foreach (@candidates) {
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
		foreach (@candidates) {
			my $iter = $tree_store_2->append(undef); #the iter is a pointer in the treestore. We use to add data.
			$tree_store_2->set ($iter,0 => $_);
		}
		@selected_names=@candidates;
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
		foreach (@candidates) {
			my $iter = $tree_store_1->append(undef); #the iter is a pointer in the treestore. We use to add data.
			$tree_store_1->set ($iter,0 => $_);
		}
		undef @selected_names;
	} );
	$table->attach_defaults($button_4, 3,4,4,5);

	return ($table);
}

