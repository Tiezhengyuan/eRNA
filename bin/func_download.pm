#! /usr/bin/perl -w
use strict;
use warnings;


####################################################
#
package sub_download;

##############################
sub check_gcc{
	my $gcc=`gcc --version 2>&1`;
	if($gcc !~ /(GCC)/i){
		die "\nError:\n\tno gcc compiler installed. Please install a gcc compiler\n";
	}else{
		if($gcc =~ /^gcc\s*\S*\s*(\d+\S+)\s*/){
			print STDERR "gcc version: $1                                      already installed, nothing to do ...\n";
		}
	} 
	#
}

########################
#download tool
sub check_wget{
	#capture STDERR as well as STDOUT
	my $wget=`wget 2>&1`;
	if($wget =~ /URL/i){
		printf ("###The download tool: wget\n");
	}else{
		die "No download tool is detected. Please install wget on your machine\n";
	}
}

########################
#download tool
sub check_gzip{
	#capture STDERR as well as STDOUT
	my $sf=`gzip -V 2>&1`;
	if($sf =~ /gzip/i){
		printf ("###The download tool: %s\n", $sf);
	}else{
		die "Error: Please install wget on your machine\n";
	}
}
##########################
#The mapper: tophat
sub get_tophat{
	my($dir)=@_;
	#dir: download
	my $dir_download=$dir.'/download/';
	mkdir($dir_download, 0755) unless -d $dir_download;
	
	#download
	my $ver='2.1.1';
	my $url='https://ccb.jhu.edu/software/tophat/downloads/tophat-'.$ver.'.Linux_x86_64.tar.gz';
	system("wget -c $url -P $dir_download "); 
	
	#uncompress
	my $local_zip= $dir_download.'tophat-'.$ver.'.Linux_x86_64.tar.gz';
	system("tar xvzf $local_zip -C $dir_download");
	
	#cp exe files 
	my $dir_exe=$dir.'/tophat/';
	mkdir($dir_exe, 0755) unless -d $dir_exe;
	my $tmp_dir=$dir_download.'tophat-'.$ver.'.Linux_x86_64';
	system("cp -r  $tmp_dir/* $dir_exe");
	#
}

########################################
#The aligner: bowtie1
sub get_bowtie1{
	my($dir)=@_;
	#dir: download
	my $dir_download=$dir.'/download/';
	mkdir($dir_download, 0755) unless -d $dir_download;
	
	#download
	my $ver_bowtie1='1.2';
	my $url_bowtie1='https://sourceforge.net/projects/bowtie-bio/files/bowtie/'.$ver_bowtie1.'.0/bowtie-'.$ver_bowtie1.'-linux-x86_64.zip';
	system("wget -c $url_bowtie1 -P $dir_download "); 
	
	#uncompress
	my $local_bowtie1= $dir_download.'bowtie-'.$ver_bowtie1.'-linux-x86_64.zip';
	system("unzip -u $local_bowtie1 -d $dir_download");
	
	#cp bowtie exe files 
	my $dir_bowtie1=$dir.'/bowtie1/';
	mkdir($dir_bowtie1, 0755) unless -d $dir_bowtie1;
	my $tmp=$dir_download.'bowtie-'.$ver_bowtie1;
	system("cp -r  $tmp/* $dir_bowtie1");
	#
}

########################################
#The aligner: bowtie2
sub get_bowtie2{
	my($dir)=@_;
	#dir: download
	my $dir_download=$dir.'/download/';
	mkdir($dir_download, 0755) unless -d $dir_download;
	
	#download
	my $ver_bowtie2='2.3.2';
	my $url_bowtie2='https://sourceforge.net/projects/bowtie-bio/files/bowtie2/'.$ver_bowtie2.'/bowtie2-'.$ver_bowtie2.'-linux-x86_64.zip';
	system("wget -c $url_bowtie2 -P $dir_download "); 
	
	#uncompress
	my $local_bowtie2= $dir_download.'bowtie2-'.$ver_bowtie2.'-linux-x86_64.zip';
	system("unzip -u $local_bowtie2 -d $dir_download");
	
	#cp bowtie exe files 
	my $dir_bowtie2=$dir.'/bowtie2/';
	mkdir($dir_bowtie2, 0755) unless -d $dir_bowtie2;
	my $tmp=$dir_download.'bowtie2-'.$ver_bowtie2;
	system("cp -r  $tmp/* $dir_bowtie2");
}
#####################################
##########################
#The asemblerer: cufflinks
sub get_cufflinks{
	my($dir)=@_;
	#dir: download
	my $dir_download=$dir.'/download/';
	mkdir($dir_download, 0755) unless -d $dir_download;
	
	#download
	my $ver='2.1.1';
	my $url='http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-'.$ver.'.Linux_x86_64.tar.gz';
	system("wget -c $url -P $dir_download "); 
	
	#uncompress
	my $local_zip= $dir_download.'cufflinks-'.$ver.'.Linux_x86_64.tar.gz';
	system("tar xvzf $local_zip -C $dir_download");
	
	#cp exe files 
	my $dir_exe=$dir.'/cufflinks/';
	mkdir($dir_exe, 0755) unless -d $dir_exe;
	my $tmp_dir=$dir_download.'cufflinks-'.$ver.'.Linux_x86_64';
	system("cp -r  $tmp_dir/* $dir_exe");
	#
}

##########################
#SAMtools
sub get_samtools{
	my($dir)=@_;
	#dir: download
	my $dir_download=$dir.'/download/';
	mkdir($dir_download, 0755) unless -d $dir_download;
	
	#download
	my $ver='1.4.1';
	my $url='https://sourceforge.net/projects/samtools/files/samtools/'.$ver.'/samtools-'.$ver.'.tar.bz2';
	system("wget -c $url -P $dir_download "); 
	
	#uncompress
	my $local_zip= $dir_download.'samtools-'.$ver.'.tar.bz2';
	system("tar xvjf $local_zip -C $dir_download");
	
	#cp exe files 
	my $dir_exe=$dir.'samtools/';
	mkdir($dir_exe, 0755) unless -d $dir_exe;
	my $tmp_dir=$dir_download.'samtools-'.$ver;
	system("cp -r  $tmp_dir/* $dir_exe");
	
	#
}
########################
1;