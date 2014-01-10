#!/usr/bin/perl
use strict;
use warnings;
use Config::General;
use File::Slurp;
#use Data::Dumper;

sub fileCheck;
sub runInParallel;

my $cfgFile="config.txt";
my $cfg=new Config::General($cfgFile);
my %cfgHash=$cfg->getall;

#print Data::Dumper->Dumper(\%cfgHash);

##########################################################################
# READ IN THE LIST OF SAMPLE NAMES AND CORRESPONDING BAM FILES
##########################################################################

my $cmd;
my $bamList=$cfgHash{bamfilelist};
fileCheck($bamList,"File containing list of bam files and corresponding sample names is required!\nCheck $cfgFile!");

open BL, "<$bamList";
my %bam2sample;
while (my $line=<BL>){
    chomp $line;
    my @tmp=split/\t/,$line;
    fileCheck($tmp[1],"Check $bamList!");
    $bam2sample{$tmp[1]}=$tmp[0];
}

##########################################################################
# SPLIT THE BAM FILES INTO CHROMOSOME LEVEL SMALLER BAM FILES
##########################################################################

open TMP, ">VSP.parallel.1";
for my $bamFile (keys %bam2sample) {
    $cmd=$cfgHash{bamtoolsexe}." split -in ".$bamFile." -reference";
    print TMP "$cmd\n";
}
close TMP;
#runInParallel("VSP.parallel.1");

##########################################################################
# GET LIST OF CHROMOSOMES
##########################################################################

my $done=0;
my @chrlist;
for my $bamFile (keys %bam2sample) {
    last if $done==1;
    $done++;
    my $bamFileName_wo_ext=$bamFile;
    $bamFileName_wo_ext=~s/\.\w+$//;
    my @chrbams = <${bamFileName_wo_ext}.REF_*.bam>;    
    for my $chr (@chrbams) {
        my @tmp=$chr=~/${bamFileName_wo_ext}.REF_(.*).bam/;
        push @chrlist,$tmp[0] unless $tmp[0] eq "unmapped";
    }
}

##########################################################################
# CALL VARIANTS AT CHROMOSOME LEVEL
##########################################################################

my $minCov=$cfgHash{Minimum_read_depth_at_a_position_to_make_a_call};
my $minReads=$cfgHash{Minimum_supporting_reads_at_a_position_to_call_variants};
my $minBaseQual=$cfgHash{Minimum_base_quality_at_a_position_to_count_a_read};
open TMP, ">VSP.parallel.2";
$done=0;
my %realSampleName;
for my $chr (@chrlist) {
    $cmd="$cfgHash{samtoolsexe} mpileup -f $cfgHash{reference}";
    my $ctr=0;
    for my $bamFile (keys %bam2sample) {
        $ctr++;
        $realSampleName{"Sample${ctr}"}=$bam2sample{$bamFile};
        my $chrBamFile=$bamFile;
        $chrBamFile=~s/\.\w+$//;
        $chrBamFile.=".REF_${chr}.bam";
        $cmd.=" $chrBamFile";
    }
    $cmd.=" | ".$cfgHash{javaexe}." -jar ".$cfgHash{varscanexe}." mpileup2snp - --min-coverage ".$minCov." --min-reads2 ".$minReads." --min-avg-qual ".$minBaseQual." --output-vcf 1 > ${chr}.vcf";
    print TMP "$cmd\n";
    $done++;
}
close TMP;
#runInParallel("VSP.parallel.2");

##########################################################################
# CONCATENATE ALL CHROMOSOME LEVEL VARIANT CALLS
##########################################################################

$cmd=$cfgHash{vcftoolsbin}."/vcf-concat";
for my $chr (@chrlist) {
    $cmd.=" ${chr}.vcf";
}
$cmd.=" > SNPs.vcf";
#system($cmd);

##########################################################################
# RENAME SAMPLE NAMES FROM Sample1, Sample2, etc. TO REAL NAMES FROM
# BAMFILELIST FILE
##########################################################################

$cmd="cat SNPs.vcf|grep -m100 ^# > header.txt";
system($cmd);
fileCheck("header.txt","Something went wrong, SNPs.vcf does not have a header!");
my $header=File::Slurp::read_file("header.txt");
my $sampleNumber=0;
for my $sampleName (keys %realSampleName) {
    my $newname=$realSampleName{$sampleName};
    $header=~s/${sampleName}/${newname}/g;
}
File::Slurp::write_file("header.txt",$header);
$cmd="cat SNPs.vcf|grep -v ^# >> header.txt";
system($cmd);
$cmd="mv header.txt SNPs.vcf";
system($cmd);
$cmd=$cfgHash{bgzipexe}." SNPs.vcf";
system($cmd);
$cmd=$cfgHash{tabixexe}." -p vcf SNPs.vcf.gz";
system($cmd);

##########################################################################
# SUBROUTINES
##########################################################################


sub fileCheck{
    my ($filename,$msg)=@_;
    $msg=$msg?$msg:"";
    if ( ! -f $filename ) {
        print "File ".$filename." does not exist!\n$msg\n";
        print "Exiting!!!\n";
        exit;
    }
}


sub runInParallel{
    my ($filename)=@_;
    $cmd=$cfgHash{parallelexe}." -j ".$cfgHash{ncpus}." < ".$filename;
    system($cmd);
}