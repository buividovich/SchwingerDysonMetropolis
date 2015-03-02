#!/usr/bin/perl

#TODO: Decode command line arguments
$app = "gw_sc";

require "../clue/os_profile.pl";
require "./apps/app_globals.pl";

$lambda   = 5.0;
$meff_sq  = 0.0;
$LT       = 2;
$nmc      = 20000000;
$mc_reporting_interval = int($nmc/100);

$suffix = sprintf("l%2.4f_nmc%i", $lambda, $nmc);

$global_observables_file = "$DataDir/G2_$suffix.dat";
$global_mc_stat_file     = "$DataDir/mc_stat_$suffix.dat";

system("rm -f -v $global_observables_file");
system("rm -f -v $global_mc_stat_file");

@all_mc_stat_files = ();
@all_observa_files = ();

require "./apps/get_mc_stat.pl";

my $cc0; my $NN0; my $v0;
if(scalar(@ARGV)==2)
{
 $cc0 = $ARGV[0];
 $NN0 = $ARGV[1];
 printf("I got cc = %2.4E, NN = %2.4E from command line\n", $cc0, $NN0); 
}
else
{
 ($cc0, $NN0, $v0) = tune_cc_NN();
 printf("Search has terminated with cc=%2.4lf, NN = %2.4lf, <nA> = %2.4lf", $cc0, $NN0, $v0);
};

lambda_descent($cc0, $NN0);

#Now join all the mc stats
#system("cat @all_mc_stat_files > $global_mc_stat_file");
#system("rm -f -v @all_mc_stat_files");
