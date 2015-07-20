#!/usr/bin/perl

require "../clue/os_profile.pl";

system("cd $GlobalCodeDir/sd_metropolis/");

$AppDir     = "$GlobalCodeDir/sd_metropolis/apps/lmsctest"; 
$DataDir    = "$GlobalDataDir/sd_metropolis/data/lmsctest";
$TmpDir     = "$GlobalCodeDir/sd_metropolis/tmp";

$executable = "$GlobalCodeDir/sd_metropolis/bin/lmsctest$ExeExt";

system("make lmsctest");

system("$OutputCleanupCmd");

$nmc        =  150000000;
$maxn       =  20000;

$mc_reporting_interval = int($nmc/100);
$pplus_tuning_interval = int($nmc/20);
 
 $cmd = $executable;
 $cmd = $cmd." --model                    1 ";
 
 $cmd = $cmd." --prod-mc-steps            $nmc ";
 $cmd = $cmd." --max-recursion-depth      $maxn ";
 $cmd = $cmd." --logs-noise-level         1 ";
 $cmd = $cmd." --mc-reporting-interval    $mc_reporting_interval ";
 $cmd = $cmd." --p-plus-tuning-interval   $pplus_tuning_interval ";
 $cmd = $cmd." --max-correlator-order     10 ";
 
 $cmd = $cmd." --cc 8.0 ";
 $cmd = $cmd." --NN 2.5  ";
 
 $cmd = $cmd." --epsilon 1.0 ";
 
 $cmd = $cmd." --exit-upon-overflow ";
 #$cmd = $cmd." --no-stack-check ";
 #$cmd = $cmd." --no-ansi-colors ";
 
 run_command("lmsctest", $cmd, "short", "$TmpDir"); 

