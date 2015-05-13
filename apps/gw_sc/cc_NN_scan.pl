#!/usr/bin/perl

require "../clue/os_profile.pl";

system("cd $GlobalCodeDir/sd_metropolis/");

$AppDir     = "$GlobalCodeDir/sd_metropolis/apps/gw_sc"; 
$DataDir    = "$GlobalDataDir/sd_metropolis/data/gw_sc";
$TmpDir     = "$GlobalCodeDir/sd_metropolis/tmp";

$executable = "$GlobalCodeDir/sd_metropolis/bin/gw_sc".$ExeExt;

#system("make gw_sc");
system("$OutputCleanupCmd");

$lambda = 3.0;

$cc_min = 15.0; $cc_max = 30.0; $cc_delta = 1.0;
$NN_min = 2.0;  $NN_max = 4.0;  $NN_delta = 0.2;

$nmc  = 50000000;

$pplus_tuning_interval = $nmc/10;
$mc_reporting_interval = $nmc/100;

$maxn = 20000;

$suffix = sprintf("l%2.4f_nmc%i", $lambda, $nmc);

$mc_stat_file = "$DataDir/mc_stat_$suffix.dat";
#system("rm -f -v $mc_stat_file");

$observables_file = "$DataDir\\G_$suffix.dat";
#system("rm -f -v $observables_file");

for($cc=$cc_min; $cc<=$cc_max; $cc += $cc_delta)
{
 for($NN=$NN_min; $NN<=$NN_max; $NN += $NN_delta)
 {
  $suffix1 = sprintf("%s_cc%2.4f_NN%2.4f", $cc, $NN);
  $job_id  = $suffix1;
  
  $cmd = $executable;
  $cmd = $cmd." --lambda                   $lambda";
  $cmd = $cmd." --prod-mc-steps            $nmc";
  $cmd = $cmd." --p-plus-tuning-interval   $pplus_tuning_interval";
  $cmd = $cmd." --max-recursion-depth      $maxn";
  $cmd = $cmd." --mc-stat-file             $mc_stat_file";
  $cmd = $cmd." --observables-file         $observables_file";
  $cmd = $cmd." --logs-noise-level         1";
  $cmd = $cmd." --mc-reporting-interval    $mc_reporting_interval";
  $cmd = $cmd." --no-stack-check ";
  $cmd = $cmd." --exit-upon-overflow ";
  
  $cmd = $cmd." --no-param-auto-tuning ";
  $cmd = $cmd." --cc                       $cc";
  $cmd = $cmd." --NN                       $NN";
  
  if($myOS eq "Linux")
  {
   $cmd = $cmd." --no-ansi-colors ";
   $cmd = $cmd." --print-errors-to-stderr ";
  };
 
  run_command($job_id, $cmd, "short", "$TmpDir");
 };
};
