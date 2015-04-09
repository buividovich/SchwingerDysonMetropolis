#!/usr/bin/perl

require "../clue/os_profile.pl";

system("cd $GlobalCodeDir/sd_metropolis/");

$AppDir     = "$GlobalCodeDir/sd_metropolis/apps/ssm"; 
$DataDir    = "$GlobalDataDir/sd_metropolis/data/ssm";
$TmpDir     = "$GlobalCodeDir/sd_metropolis/tmp";

$executable = "$GlobalCodeDir/sd_metropolis/bin/ssm";

system("make ssm");
system("$OutputCleanupCmd");

$lambda = 5.0;

$cc_min = 14.0; $cc_max = 18.0; $cc_delta = 0.4;
$NN_min = 0.2;  $NN_max = 2.4;  $NN_delta = 0.2;

$nmc  = 5000000;
$maxn = 20000;

$suffix = sprintf("l%2.4f_nmc%i", $lambda, $nmc);

$mc_stat_file = "$DataDir/mc_stat_$suffix.dat";
system("rm -f -v $mc_stat_file");

$observables_file = "$DataDir\\G_$suffix.dat";
system("rm -f -v $observables_file");

for($cc=$cc_min; $cc<=$cc_max; $cc += $cc_delta)
{
 for($NN=$NN_min; $NN<=$NN_max; $NN += $NN_delta)
 {
  $suffix1 = sprintf("%s_cc%2.4f_NN%2.4f", $cc, $NN);
  $job_id  = $suffix1;
  
  $cmd = $executable;
  $cmd = $cmd." --lambda                   $lambda";
  $cmd = $cmd." --prod-mc-steps            $nmc";
  $cmd = $cmd." --max-recursion-depth      $maxn";
  $cmd = $cmd." --mc-stat-file             $mc_stat_file";
  $cmd = $cmd." --observables-file         $observables_file";
  $cmd = $cmd." --logs-noise-level         1";
  $cmd = $cmd." --mc-reporting-interval    50000";
  $cmd = $cmd." --no-stack-check ";
  
  $cmd = $cmd." --no-param-auto-tuning ";
  $cmd = $cmd." --cc                       $cc";
  $cmd = $cmd." --NN                       $NN";
  
  if($myOS eq "Linux")
  {
   $cmd = $cmd." --no-stack-check ";
   $cmd = $cmd." --no-ansi-colors ";
   $cmd = $cmd." --print-errors-to-stderr ";
  };
 
  run_command($job_id, $cmd, "short", "$TmpDir");
 };
};
