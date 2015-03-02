#!/usr/bin/perl

require "../clue/os_profile.pl";

system("cd $GlobalCodeDir/sd_metropolis/");

$AppDir     = "$GlobalCodeDir/sd_metropolis/apps/ssm"; 
$DataDir    = "$GlobalDataDir/sd_metropolis/data/ssm";
$TmpDir     = "$GlobalCodeDir/sd_metropolis/tmp";

$executable = "$GlobalCodeDir/sd_metropolis/bin/ssm";

system("make ssm");
system("$OutputCleanupCmd");

$lambda_min = 0.005;
$lambda_max = 0.1;

if($myOS eq "Linux")
{
 $dlambda    = 0.005;
 $nmc  = 2000000000;
 $maxn = 20000;
 @masses = (0.5, 1.0, 2.0); 
}
else
{
 $dlambda    = 0.005;
 $nmc  = 10000000;
 $maxn = 20000;
 @masses = (0.5, 1.0, 2.0);
};

foreach $meff_sq (@masses)
{
 $mc_reporting_interval = int($nmc/100);

 $suffix0 = sprintf("m%2.4f_nmc%i", $meff_sq, $nmc);

 $mc_stat_file = "$DataDir/mc_stat_$suffix0.dat";
 system("rm -f -v $mc_stat_file");

 $action_stat_file = "$DataDir/action_stat_$suffix0.dat";
 system("rm -f -v $action_stat_file");

 $observables_file = "$DataDir/G_$suffix0.dat";
 system("rm -f -v $observables_file");

 for($lambda=$lambda_min; $lambda<=$lambda_max; $lambda += $dlambda)
 {
  $suffix1 = sprintf("%s_l%2.4f", $suffix0, $lambda);
  $job_id  = $suffix1;
 
  $cmd = $executable;
  $cmd = $cmd." --lambda                   $lambda";
  $cmd = $cmd." --meff-sq                  $meff_sq";
  $cmd = $cmd." --prod-mc-steps            $nmc";
  $cmd = $cmd." --max-recursion-depth      $maxn";
  $cmd = $cmd." --mc-stat-file             $mc_stat_file";
  $cmd = $cmd." --action-stat-file         $action_stat_file";
  $cmd = $cmd." --observables-file         $observables_file";
  $cmd = $cmd." --logs-noise-level         1";
  $cmd = $cmd." --mc-reporting-interval    $mc_reporting_interval";
  
  $cmd = $cmd." --no-stack-check ";
  
  if($myOS eq "Linux")
  {
   $cmd = $cmd." --no-stack-check ";
   $cmd = $cmd." --no-ansi-colors ";
   $cmd = $cmd." --print-errors-to-stderr ";
  };
 
  run_command($job_id, $cmd, "short", "$TmpDir");
 };
};

system("$QueueWatchCmd");

