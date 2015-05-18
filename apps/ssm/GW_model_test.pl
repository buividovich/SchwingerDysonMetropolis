#!/usr/bin/perl

require "../clue/os_profile.pl";

system("cd $GlobalCodeDir/sd_metropolis/");

$AppDir     = "$GlobalCodeDir/sd_metropolis/apps/ssm"; 
$DataDir    = "$GlobalDataDir/sd_metropolis/data/ssm";
$TmpDir     = "$GlobalCodeDir/sd_metropolis/tmp";

$executable = "$GlobalCodeDir/sd_metropolis/bin/ssm$ExeExt";

system("make ssm");

system("$OutputCleanupCmd");

$lambda_min =  60.0;
$lambda_max =  100.0;
$dlambda    =  5.0;
$nmc        =  500000000;
$maxn       =  20000;
$meff_sq    = -2.0;
$LS         = 16;

$mc_reporting_interval = int($nmc/100);
$pplus_tuning_interval = int($nmc/20);

$suffix0 = sprintf("gwsc_test_m%2.2f_nmc%i", $meff_sq, $nmc);

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
 
 $stack_stat_file = "$DataDir/stack_stat_$suffix1.dat";
 
 $cmd = $executable;
 $cmd = $cmd." --DIM                      2 ";
 $cmd = $cmd." --LT                       $LS ";
 $cmd = $cmd." --LS                       $LS ";
 $cmd = $cmd." --lambda                   $lambda ";
 $cmd = $cmd." --meff-sq                  $meff_sq ";
 $cmd = $cmd." --prod-mc-steps            $nmc ";
 $cmd = $cmd." --max-recursion-depth      $maxn ";
 $cmd = $cmd." --mc-stat-file             $mc_stat_file ";
 $cmd = $cmd." --action-stat-file         $action_stat_file ";
 $cmd = $cmd." --observables-file         $observables_file ";
 $cmd = $cmd." --logs-noise-level         1 ";
 $cmd = $cmd." --mc-reporting-interval    $mc_reporting_interval ";
 $cmd = $cmd." --p-plus-tuning-interval   $pplus_tuning_interval ";
 $cmd = $cmd." --stack-stat-file          $stack_stat_file ";
 
 $cmd = $cmd." --exit-upon-overflow ";
 $cmd = $cmd." --no-stack-check ";
  
 run_command($job_id, $cmd, "short", "$TmpDir");
};

system("$QueueWatchCmd");

