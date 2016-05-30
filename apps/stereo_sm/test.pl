#!/usr/bin/perl

require "../clue/os_profile.pl";

system("cd $GlobalCodeDir/sd_metropolis/");

$AppDir     = "$GlobalCodeDir/sd_metropolis/apps/smlm"; 
$DataDir    = "$GlobalDataDir/sd_metropolis/data/smlm";
$TmpDir     = "$GlobalCodeDir/sd_metropolis/tmp";

$executable = "$GlobalCodeDir/sd_metropolis/bin/smlm$ExeExt";

system("make smlm");

system("$OutputCleanupCmd");

$DIM        =  1;
$LS         =  2;
$lambda     =  0.6;
$alpha      =  0.2;
$cc         =  3.0;
$NN         =  1.5;
$nmc        =  1500000000;

$meff_sq    =  0.0;
$maxn       =  20000;

$mc_reporting_interval = int($nmc/100);
$pplus_tuning_interval = int($nmc/20);

for($LS=2; $LS<=2000; $LS *= 10)
{
 $suffix0 = sprintf("nmc%i_a%2.2lf_l%2.2lf_s%i", $nmc, $alpha, $lambda, $LS);

 $mc_stat_file = "$DataDir/mc_stat_$suffix0.dat";
 system("rm -f -v $mc_stat_file");

 $action_stat_file = "$DataDir/action_stat_$suffix0.dat";
 system("rm -f -v $action_stat_file");

 $observables_file = "$DataDir/G2_$suffix0.dat";
 system("rm -f -v $observables_file");

for($itrial=0; $itrial<10; $itrial ++)
{
 $suffix1 = sprintf("%s_trial%i", $suffix0, $itrial);
 $job_id  = $suffix1;
 
 $mean_link = 1.0; #($lambda<4.0? 1.0 - $lambda/8 : 2/$lambda);
 
 $stack_stat_file = "$DataDir/stack_stat_$suffix1.dat";
 
 $cmd = $executable;
 $cmd = $cmd." --DIM                      $DIM ";
 $cmd = $cmd." --LT                       $LS ";
 $cmd = $cmd." --LS                       $LS ";
 $cmd = $cmd." --lambda                   $lambda ";
 $cmd = $cmd." --meff-sq                  $meff_sq ";
 $cmd = $cmd." --prod-mc-steps            $nmc ";
 $cmd = $cmd." --max-recursion-depth      $maxn ";
 $cmd = $cmd." --logs-noise-level         1 ";
 $cmd = $cmd." --mc-reporting-interval    $mc_reporting_interval ";
 $cmd = $cmd." --p-plus-tuning-interval   $pplus_tuning_interval ";
 
 $cmd = $cmd." --cc $cc ";
 $cmd = $cmd." --NN $NN  ";
 $cmd = $cmd." --alpha $alpha ";
 $cmd = $cmd." --max-correlator-order 10 ";
 
 
 #Extremely important: mean link value!!!
 $cmd = $cmd." --mean-link                 $mean_link ";

 $cmd = $cmd." --observables-file          $observables_file "; 
 $cmd = $cmd." --mc-stat-file              $mc_stat_file ";
 #$cmd = $cmd." --action-stat-file         $action_stat_file ";
 #$cmd = $cmd." --stack-stat-file          $stack_stat_file ";
 
 $cmd = $cmd." --exit-upon-overflow ";
 $cmd = $cmd." --no-stack-check ";
 #$cmd = $cmd." --no-ansi-colors ";
 
 run_command($job_id, $cmd, "short", "$TmpDir"); 
};
};

system("$QueueWatchCmd");

