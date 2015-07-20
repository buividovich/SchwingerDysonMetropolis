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

$lambda_min =  60.0;
$lambda_max =  60.0;
$dlambda    =  2.5;
$nmc        =  30000000;
$maxn       =  20000;
$meff_sq    =  -2.0*$DIM;

$mc_reporting_interval = int($nmc/100);
$pplus_tuning_interval = int($nmc/20);

#$suffix0 = sprintf("gw_d%i_s%i_m%2.2f_nmc%i", $DIM, $LS, $meff_sq, $nmc);
$suffix0 = sprintf("gw_m%2.2f_nmc%i", $meff_sq, $nmc);

$mc_stat_file = "$DataDir/mc_stat_$suffix0.dat";
#system("rm -f -v $mc_stat_file");

$action_stat_file = "$DataDir/action_stat_$suffix0.dat";
#system("rm -f -v $action_stat_file");

$observables_file = "$DataDir/G_$suffix0.dat";
#system("rm -f -v $observables_file");

for($lambda=$lambda_min; $lambda<=$lambda_max; $lambda += $dlambda)
{
 $suffix1 = sprintf("%s_l%2.4f", $suffix0, $lambda);
 $job_id  = $suffix1;
 
 $mean_link = ($lambda<4? 1.0 - $lambda/8 : 2/$lambda);
 
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
 
 $cmd = $cmd." --cc 8.0 ";
 $cmd = $cmd." --NN 2.5  ";
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

system("$QueueWatchCmd");

