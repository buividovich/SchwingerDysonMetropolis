require "../clue/os_profile.pl";

system("cd $GlobalCodeDir/sd_metropolis/");

$AppDir     = "$GlobalCodeDir/sd_metropolis/apps/hmm"; 
$DataDir    = "$GlobalDataDir/sd_metropolis/data/hmm";
$TmpDir     = "$GlobalDataDir/sd_metropolis/tmp";

#drmemory -show_reachable 
$executable = "$GlobalCodeDir/sd_metropolis/bin/hmm$ExeExt";

system("make hmm");

system("$OutputCleanupCmd");

$lambda_min =  0.0;
$lambda_max =  0.0;
$dlambda    =  0.005;

$nmc  = 20000000;
$maxn = 20000;
$pplus_tuning_interval = int($nmc/10);
$mc_reporting_interval = int($nmc/100);

$mc_stat_file = "$DataDir/mc_stat_nmc$nmc.dat";
system("rm -f -v $mc_stat_file");

$observables_file = sprintf("%s/G_nmc%i.dat", $DataDir, $nmc);
system("rm -f -v $observables_file");

for($lambda=$lambda_min; $lambda<=$lambda_max; $lambda += $dlambda)
{
 $suffix1 = sprintf("l%2.4f_nmc%i", $lambda, $nmc);
 $job_id  = $suffix1;
 
 $stack_stat_file  =  sprintf("%s/stack_stat_%s.dat", $DataDir, $suffix1);
 $ns_history_file  =  sprintf("%s/ns_history_%s.dat", $DataDir, $suffix1);

 $cmd = $executable;
 $cmd = $cmd." --lambda                     $lambda ";
 $cmd = $cmd." --prod-mc-steps              $nmc ";
 $cmd = $cmd." --max-recursion-depth        $maxn ";
 $cmd = $cmd." --mc-stat-file               $mc_stat_file ";
 $cmd = $cmd." --observables-file           $observables_file ";
 $cmd = $cmd." --logs-noise-level           1 ";
 $cmd = $cmd." --mc-reporting-interval      $mc_reporting_interval ";
 $cmd = $cmd." --p-plus-tuning-interval     $pplus_tuning_interval ";
 
 $cmd = $cmd." --cc                         3.46 ";
 $cmd = $cmd." --NN                         2.0 "; #was 1.29
 $cmd = $cmd." --max-correlator-order       7 ";
 
 $cmd = $cmd." --max-genus                  7 ";
 $cmd = $cmd." --genus-f-exponent           1.4 ";
 $cmd = $cmd." --genus-A                    0.3 ";
 
#$cmd = $cmd." --stack-stat-file            $stack_stat_file ";
#$cmd = $cmd." --ns-history-file            $ns_history_file ";
 
 $cmd = $cmd." --exit-upon-overflow ";
 $cmd = $cmd." --no-stack-check ";
#$cmd = $cmd." --no-ansi-colors ";
 
 print $cmd;
 
 run_command($job_id, $cmd, "short", "$TmpDir");
};
