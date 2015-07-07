require "../clue/os_profile.pl";

system("cd $GlobalCodeDir/sd_metropolis/");

$AppDir     = "$GlobalCodeDir/sd_metropolis/apps/hmm"; 
$DataDir    = "$GlobalDataDir/sd_metropolis/data/hmm";
$TmpDir     = "$GlobalDataDir/sd_metropolis/tmp";

#drmemory -show_reachable 
$executable = "$GlobalCodeDir/sd_metropolis/bin/hmm$ExeExt";

system("make hmm");

system("$OutputCleanupCmd");

@lambdas = ( 0.08,  0.07,  0.06,  0.05,  0.04,  0.03,  0.02,  0.01,  0.00);
@ccs     = ( 2.80,  2.60,  2.50,  2.45,  2.40,  3.00,  4.30,  4.40,  4.40);
@nns     = ( 0.80,  0.80,  0.80,  0.80,  0.90,  0.90,  0.90,  2.00,  3.00);
@fes     = ( 5.40,  2.50,  0.80,  0.30,  0.20,  0.05,  0.01, 0.001, 0.001);
@As      = ( 0.20,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.20,  0.20);
@Bs      = (0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005,  1.00,  1.00);

$nmc  = 20000000;
$maxn = 20000;
$pplus_tuning_interval = int($nmc/20);
$mc_reporting_interval = int($nmc/100);

$mc_stat_file = "$DataDir/mc_stat_nmc$nmc.dat";
#system("rm -f -v $mc_stat_file");

$observables_file = sprintf("%s/G_nmc%i.dat", $DataDir, $nmc);
#system("rm -f -v $observables_file");

for($i=0; $i<2000; $i++)
{
 $il = 0;
 $lambda = $lambdas[$il];
 $suffix1 = sprintf("l%2.4f_nmc%i", $lambda, $nmc);
 $job_id  = $suffix1;
 
 $stack_stat_file  =  sprintf("%s/stack_stat_%s.dat", $DataDir, $suffix1);
 $ns_history_file  =  sprintf("%s/ns_history_%s.dat", $DataDir, $suffix1);
 $gmn_hist_file    =  sprintf("%s/gmn_hist_%s.dat",   $DataDir, $suffix1);
 $fcn_data_file    =  sprintf("%s/fcn_data_%s.dat",   $DataDir, $suffix1);

 $cmd = $executable;
 $cmd = $cmd." --lambda                     $lambda ";
 $cmd = $cmd." --prod-mc-steps              $nmc ";
 $cmd = $cmd." --max-recursion-depth        $maxn ";
 $cmd = $cmd." --mc-stat-file               $mc_stat_file ";
 $cmd = $cmd." --observables-file           $observables_file ";
 $cmd = $cmd." --logs-noise-level           1 ";
 $cmd = $cmd." --mc-reporting-interval      $mc_reporting_interval ";
 $cmd = $cmd." --p-plus-tuning-interval     $pplus_tuning_interval ";
 
 $cmd = $cmd." --cc                         $ccs[$il] ";
 $cmd = $cmd." --NN                         $nns[$il] ";
 $cmd = $cmd." --max-correlator-order       7 ";
 
 $cmd = $cmd." --max-genus                  7 ";
 $cmd = $cmd." --genus-f-exponent           $fes[$il] ";
 $cmd = $cmd." --genus-A                    $As[$il] ";
 $cmd = $cmd." --genus-B                    $Bs[$il] ";
 
#$cmd = $cmd." --stack-stat-file            $stack_stat_file ";
#$cmd = $cmd." --ns-history-file            $ns_history_file ";

#$cmd = $cmd." --gmn-hist-file              $gmn_hist_file ";
#$cmd = $cmd." --fcn-data-file              $fcn_data_file ";
 
 $cmd = $cmd." --exit-upon-overflow ";
 $cmd = $cmd." --no-stack-check ";
#$cmd = $cmd." --no-ansi-colors ";
 
 print $cmd;
 
 run_command($job_id, $cmd, "short", "$TmpDir");
};
