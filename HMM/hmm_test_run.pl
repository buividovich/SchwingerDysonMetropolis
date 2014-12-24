$BaseDir    = "G:\\LAT\\sd_metropolis\\HMM"; 
$TmpDir     = "C:\\Temp";
$DataDir    = "$BaseDir\\data";

$executable = "$BaseDir\\hmm.exe";

system("cd $BaseDir");
system("make clean");
system("make all");

$lambda_min = -0.1;
$lambda_max =  0.1;
$dlambda    =  0.002;

$nmc  = 5000000;
$maxn = 20000;

$mc_stat_file = "$DataDir\\mc_stat_nmc$nmc.dat";
system("rm -f -v $mc_stat_file");

for($lambda=$lambda_min; $lambda<=$lambda_max; $lambda += $dlambda)
{
 $ns_hist_file     = sprintf("%s\\ns_hist_l%2.4f_nmc%i.dat", $DataDir, $lambda, $nmc);
 $observables_file = sprintf("%s\\G_l%2.4f_nmc%i.dat", $DataDir, $lambda, $nmc);
 #$ns_history_file  =  sprintf("%s\\ns_history_l%2.4f_nmc%i.dat", $TmpDir, $lambda, $nmc);

 $cmd = $executable;
 $cmd = $cmd." --param-auto-tuning";
 $cmd = $cmd." --lambda           $lambda";
 $cmd = $cmd." --nmc              $nmc";
 $cmd = $cmd." --maxn             $maxn";
 $cmd = $cmd." --mc-stat-file     $mc_stat_file";
 $cmd = $cmd." --ns-hist-file     $ns_hist_file";
 $cmd = $cmd." --observables-file $observables_file";
 #$cmd = $cmd." --ns-history-file  $ns_history_file";
 
 system("$cmd");
};
