$BaseDir    = "G:\\LAT\\sd_metropolis\\apps\\hmm"; 
$TmpDir     = "C:\\Temp";
$DataDir    = "G:\\LAT\\sd_metropolis\\data\\hmm";

$executable = "G:\\LAT\\sd_metropolis\\bin\\hmm.exe";

system("cd $BaseDir\\..\\..");
system("make hmm");
system("cd $BaseDir");

$lambda_min = -0.09;
$lambda_max =  0.09;
$dlambda    =  0.002;

$nmc  = 5000000;
$maxn = 20000;

$mc_stat_file = "$DataDir\\mc_stat_nmc$nmc.dat";
system("rm -f -v $mc_stat_file");

$observables_file = sprintf("%s\\G_nmc%i.dat", $DataDir, $nmc);
system("rm -f -v $observables_file");

for($lambda=$lambda_min; $lambda<=$lambda_max; $lambda += $dlambda)
{
 $ns_history_file  =  sprintf("%s\\ns_history_l%2.4f_nmc%i.dat", $TmpDir, $lambda, $nmc);

 $cmd = $executable;
 $cmd = $cmd." --lambda              $lambda";
 $cmd = $cmd." --prod-mc-steps       $nmc";
 $cmd = $cmd." --max-recursion-depth $maxn";
 $cmd = $cmd." --mc-stat-file        $mc_stat_file";
 $cmd = $cmd." --observables-file    $observables_file";
 $cmd = $cmd." --ns-history-file     $ns_history_file";
 $cmd = $cmd." --logs-noise-level         0";
 $cmd = $cmd." --mc-reporting-interval    50000";
 
 system("$cmd");
};
