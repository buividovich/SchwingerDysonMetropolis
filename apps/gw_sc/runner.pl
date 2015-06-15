$BaseDir    = "G:\\LAT\\sd_metropolis\\apps\\gw_sc"; 
$TmpDir     = "C:\\Temp";
$DataDir    = "G:\\LAT\\sd_metropolis\\data\\gw_sc";

$executable = "G:\\LAT\\sd_metropolis\\bin\\gw_sc.exe";

system("cd $BaseDir\\..\\..");
system("make gw_sc");
system("cd $BaseDir");

$lambda_min = 5.0;
$lambda_max =  10.0;
$dlambda    =  0.5;

$nmc  = 20000000;
$maxn = 20000;

$mc_stat_file = "$DataDir\\mc_stat_nmc$nmc.dat";
system("rm -f -v $mc_stat_file");

$observables_file = sprintf("%s\\G_nmc%i.dat", $DataDir, $nmc);
system("rm -f -v $observables_file");


for($lambda=$lambda_min; $lambda<=$lambda_max; $lambda += $dlambda)
{
 $ns_history_file  =  sprintf("%s\\ns_history_l%2.4f_nmc%i.dat", $TmpDir, $lambda, $nmc);

 $gtotal_file = sprintf("%s\\T_l%2.4f_nmc%i.dat", $DataDir, $lambda, $nmc);
 system("rm -f -v $gtotal_file");

 $cmd = $executable;
 $cmd = $cmd." --lambda                   $lambda";
 $cmd = $cmd." --prod-mc-steps            $nmc";
 $cmd = $cmd." --max-recursion-depth      $maxn";
 $cmd = $cmd." --max-correlator-order     10";
 $cmd = $cmd." --mc-stat-file             $mc_stat_file";
 $cmd = $cmd." --gtotal-file              $gtotal_file";
 $cmd = $cmd." --observables-file         $observables_file";
 $cmd = $cmd." --ns-history-file          $ns_history_file";
 $cmd = $cmd." --logs-noise-level         1";
 $cmd = $cmd." --mc-reporting-interval    50000";
 $cmd = $cmd." --no-stack-check ";
 $cmd = $cmd." --no-ansi-colors ";
 
 system("$cmd");
};
