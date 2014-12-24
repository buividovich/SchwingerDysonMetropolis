$BaseDir    = "G:\\LAT\\sd_metropolis\\GW"; 
$TmpDir     = "C:\\Temp";
$DataDir    = "$BaseDir\\data";

$executable = "$BaseDir\\gw.exe";

system("cd $BaseDir");
#system("make clean");
system("make all");

$mode  = 0;

if($mode==0)
{
 $lambda_min = 4.0;
 $lambda_max = 10.0;
 $dlambda    = 0.5;
};

if($mode==0)
{
 $lambda_min = 0.2;
};

$nmc  = 5000000;
$maxn = 10000;

$suffix = sprintf("mode%i_nmc%i", $mode, $nmc);
$mc_stat_file     = sprintf("%s\\mc_stat_%s.dat", $DataDir, $suffix);
$observables_file = sprintf("%s\\obs_%s.dat",     $DataDir, $suffix);
system("rm -f -v $mc_stat_file");
system("rm -f -v $observables_file");

for($lambda=$lambda_min; $lambda<=$lambda_max; $lambda += $dlambda)
{
 #$ns_hist_file     = sprintf("%s\\ns_hist_l%2.4f_nmc%i.dat", $DataDir, $lambda, $nmc);
 #$ns_history_file  =  sprintf("%s\\ns_history_l%2.4f_nmc%i.dat", $TmpDir, $lambda, $nmc);

 $cmd = $executable;
 $cmd = $cmd." --param-auto-tuning";
 $cmd = $cmd." --mode             $mode";
 $cmd = $cmd." --lambda           $lambda";
 $cmd = $cmd." --nmc              $nmc";
 $cmd = $cmd." --maxn             $maxn";
 $cmd = $cmd." --mc-stat-file     $mc_stat_file";
 #$cmd = $cmd." --ns-hist-file     $ns_hist_file";
 $cmd = $cmd." --observables-file $observables_file";
 #$cmd = $cmd." --ns-history-file  $ns_history_file";
 
 system("$cmd");
};
