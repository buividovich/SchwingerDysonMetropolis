$BaseDir    = "G:\\LAT\\sd_metropolis\\apps\\hmm"; 
$TmpDir     = "C:\\Temp";
$DataDir    = "G:\\LAT\\sd_metropolis\\data\\hmm";

$executable = "G:\\LAT\\sd_metropolis\\bin\\hmm.exe";

system("cd $BaseDir\\..\\..");
system("make hmm");
system("cd $BaseDir");

$lambda = 0.08;

$nmc  = 20000000;
$maxn = 20000;

$mc_stat_file = sprintf("%s\\mc_stat_acctest_nmc%i_l%2.4f.dat", $DataDir, $nmc, $lambda);
system("rm -f -v $mc_stat_file");

for($pplus=0.05; $pplus<=0.95; $pplus += 0.05)
{
 $cmd = $executable;
 $cmd = $cmd." --lambda              $lambda";
 $cmd = $cmd." --prod-mc-steps       $nmc";
 $cmd = $cmd." --max-recursion-depth $maxn";
 $cmd = $cmd." --mc-stat-file        $mc_stat_file";
 $cmd = $cmd." --logs-noise-level         0";
 $cmd = $cmd." --mc-reporting-interval    50000";
 $cmd = $cmd." --p-plus              $pplus";
 
 system("$cmd");
};
