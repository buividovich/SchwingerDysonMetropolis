#!/usr/bin/perl

sub run_app{
 my @mc_stat_files = ();
 my $cmd0 = $executable;
 
 $cmd0 = $cmd0." --lambda                   $lambda ";
 $cmd0 = $cmd0." --meff-sq                  $meff_sq ";
 $cmd0 = $cmd0." --LT                       $LT ";
 $cmd0 = $cmd0." --prod-mc-steps            $nmc ";
 $cmd0 = $cmd0." --mc-reporting-interval    $mc_reporting_interval ";
 $cmd0 = $cmd0." --observables-file         $observables_file ";

 my $njobs = (int(scalar(@_)/2)>1? int(scalar(@_)/2) : 1);
 print "I am going to run $njobs jobs...\n";
 
 for(my $ijob=0; $ijob<$njobs; $ijob++)
 {
  my $cmd = $cmd0;
  if(scalar(@_)>=2)
  {
   $cmd = $cmd." --no-param-auto-tuning ";
   $cmd = $cmd." --cc $_[2*$ijob] ";
   $cmd = $cmd." --NN $_[2*$ijob+1] ";
   $mc_stat_file = sprintf("%s/l%2.4lf_cc%2.4lf_NN%2.4lf.dat", $TmpDir, $lambda, $_[2*$ijob+0], $_[2*$ijob+1]);
  }
  else
  {
   $mc_stat_file = sprintf("%s/l%2.4lf.dat", $TmpDir, $lambda, $_[0], $_[1]);
  };
  $cmd = $cmd." --mc-stat-file $mc_stat_file ";
 
  system("rm -f -v $mc_stat_file");
 
  system("$NohupCmdStart $cmd $NohupCmdEnd");
  
  push(@mc_stat_files, $mc_stat_file);
 };
 
 #Now wait until all the mc_stat_files exist and all the locks are removed
 my $all_files_created = 1;
 my @lock_files = (); 
 my $wait_count = 0;
 do{
  sleep(1);
  $all_files_created = 1;
  for(my $ijob=0; $ijob<scalar(@mc_stat_files); $ijob++)
  {
   $all_files_created = $all_files_created & (-e $mc_stat_files[$ijob]);
  };
  @lock_files = glob("$DataDir/*.lock");
  $wait_count ++;
 }while(($all_files_created==0 | scalar(@lock_files)!=0) & $wait_count<100);
 
 #And now read the numerical values from all the files
 my @mc_stat_data = ();
 for(my $ijob=0; $ijob<scalar(@mc_stat_files); $ijob++)
 {
  open(MC_STAT,"<$mc_stat_files[$ijob]");
  $row_count = 0;
  while(<MC_STAT>)
  {
   @cols = split;
   if($row_count==0 && scalar(@cols)>=10)
   {
    push(@mc_stat_data, $cols[2]); #cc
    push(@mc_stat_data, $cols[3]); #NN
    push(@mc_stat_data, $cols[8]); #nA
    push(@mc_stat_data, $cols[9]); #dnA
   }
   else
   {
    print "\n\t >>>>>>>> Something is strange: 1st line in file $mc_stat_file has less than 9 columns\n\n";
   };
   $row_count ++;
  };
  if($row_count>1)
  {
   print "\n\t >>>>>>>> Something is strange: more than 1 line in file $mc_stat_file\n\n";
  };
  close(MC_STAT); 
 };
 
 return @mc_stat_data;
};

sub tune_cc_NN
{
 #initial run with auto-tuning to get initial values
 my @mc_stat_data = run_ssm();
 my $cc0 = $mc_stat_data[0];
 my $NN0 = $mc_stat_data[1];
 my $v0  = $mc_stat_data[2];

 my $sfinish = 0;
 while(!$sfinish)
 {
  @mc_stat_data = run_ssm(0.95*$cc0, $NN0, 1.05*$cc0, $NN0, $cc0, 0.95*$NN0, $cc0, 1.05*$NN0);
  my $minsdir = 0;
  my $minv0   = $mc_stat_data[4*0 + 2];
  for(my $sdir=0; $sdir<4; $sdir++)
  {
   printf("At cc = %2.4lf, NN = %2.4lf, nA = %2.4lf +/- %2.4lf\n", $mc_stat_data[4*$sdir + 0], $mc_stat_data[4*$sdir + 1], $mc_stat_data[4*$sdir + 2], $mc_stat_data[4*$sdir + 3]);
   if($mc_stat_data[4*$sdir + 2]<$minv0)
   {
    $minv0 = $mc_stat_data[4*$sdir + 2];
    $minsdir = $sdir;
   };
  };
  if($minv0 < $v0)
  {
   print "\n\t => New search direction: $minsdir, (nA ~ $minv0)\n\n";
   $cc0 = $mc_stat_data[4*$minsdir + 0];
   $NN0 = $mc_stat_data[4*$minsdir + 1];
   $v0  = $mc_stat_data[4*$minsdir + 2];
  }
  else
  {
   $sfinish = 1;
  };
 };
 
 return ($cc0, $NN0, $v0); 
};

sub lambda_descent
{
 my $v0, $dv0;
 my $cc0 = $_[0];
 my $NN0 = $_[1];
 printf("\n\n ##### Starting lambda descent with cc=%2.4lf, NN=%2.4lf\n\n", $cc0, $NN0);
 do{
  $lambda *= 0.9;
  my ($cc0, $NN0, $v0) = tune_cc_NN();
  printf("With lambda = %2.4lf, <nA> = %2.4lf\n", $lambda, $v0);
 }while($v0<0.8);
 
 return $v0;
}

return 1;
