#!/usr/bin/perl

system("cd $GlobalCodeDir/sd_metropolis/");

$AppDir     = "$GlobalCodeDir/sd_metropolis/apps/$app"; 
$DataDir    = "$GlobalDataDir/sd_metropolis/data/$app";
$TmpDir     = "$GlobalCodeDir/sd_metropolis/tmp";

$executable = "$GlobalCodeDir/sd_metropolis/bin/$app";
$executable = $executable." --no-stack-check ";
$executable = $executable." --logs-noise-level 1 ";
$executable = $executable." --no-ansi-colors ";
$executable = $executable." --io-sleep-time 0.5 ";
$executable = $executable." --io-write-attempts 20 ";

if($myOS eq "Linux")
{
  $executable = $executable." --print-errors-to-stderr ";
};

#system("make $app");
system("$OutputCleanupCmd");

return 1;

