$lambda = 3.23;

while(1==1)
{
foreach $a (0.5, 1.0, 1.5)
{
 $alpha  = $a*$lambda/8.0;
 foreach $cc (0.5, 1.0, 1.5)
 {
  foreach $NN (0.5, 1.0, 1.5)
  {
   $label = sprintf("a%2.2f_c%2.2f_n%2.2f", $a, $cc, $NN);
   system("G:\\LAT\\sd_metropolis\\bin\\stereo_sm.exe --logs-noise-level 1 --exit-upon-overflow \
     --lambda              $lambda             \
     --alpha               $alpha              \
     --cc                  $cc                 \
     --NN                  $NN                 \
     --max-order           9                   \
     --number-mc-steps     20000000            \
     --DIM                 2                   \
     --LT                  32                  \
     --LS                  32                  \
     --data-dir            ./data/stereo_sm/   \
     --label               $label              \
     --save-sampling-hist                      ");
  };
 };
};
};

