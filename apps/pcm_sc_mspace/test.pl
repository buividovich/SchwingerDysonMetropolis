
$lambda = 4.0;
$beta   = 1/$lambda;

$counter = 0;

while(1==1 || $counter<1)
{
 foreach $m0 (2.6, 2.8)
 {
  $meff_sq = -4.0 + $m0;
  system("G:\\LAT\\sd_metropolis\\bin\\pcm_sc_mspace.exe --logs-noise-level 1 --exit-upon-overflow \
   --lambda              $lambda                     \
   --meff-sq             $meff_sq                    \
   --alpha               $beta                       \
   --cc                  4.0                         \
   --NN                  1.0                         \
   --max-order           6                           \
   --number-mc-steps     50000000                    \
   --DIM                 2                           \
   --LT                  8                           \
   --LS                  8                           \
   --data-dir            ./data/pcm_sc_mspace/       \
   --resummation                                     \
   --save-sampling-hist                              ");
 };
 $counter ++;  
};
