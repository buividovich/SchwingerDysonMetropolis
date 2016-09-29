$lambda = 6.0;
$beta   = 1/$lambda;

$counter = 0;

$NN = 1.0;
$a  = $lambda;

while(1==1 || $counter<1)
{ 
 foreach $cc (1.0)
 {
  $label = sprintf("a%2.2f_c%2.2f_n%2.2f", $a, $cc, $NN);
  $alpha = $beta*$a;
  
  system("G:\\LAT\\sd_metropolis\\bin\\pcm_sc_cspace.exe --logs-noise-level 1 --exit-upon-overflow \
  --lambda              $lambda                     \
  --alpha               $alpha                      \
  --cc                  $cc                         \
  --NN                  $NN                         \
  --max-order           5                           \
  --number-mc-steps     50000000                    \
  --DIM                 2                           \
  --LT                  8                           \
  --LS                  8                           \
  --data-dir            ./data/pcm_sc_cspace/       \
  --label               $label                      \
  --save-sampling-hist                              ");
 };
}; 

