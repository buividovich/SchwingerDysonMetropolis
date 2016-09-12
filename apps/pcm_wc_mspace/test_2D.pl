
#@lambdas = ( 5.0,  4.0,  3.57,  3.45,  3.33,  3.23);
#@lambdas = (4.0,  3.23);

@lambdas = (3.23);

$counter = 0;

while(0==0 || $counter<1)
{
 for($i=0; $i<scalar(@lambdas); $i++)
 {
  $lambda = $lambdas[$i];
  
  $alpha  = 1.0*($lambda/8.0);
  $cc     = 1.0;
  $NN     = 0.5;
  
  system("G:\\LAT\\sd_metropolis\\bin\\pcm_wc_mspace.exe --logs-noise-level 1 --exit-upon-overflow \
  --lambda              $lambda                 \
  --alpha               $alpha                  \
  --cc                  $cc                     \
  --NN                  $NN                     \
  --max-order           11                      \
  --number-mc-steps     20000000                \
  --DIM                 2                       \
  --LT                  32                      \
  --LS                  32                      \
  --data-dir            ./data/pcm_wc_mspace/   \
  --save-sampling-hist                          \
  --save-correlators                            ");
 };
 $counter ++;
};

