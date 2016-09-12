
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
  
  #system("G:\\LAT\\sd_metropolis\\bin\\stereo_sm.exe --lambda $lambdas[$i] --alpha $alphas[$i] --average-seq-len 1.5 --average-num-seq 1.1 --no-stack-check --max-alpha-order 7 --prod-mc-steps 10000000 --DIM 2 --LT 32 --LS 32 --data-dir ./data/stereo_sm --logs-noise-level 1 ");
  system("G:\\LAT\\sd_metropolis\\bin\\stereo_sm.exe --logs-noise-level 1 --exit-upon-overflow \
  --lambda              $lambda             \
  --alpha               $alpha              \
  --cc                  $cc                 \
  --NN                  $NN                 \
  --max-order           11                  \
  --number-mc-steps     20000000            \
  --DIM                 2                   \
  --LT                  32                  \
  --LS                  32                  \
  --data-dir            ./data/stereo_sm/   \
  --save-sampling-hist                      \
  --save-correlators                        ");
 };
 $counter ++;
};

#system("drmemory -show_reachable G:\\LAT\\sd_metropolis\\bin\\stereo_smdbg.exe --no-ansi-colors --logs-noise-level 1 \
