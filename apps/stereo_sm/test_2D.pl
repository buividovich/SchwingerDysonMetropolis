
#@lambdas = ( 5.0,  4.0,  3.57,  3.45,  3.33,  3.23);
#@alphas  = (0.15, 0.14, 0.133,  0.13, 0.127, 0.125);

@lambdas = (4.0,  3.23);
@alphas  = (0.16, 0.14);

$counter = 0;

while(0==0 || $counter<1)
{
 for($i=0; $i<scalar(@lambdas); $i++)
 {
  $lambda = $lambdas[$i];
  $alpha  = $lambda/8.0;
  #system("G:\\LAT\\sd_metropolis\\bin\\stereo_sm.exe --lambda $lambdas[$i] --alpha $alphas[$i] --average-seq-len 1.5 --average-num-seq 1.1 --no-stack-check --max-alpha-order 7 --prod-mc-steps 10000000 --DIM 2 --LT 32 --LS 32 --data-dir ./data/stereo_sm --logs-noise-level 1 ");
  system("G:\\LAT\\sd_metropolis\\bin\\stereo_sm.exe --no-stack-check --logs-noise-level 1 --exit-upon-overflow \
  --no-param-auto-tuning                    \
  --lambda              $lambda             \
  --alpha               $alpha              \
  --cc                  1.0                 \
  --NN                  1.0                 \
  --max-alpha-order     7                   \
  --prod-mc-steps       20000000            \
  --DIM                 2                   \
  --LT                  32                  \
  --LS                  32                  \
  --data-dir            ./data/stereo_sm    ");
 };
 $counter ++;
};

#system("drmemory -show_reachable G:\\LAT\\sd_metropolis\\bin\\stereo_smdbg.exe --no-ansi-colors --no-stack-check   --logs-noise-level 1 \
