
@lambdas = ( 5.0,  4.0,  3.57,  3.45,  3.33,  3.23);
@alphas  = (0.15, 0.14, 0.133,  0.13, 0.127, 0.125);

while(0==0)
{
 for($i=0; $i<scalar(@lambdas); $i++)
 {
  system("G:\\LAT\\sd_metropolis\\bin\\stereo_sm.exe --lambda $lambdas[$i] --alpha $alphas[$i] --average-seq-len 1.5 --average-num-seq 1.1 --no-stack-check --max-alpha-order 7 --prod-mc-steps 5000000 --DIM 2 --LT 32 --LS 32");
 };
};
