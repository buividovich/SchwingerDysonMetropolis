while true
do
 G:\\LAT\\sd_metropolis\\bin\\pcm_sc_cspace.exe --logs-noise-level 1 --exit-upon-overflow \
  --lambda              5.0                         \
  --alpha               0.20                        \
  --cc                  1.0                         \
  --NN                  0.5                         \
  --max-order           2                           \
  --number-mc-steps     50000000                    \
  --DIM                 1                           \
  --LT                  3                           \
  --LS                  1                           \
  --data-dir            ./data/pcm_sc_cspace/       \
  --save-sampling-hist
done
