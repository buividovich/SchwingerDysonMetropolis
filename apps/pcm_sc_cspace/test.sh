while true
do
 G:\\LAT\\sd_metropolis\\bin\\pcm_sc_mspace.exe --logs-noise-level 1 --exit-upon-overflow \
  --lambda              4.0                 \
  --meff-sq            -1.0                 \
  --alpha               0.25                \
  --cc                  1.0                 \
  --NN                  1.0                 \
  --max-order           4                   \
  --number-mc-steps     50000000            \
  --DIM                 1                   \
  --LT                  2                   \
  --LS                  1                   \
  --data-dir            ./data/sc_sm/       \
  --resummation                             \
  --save-sampling-hist                      
done
