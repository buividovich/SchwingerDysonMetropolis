
LTs="4 8 12 16 20 24 28 32 34 36 38 40 42 48 54 60 66 72 90"
lambdas="3.0120 3.1000 3.2300 3.3300 3.4500 3.5700 4.0000 5.0000"

rm -v /g/LAT/sd_metropolis/data/pcm_wc_mspace/*summary*.dat

for LT in $LTs
do
  suffix="d2_t"$LT"_s108_l3.0120"
  echo -e "\x1b[1;36m $suffix \x1b[0m"
  /g/LAT/sd_metropolis/bin/pcmwcdp.exe --data-suffix $suffix --scan-suffix LT --scan-label $LT
done

for lambda in $lambdas
do
  suffix="d2_t108_s108_l"$lambda
  echo -e "\x1b[1;36m $suffix \x1b[0m"
  /g/LAT/sd_metropolis/bin/pcmwcdp.exe --data-suffix $suffix --scan-suffix lambda --scan-label $lambda
done

