
LTs="4 8 12 16 20 24 28 32 34 36 38 40 42 48 54 60 66 72 90"
lambdas="3.0120 3.1000 3.2300 3.3300 3.4500 3.5700 4.0000 5.0000"

rm -v /g/LAT/sd_metropolis/data/pcm_wc_mspace/*summary*.dat

#for LT in $LTs
#do
#  suffix="d2_t"$LT"_s108_l3.0120"
#  echo -e "\x1b[1;36m $suffix \x1b[0m"
#  /g/LAT/sd_metropolis/bin/pcmwcdp.exe \
#  --data-suffix                $suffix \
#  --scan-suffix                     LT \
#  --scan-label                     $LT \
#  --max-order                       11 \
#  --LS                             108 \
#  --calculate-correlators              \
#  --cluster-label                  idc \
#  --cluster-label                 idc1 \
#  --cluster-label                itep1 \
#  --cluster-label                itep2
#done
#
#for lambda in $lambdas
#do
#  suffix="d2_t108_s108_l"$lambda
#  echo -e "\x1b[1;36m $suffix \x1b[0m"
#  /g/LAT/sd_metropolis/bin/pcmwcdp.exe \
#  --data-suffix                $suffix \
#  --scan-suffix                 lambda \
#  --scan-label                 $lambda \
#  --max-order                       11 \
#  --LS                             108 \
#  --calculate-correlators              \
#  --cluster-label                  idc \
#  --cluster-label                 idc1 \
#  --cluster-label                itep1 \
#  --cluster-label                itep2
#done

#And now the special large-volume run

lambda=3.1000
suffix="d2_t256_s256_l"$lambda
echo -e "\x1b[1;36m $suffix \x1b[0m"
/g/LAT/sd_metropolis/bin/pcmwcdp.exe \
--data-suffix                $suffix \
--scan-suffix                 lambda \
--scan-label                 $lambda \
--max-order                       11 \
--LS                             256 \
--cluster-label                vol00 \
--cluster-label                vol01 \
--cluster-label                vol02 \
--cluster-label                vol03 \
--cluster-label                vol04 \
--cluster-label                vol05 \
--cluster-label                vol06 \
--cluster-label                vol07 \
--cluster-label                vol08 \
--cluster-label                vol09 \
--cluster-label                vol11 \
--cluster-label                vol12 \
--cluster-label                vol13 \
--cluster-label                vol14 \
--cluster-label                vol15 \
--cluster-label                vol16 \
--cluster-label                vol17 \
--cluster-label                vol18 \
--cluster-label                vol19 \
--cluster-label                vol20 \
--cluster-label                vol21 \
--cluster-label                vol22 \
--cluster-label                vol23 \
--cluster-label                vol24 \
--cluster-label                vol25 \
--cluster-label                vol26 \
--cluster-label                vol27 \
--cluster-label                vol28 \
--cluster-label                vol29 \
--cluster-label                vol30 \
--cluster-label                vol31 \
--cluster-label                vol32 \
--cluster-label                vol33 \
--cluster-label                vol34 \
--cluster-label                vol35 \
--cluster-label                vol36 \
--cluster-label                vol37 \
--cluster-label                vol38 \
--cluster-label                vol39

#--calculate-correlators              \
