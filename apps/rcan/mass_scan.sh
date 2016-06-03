rm -f -v ./data/rcan/mscan_l2.3100.dat

for i in $(seq -1.5 0.005 0.5)
do
 ./bin/rcan1.exe --lambda 2.31 --meff-sq $i --mmax 12 --append-mode --no-check
done
