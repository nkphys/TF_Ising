dt=0.001
site_1=15
file="SzSzCorr_dt0.001.txt"

grep "^0   ${site_1} " ${file}

for ts in {0..2002}
do 

time=$(echo "${dt}*${ts}" | bc -l)
time_3dp=$(printf "%1.3f" ${time})
#time_2dp=$(printf "%1.2f" ${time})
grep "^${time_3dp}   ${site_1} " ${file}
#grep "^${time_2dp}   ${site_1} " ${file}

done
