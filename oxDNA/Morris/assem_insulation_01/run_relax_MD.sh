
#!/bin/sh
#PBS -lwalltime=00:20:00
#PBS -lselect=1:ncpus=1:mem=10GB

cd /rds/general/user/tc2523/home/assem_insulation/
inputfile="input_relax_MD.txt"

runcommand="/rds/general/user/tc2523/home/assem_insulation/  $inputfile"
timeout 15m $runcommand

