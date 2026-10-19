date
echo "--------------- this script ----------------"
cat $0
echo "--------------- this script, done ----------"
echo
module list

export OMP_NUM_THREADS=2
echo OMP_NUM_THREAD: $OMP_NUM_THREADS
export OMP_PROC_BIND=spread
env

echo "--------------- main.yaml ------------------"
cat main.yaml
echo "--------------- main.yaml, done-------------"
echo
INPUT=test_alt_Spectrum_Domainwall_5din_Hadron2ptFunction.yaml
echo $INPUT
echo "--------------- INPUT ----------------------"
cat $INPUT
echo "--------------- INPUT, done ----------------"
echo
BIN=../benchmark/domainwall/bridge.elf
ls -l $BIN
date
mpirun -bind-to none $BIN alt_accel
date
echo done


