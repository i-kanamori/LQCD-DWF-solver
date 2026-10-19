#!/bin/bash
#SBATCH --job-name=domainwall
#SBATCH --partition=qc-gh200
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 01:00:00

module load system/qc-gh200 nvhpc-hpcx/24.9

echo
echo "==================================================================="
echo
NX=8
NY=8
NZ=8
NT=8
cat main_template.yaml |sed -e "s/xxx_lattice_size_xxx/$NX,$NY,$NZ,$NT/" -e "s/xxx_grid_size_xxx/1,1,1,1/" > main.yaml
output=run_${NX}x${NY}x${NZ}x${NT}.log
sh ./run.sh >& $output

echo
echo "==================================================================="
echo
NX=16
NY=16
NZ=16
NT=16
cat main_template.yaml |sed -e "s/xxx_lattice_size_xxx/$NX,$NY,$NZ,$NT/" -e "s/xxx_grid_size_xxx/1,1,1,1/" > main.yaml
output=run_${NX}x${NY}x${NZ}x${NT}.log
sh ./run.sh >& $output

echo
echo "==================================================================="
echo
NX=24
NY=24
NZ=24
NT=24
cat main_template.yaml |sed -e "s/xxx_lattice_size_xxx/$NX,$NY,$NZ,$NT/" -e "s/xxx_grid_size_xxx/1,1,1,1/" > main.yaml
output=run_${NX}x${NY}x${NZ}x${NT}.log
sh ./run.sh >& $output

echo
echo "==================================================================="
echo
NX=32
NY=32
NZ=32
NT=32
cat main_template.yaml |sed -e "s/xxx_lattice_size_xxx/$NX,$NY,$NZ,$NT/" -e "s/xxx_grid_size_xxx/1,1,1,1/" > main.yaml
output=run_${NX}x${NY}x${NZ}x${NT}.log
sh ./run.sh >& $output

echo
echo "==================================================================="
echo
NX=32
NY=32
NZ=32
NT=64
cat main_template.yaml |sed -e "s/xxx_lattice_size_xxx/$NX,$NY,$NZ,$NT/" -e "s/xxx_grid_size_xxx/1,1,1,1/" > main.yaml
output=run_${NX}x${NY}x${NZ}x${NT}.log
sh ./run.sh >& $output


exit

echo
echo "==================================================================="
echo
NX=64
NY=64
NZ=64
NT=64
cat main_template.yaml |sed -e "s/xxx_lattice_size_xxx/$NX,$NY,$NZ,$NT/" -e "s/xxx_grid_size_xxx/1,1,1,1/" > main.yaml
output=run_${NX}x${NY}x${NZ}x${NT}.log
sh ./run.sh >& $output

