#!/bin/bash
#SBATCH --job-name=domainwall
#SBATCH --partition=qc-gh200
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -t 00:20:00

module load system/qc-gh200  nvhpc-hpcx-cuda12/24.3 

echo
echo "==================================================================="
echo
cat main_template.yaml |sed -e "s/xxx_lattice_size_xxx/32,8,8,48/" -e "s/xxx_grid_size_xxx/1,1,1,4/" > main.yaml
sh ./run.sh

exit

echo
echo "==================================================================="
echo
cat main_template.yaml |sed -e "s/xxx_lattice_size_xxx/64,32,16,48/" -e "s/xxx_grid_size_xxx/1,1,1,4/" > main.yaml
sh ./run.sh


echo
echo "==================================================================="
echo
cat main_template.yaml |sed -e "s/xxx_lattice_size_xxx/32,32,32,256/" -e "s/xxx_grid_size_xxx/1,1,1,4/" > main.yaml
sh ./run.sh
