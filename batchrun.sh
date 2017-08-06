#!/bin/bash
# cores per node
cpn=20
# list of executables
for x in mc_p ge_p ge_pg; do
# for each number of procs
 for i in 2 4 8 16 32 64 128; do
            if [[ $i -lt $cpn ]]; then
                ppn=$i
            else
                ppn=$cpn
            fi
            #replace parameters in template sbatch file
            sed -i "s/\#SBATCH --job-name=.*/\#SBATCH --job-name=$x/;s/\#SBATCH --ntasks=.*/\#SBATCH --ntasks=$i/;s/\#SBATCH --ntasks-per-node=.*/\#SBATCH --ntasks-per-node=$ppn/;s/^N_PROC=.*/N_PROC=$i/;s/^BIN_NAME=.*/BIN_NAME=$x/" comp.sbatch
            #submit job
            sbatch comp.sbatch
 done
done

