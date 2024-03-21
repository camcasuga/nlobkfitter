#!/bin/sh
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 2:00:00
#SBATCH -o log.txt
#SBATCH --partition=small
#SBATCH --account=lappi
#SBATCH --mem-per-cpu=10M

#### Argument order: SCHEME BK RC useImprovedZ2Bound useBoundLoop [Qs0 C^2 gamma] X0_if X0_bk e_c Q0sq Y0 eta0 old_sigma0 mass_charm mass_bottom | tee outputfile.dat
#### ${1} = Qs0^2;  ${2} = C2;   ${3} = sigma0/2;  ${4} = gamma; ${5} = output file name

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/users/caaucasu/nlobkfitter/Minuit2/lib
export LIBRARY_PATH=${LIBRARY_PATH}:/users/caaucasu/nlobkfitter/Minuit2/lib
module add cmake/3.23.1
module add gsl

OMP_NUM_THREADS=12 CUBACORES=0 ./build/bin/nlofit unlpc resumbk pdrc z2imp unb ${1} ${2} ${4} 1.0 1.0 1 1 0 0 ${3} 1.4 4.9 | tee ./sigmars_nlo1/${5}.dat
