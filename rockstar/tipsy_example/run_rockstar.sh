#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks=80
#SBATCH --account=rrg-wadsley-ab
#SBATCH --time=02:00:00

cd $SLURM_SUBMIT_DIR

module load NiaEnv/2019b gcc/8.3.0 openmpi/4.0.1

rm test_output.dat

./rockstar-galaxies -c rockstar-galaxies.cfg &> server.dat &
mpirun -N 1 ./rockstar-galaxies -c auto-rockstar.cfg >> test_output.dat
