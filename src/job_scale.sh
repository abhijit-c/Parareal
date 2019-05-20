#!/bin/sh
#
#SBATCH --verbose
#SBATCH --job-name=ac6361_parareal
#SBATCH --output=parareal_%j.out
#SBATCH --err=parareal_%j.err
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=4GB
#SBATCH --mail-type=END
#SBATCH --mail-user=ac6361@nyu.edu
 
srun ./linear1d_test.out
