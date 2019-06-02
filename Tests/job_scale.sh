#!/bin/sh
#
#SBATCH --verbose
#SBATCH --job-name=ac6361_parareal
#SBATCH --output=parareal_%j.txt
#SBATCH --err=parareal_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=4GB
#SBATCH --mail-type=END
#SBATCH --mail-user=ac6361@nyu.edu
 
srun ./strongscale_test.out
srun ./weakscale_test.out
