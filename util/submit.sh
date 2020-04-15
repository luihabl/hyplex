#!/bin/bash

#SBATCH --comment="Testing Hyplex complete run"         # An arbitrary comment on the job
#SBATCH --job-name="hf01"                               # Give this job an arbitrary name
#SBATCH --mail-type=ALL                                 # When to send mail  (BEGIN, END, FAIL, REQUEUE, ALL)
#SBATCH --mail-user=lui.habl@lpp.polytechnique.fr       # Where to send mail.  
#SBATCH --workdir="/home/LPP/lui.habl/hyplex"           # Set current directory before job starts - LOGIN is your login
#SBATCH --error="output/log/%j.err"                         # Direct STDERR here (file identifier), %j is substituted for the job number
#SBATCH --output="output/log/%j.out"                        # Direct STDOUT here (file identifier), %j is substituted for the job number
#SBATCH --verbose                                       # Increase informational messages
#SBATCH --ntasks=1                                      # Number of core (max 64) for your parallel job
#SBATCH --cpus-per-task=3                               # Number of CPUs per task

# Put your email here to receive execution location by mail
echo "" |mail -s "Your job $SLURM_JOB_ID is running on $HOSTNAME" lui.habl@lpp.polytechnique.fr

# Here you have to charge your module

module load openmpi/gcc/64/3.1.4

# Here the code to execute with the SRUN command

mpirun --mca btl self,openib,vader -np $SLURM_NTASKS ./run

# Make sure you return zero as exit code otherwiser SLURM will report your job
# as failed.
exit 0