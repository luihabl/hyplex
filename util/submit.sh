#!/bin/bash

#SBATCH --mail-type=ALL                                 # When to send mail  (BEGIN, END, FAIL, REQUEUE, ALL)
#SBATCH --mail-user=lui.habl@lpp.polytechnique.fr       # Where to send mail.  
#SBATCH --workdir="/home/LPP/lui.habl/hyplex"           # Set current directory before job starts - LOGIN is your login
#SBATCH --error="output/log/%j.err"                     # Direct STDERR here (file identifier), %j is substituted for the job number
#SBATCH --output="output/log/%j.out"                    # Direct STDOUT here (file identifier), %j is substituted for the job number
#SBATCH --verbose                                       # Increase informational messages
#SBATCH --ntasks=12                                      # Number of tasks (max 64) for your parallel job

# Put your email here to receive execution location by mail
echo "" |mail -s "Your job $EXP_NAME ($SLURM_JOB_ID) is running on $HOSTNAME" lui.habl@lpp.polytechnique.fr

# Here you have to charge your module

module purge
module load slurm/14.03.0
module load gcc/6.3.0
module load make/4.2
module load openmpi/gcc/64/3.1.4

# Here the code to execute with the SRUN command

mpirun --bind-to none --mca btl openib,self,vader -np $SLURM_NTASKS ./run --config=input/config/config-$EXP_NAME.yaml

exit 0
