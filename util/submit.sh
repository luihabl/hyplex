#!/bin/bash

#SBATCH --comment="Running nominal and dc cases"        # An arbitrary comment on the job
#SBATCH --job-name="nom-dc"                             # Give this job an arbitrary name
#SBATCH --mail-type=ALL                                 # When to send mail  (BEGIN, END, FAIL, REQUEUE, ALL)
#SBATCH --mail-user=lui.habl@lpp.polytechnique.fr       # Where to send mail.  
#SBATCH --workdir="/home/LPP/lui.habl/hyplex"           # Set current directory before job starts - LOGIN is your login
#SBATCH --error="output/log/%j.err"                     # Direct STDERR here (file identifier), %j is substituted for the job number
#SBATCH --output="output/log/%j.out"                    # Direct STDOUT here (file identifier), %j is substituted for the job number
#SBATCH --verbose                                       # Increase informational messages

#SBATCH --ntasks=4                                      # Number of tasks (max 64) for your parallel job
#SBATCH --cpus-per-task=2                               # number of cores per task (e.g. OpenMP)

# Put your email here to receive execution location by mail
echo "" |mail -s "Your job $SLURM_JOB_ID is running on $HOSTNAME" lui.habl@lpp.polytechnique.fr

# Here you have to charge your module

module purge
module load slurm/14.03.0
module load gcc/6.3.0
module load make/4.2
module load openmpi/gcc/64/3.1.4

# Here the code to execute with the SRUN command

mpirun --bind-to none --mca btl openib,self,vader -np $SLURM_NTASKS ./run --config=input/config/config-nominal.yaml


#mpirun --mca btl self,openib,vader -np $SLURM_NTASKS ./run -n 1 &
#mpirun --mca btl self,openib,vader -np $SLURM_NTASKS ./run -n 2 &
#mpirun --bind-to none --mca btl openib,self,vader -np 1 ./run

#parallel --max-procs=${SLURM_CPUS_PER_TASK} "$cmd,subjob={1}; sleep 30" ::: {1..5}

#for sj in {0..2}
#do
#    mpirun --bind-to none --mca btl openib,self,vader -np 1 ./run -n $sj > output/log/${SLURM_JOB_ID}_${sj}.out &
#done
#wait

#mpi_cmd="mpirun --bind-to none --mca btl openib,self,vader -np 1"
#parallel_cmd="parallel --delay .2 -j $SLURM_NTASKS"
#$parallel_cmd "$mpi_cmd ./run{1} -n > output/log/${SLURM_JOB_ID}_{1}.out" ::: {0..2}

# mpirun --bind-to none --mca btl openib,self,vader -np $SLURM_NTASKS ./run
# Make sure you return zero as exit code otherwiser SLURM will report your job
# as failed.
exit 0
