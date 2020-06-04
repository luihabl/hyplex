EXP_NAME=$1
sbatch --job-name=${EXP_NAME} --export=EXP_NAME=${EXP_NAME} submit.sh
