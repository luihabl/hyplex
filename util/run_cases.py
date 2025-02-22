from pathlib import Path
import os, argparse

def run_all(dir_with_files, run_file):

    run_file = Path(run_file)
    dir_with_files = Path(dir_with_files).absolute()
    config_files = list(dir_with_files.glob('*.yaml')) 

    for case in config_files:
        print(f"running: {case}")
        os.chdir(run_file.parent)
        print(f'mpirun -np 1 {run_file.name} -c {case}')
        os.system(f'mpirun -np 1 {run_file.name} -c {case}')
        os.system('sleep 1')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run Hyplex in sequence with all the config files in a directory.')
    parser.add_argument("-c", "--config", help="directory with the config files.", default="config/")
    parser.add_argument("-r", "--run", help="Hyplex executable", default= "../run")
    args = parser.parse_args()

    run_all(args.config, args.run)






