import os
import subprocess

# INPUT_DIRECTORY = 'data/750ns'
# OUTPUT_DIRECTORY = 'results/750ns'
INPUT_DIRECTORY = 'data'
OUTPUT_DIRECTORY = 'results'

def execute_command(command):
    try:
        subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"Command execution failed: {e}")

for dirpath, dirnames, filenames in os.walk(INPUT_DIRECTORY):    
    for dirname in dirnames:
        if dirname.startswith('apo_'):
            dir_in = os.path.join(dirpath, dirname)
            dir_out = os.path.join(OUTPUT_DIRECTORY, dirname)
            command = f'analyse.py -i {dir_in} -o {dir_out}'
            print(f'-'*80)
            print(f'-------- {command} -------- ')
            print(f'-'*80)
            execute_command(command)

