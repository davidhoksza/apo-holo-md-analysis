import argparse
import pymol
from pymol import cmd


def main(args):    

    pymol.finish_launching(['pymol', '-c', '-Q'])
    
    cmd.reinitialize()    
    cmd.load(args.gro, 'trj')
    cmd.load_traj(args.trj, 'trj', stop=100)

    print(cmd.count_states('trj'))
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    # Required arguments
    parser.add_argument("--trj", type=str, help="Trajectory file")
    parser.add_argument("--gro", type=str, help="Gromac file")
    
    
    # Parse arguments
    args = parser.parse_args()
    
    # Call the main function with parsed arguments
    main(args)