import argparse
from pymol import cmd

PROJECTDIR=f"/storage/praha1/home/davidhoksza/projects/apoholo-trj-comp"

PATH_DATA = f"{PROJECTDIR}/data/G01/"
PATH_OUT = f"{PROJECTDIR}/out/"

class PATH:
    APO_GRO = f'{PATH_DATA}apo_2FJY/conf_wh20.gro'
    APO_TRJ = f'{PATH_DATA}apo_2FJY/traj_350ns_w_protein.xtc'
    APO_IX = f'{PATH_DATA}apo_2FJY/index.ndx'
    HOLO_GRO = f'{PATH_DATA}holo_2P70/conf_wh20.gro'
    HOLO_TRJ = f'{PATH_DATA}holo_2P70/traj_350ns_w_protein.xtc'
    HOLO_IX = f'{PATH_DATA}holo_2P70/index.ndx'

def main(args):    

    #pymol.finish_launching(['pymol', '-c', '-Q'])
    
    #cmd.reinitialize()        
    cmd.load(PATH.APO_GRO, 'trj')
    cmd.load_traj(PATH.APO_TRJ, 'trj')
    cnt_apo = cmd.count_states('trj')
    cmd.reinitialize()    
    cmd.load(PATH.HOLO_GRO, 'trj')
    cmd.load_traj(PATH.HOLO_TRJ, 'trj')
    cnt_holo = cmd.count_states('trj')

    window = int(args.window)

    f = open(args.output, 'w')
    i = 0
    while i < cnt_apo:        
        i1 = min(i+window, cnt_apo)
        j = 0
        while j < cnt_holo:
            j1 = min(j+window, cnt_holo)
            f.write(
                ('qsub -l select=1:ncpus=1:mem=6gb:scratch_local=10gb -l walltime=01:00:00 -v "'
                 f'apo_gro={PATH.APO_GRO},'
                 f'apo_trj={PATH.APO_TRJ},'
                 f'apo_ndx={PATH.APO_GRO},'
                 f'apo_ndx={PATH.APO_IX},'
                 f'apo_ix_range={i}-{i1},'
                 f'holo_gro={PATH.HOLO_GRO},'
                 f'holo_trj={PATH.HOLO_TRJ},'
                 f'holo_ndx={PATH.HOLO_IX},'
                 f'holo_ix_range={j}-{j1},'
                 f'blosum={PROJECTDIR}/BLOSUM62,'
                 f'output={PATH_OUT}rmsd_{i}-{i1}_{j}-{j1}'
                 '" job.sh\n')
                )
            j += window    
        i += window

    f.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    # Required arguments
    parser.add_argument("-w", "--window", type=int, help="Window size")
    parser.add_argument("-o", "--output", type=str, help="Output file")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Call the main function with parsed arguments
    main(args)

