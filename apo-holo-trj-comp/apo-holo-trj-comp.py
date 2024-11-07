import argparse
import structure
from pymol import cmd
import numpy as np
import re


def main(args):
    
    class PATH:
        APO_GRO = args.apo_gro
        APO_TRJ = args.apo_trj
        APO_IX = args.apo_ndx
        HOLO_GRO = args.holo_gro
        HOLO_TRJ = args.holo_trj
        HOLO_IX = args.holo_ndx

    range_apo = [int(x) for x in args.apo_ix_range.split('-')]
    range_holo = [int(x) for x in args.holo_ix_range.split('-')]

    cmd.reinitialize()
    cmd.load(PATH.APO_GRO, 'trj-apo')
    cmd.load_traj(PATH.APO_TRJ, 'trj-apo')
    cmd.load(PATH.HOLO_GRO, 'trj-holo')
    cmd.load_traj(PATH.HOLO_TRJ, 'trj-holo')

    cmd.align(mobile='trj-holo', target='trj-apo', mobile_state=1, target_state=1, cycles=10, matrix=args.blosum) #the matrix parameter is here because PyMOL module on MetaCentrum has wrong mappings and can't loacate the blossum matrix; in this way, it uses the one provided in the local directory
    cmd.intra_fit('trj-apo')
    cmd.intra_fit('trj-holo')

    # Extract pocket indexes
    f = open(PATH.APO_IX, 'r').read()
    m = re.findall(r'r_[0-9]+_*', f)
    pocket_ixs_apo = [r.split('_')[1] for r in m]

    f = open(PATH.HOLO_IX, 'r').read()
    m = re.findall(r'r_[0-9]+_*', f)
    pocket_ixs_holo = [r.split('_')[1] for r in m]

    assert set(pocket_ixs_apo) == set(pocket_ixs_holo)

    
    #Create pocket selections
    pocket_res_sel = ' or '.join([f"(resi {r})" for r in pocket_ixs_apo])
    def_apo = f"trj-apo and ({pocket_res_sel})"
    def_holo = f"trj-holo and ({pocket_res_sel})"
    cmd.create("pocket-apo", def_apo)
    cmd.create("pocket-holo", def_holo)
    
    pocket_trj_apo = []
    for i in range(range_apo[0], range_apo[1]):
        pocket_trj_apo.append(structure.Structure(cmd.get_model("pocket-apo", i+1), def_apo, i+1))
    pocket_trj_holo = []
    for i in range(range_holo[0], range_holo[1]):
        pocket_trj_holo.append(structure.Structure(cmd.get_model("pocket-holo", i+1), def_holo, i+1))

    #cnt_states_apo = cmd.count_states('trj-apo')
    #cnt_states_holo = cmd.count_states('trj-holo')

    rmsd_apo_holo = np.empty(shape=(range_apo[1]-range_apo[0], range_holo[1]-range_holo[0]))

    for i in range(len(pocket_trj_apo)):
        #i_real = range_apo[0] + i
        for j in range(i, len(pocket_trj_holo)):
            #j_real = range_holo[0] + j

            aux = pocket_trj_apo[i].rmsd(pocket_trj_holo[j])

            # rmsd_apo_holo[i_real,j_real] = aux
            # rmsd_apo_holo[j_real,i_real] = aux
            rmsd_apo_holo[i,j] = aux
            rmsd_apo_holo[j,i] = aux

    np.save(args.output, rmsd_apo_holo)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    # Required arguments
    parser.add_argument("--apo-trj", type=str, help="Trajectory file of apo structure")
    parser.add_argument("--apo-gro", type=str, help="Gromac file of apo structure")
    parser.add_argument("--apo-ndx", type=str, help="Index file of apo structure")
    parser.add_argument("--apo-ix-range", type=str, help="Range of frames to process. The format needs to be: ix_from-ix_to")

    parser.add_argument("--holo-trj", type=str, help="Trajectory file of holo structure")
    parser.add_argument("--holo-gro", type=str, help="Gromac file of holo structure")
    parser.add_argument("--holo-ndx", type=str, help="Index file of holo structure")
    parser.add_argument("--holo-ix-range", type=str, help="Range of frames to process. The format needs to be: ix_from-ix_to")

    parser.add_argument("--blosum", type=str, help="Path to blossum matrix to be used when superposing")

    parser.add_argument("--output", type=str, help="Path to the output npy file containing a matrix of dimension (apo-ix-range[1]-apo-ix-range[0], holo-ix-range[1]-holo-ix-range[0])")
    
    
    # Parse arguments
    args = parser.parse_args()
    
    # Call the main function with parsed arguments
    main(args)