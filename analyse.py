import structure
import sys
import pymol
import numpy as np
import json 
import argparse
import logging
import typing
import os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px

from timeit import default_timer as timer

def init_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
    
COLORS = {
    'apo': 'blue',
    'holo': 'green',    
    'trj': 'yellow',    
    'lig': 'white'
}
POCKET_THRESHOLD = 6

class PATH:
    GRO = ''
    TRJ = ''
    APO = ''
    HOLO = ''
CHAIN_APO = ""
CHAIN_HOLO = ""
LIGAND = ""

MAX_DIST_THRESHOLD = 16

# different types of distances between a pocket and pockets in trajectory
class PocketDistances:
    rmsd_all_atoms:np.ndarray[typing.Any, float] = None # RMSD distance based on all atoms
    rmsd_ca_com:np.ndarray[typing.Any, float] = None # RMSD of pockets in terms of C-alpa and center of mass points
    l2_apo:np.ndarray[typing.Any, np.ndarray] = None #Euclidean distances of all atoms - one array per pocket -> 2D array
    l2_apo_max:np.ndarray[typing.Any, float] = None # For each pocket stores max distance
    sasa_trj_relative_to_start:np.ndarray[typing.Any, float] = None
    rmsd_trj_relative_to_start_all_atoms:np.ndarray[typing.Any, float] = None
    
def align():
    pymol.cmd.intra_fit('trj')
    pymol.cmd.align(mobile='apo', target='trj', target_state=1)
    pymol.cmd.align(mobile='holo', target='trj', target_state=1)

def create_ligand():
    pymol.cmd.create('lig', f'holo and chain {CHAIN_HOLO} and resn {LIGAND}')
    pymol.cmd.remove(f'holo and resn {LIGAND}')

def color():
    for sel in ['apo', 'holo', 'trj', 'lig']:
        pymol.cmd.color(COLORS[sel], sel)

def create_pockets(pocket_threshold):
    def_holo = f'byres(holo and polymer.protein within {pocket_threshold} of resn {LIGAND})'
    pymol.cmd.create('holo-pocket', def_holo)

    holo_pocket_res_sel = ' or '.join([f"(resi {r.resi} and resn {r.resn})" for r in structure.Structure(pymol.cmd.get_model('holo-pocket'), def_holo).get_residues()])

    def_trj = f"trj and ({holo_pocket_res_sel})"
    def_apo = f"apo and ({holo_pocket_res_sel})"
    pymol.cmd.create("trj-pocket", def_trj)
    pymol.cmd.create("apo-pocket", def_apo)

    pymol.cmd.show_as('sticks', 'holo-pocket')
    pymol.cmd.show_as('sticks', 'trj-pocket')
    pymol.cmd.show_as('sticks', 'apo-pocket')

    return def_apo, def_holo, def_trj

def get_pockets(def_apo:str, def_holo:str, def_trj:str) -> tuple[structure.Structure, structure.Structure, list[structure.Structure]]:
    pocket_apo = structure.Structure(pymol.cmd.get_model('apo-pocket'), def_apo)
    pocket_holo = structure.Structure(pymol.cmd.get_model('holo-pocket'),def_holo)

    cnt_states = pymol.cmd.count_states('trj')
    pocket_trj = []
    for i in range(cnt_states):
        pocket_trj.append(structure.Structure(pymol.cmd.get_model("trj-pocket", i+1), def_trj, i+1))

    #assert(len(pocket_apo.atoms) == len(pocket_holo.atoms))
    assert(len(pocket_apo.atoms) == len(pocket_trj[0].atoms))

    return pocket_apo, pocket_holo, pocket_trj

def get_coords_ca_com(atoms:list[structure.Atom])->np.ndarray:
    # Retrieves coordinates of the CA residue and center of mass
    ca = None
    for a in atoms:
        if a.name == 'CA':
            ca = a.coord
    com = np.mean(structure.get_coords(atoms), axis=0)
    return np.array([ca, com])

def compute_pocket_distances(pocket_apo:structure.Structure, pocket_trj:list[structure.Structure]):
    
    dists = PocketDistances()
    
    logging.info("\tComputing RMSDs to APO")
    dists.rmsd_all_atoms = [pocket_apo.rmsd(pt) for pt in pocket_trj]
    dists.rmsd_ca_com = [pocket_apo.rmsd(pt, get_coords_from_residue=get_coords_ca_com) for pt in pocket_trj]
    
    logging.info("\tComputing L2 to APO")
    dists.l2_apo = [pocket_apo.get_distances(pt) for pt in pocket_trj]    
    dists.l2_apo_max = [np.max(dists) for dists in dists.l2_apo]
    dists.rmsd_trj_relative_to_start_all_atoms = [pocket_trj[0].rmsd(pt) for pt in pocket_trj]

    logging.info("\tComputing SASA to trj1")
    sasa1 = pocket_trj[0].get_sasa(pymol.cmd)
    dists.sasa_trj_relative_to_start = [sasa1 - p.get_sasa(pymol.cmd) for p in pocket_trj]    
    dists.sasa_pct_trj_relative_to_start = [(1-sasa1/p.get_sasa(pymol.cmd))*100 for p in pocket_trj ]

    return dists

def compute_sasa(pocket_apo, pocket_holo, pocket_trj):    
    logging.info("\tSetting dot_solvent")
    pymol.cmd.set('dot_solvent', 1)
    logging.info("\tComputing SASA")
    pocket_apo.get_sasa(pymol.cmd)
    pocket_holo.get_sasa(pymol.cmd)
    for pocket in pocket_trj:
        pocket.get_sasa(pymol.cmd)

def dist_outlier(distances:np.ndarray, threshold:float=4)->bool:
    for d in distances:
        if d > threshold:
            return True        
    return False

def plot_distances_table(dists: PocketDistances, fig, row, col):
    sort_ix_all = np.argsort(dists.rmsd_all_atoms)
    sort_ix_ca_com = np.argsort(dists.rmsd_ca_com)
    fig.add_trace(
        go.Table(
            header=dict(
                values=["Order", "Index - RMSD All", "Distance - RMSD All", "Index - RMSD CA-CoM", "Distance - RMSD CA-CoM"],
                #font=dict(size=10),
                align="left"
            ),
            cells=dict(
                values=[list(range(1,len(sort_ix_all)+1)), sort_ix_all, [dists.rmsd_all_atoms[ix] for ix in sort_ix_all], sort_ix_ca_com, [dists.rmsd_ca_com[ix] for ix in sort_ix_ca_com]],
                # fill=dict(color=['white', 'white', 'lightgray', 'lightgray']),
                align = "left")
        ),
        row=row, col=col
    )
    return row+1

def plot_min_distances_table(dists: PocketDistances, fig, row, col):
    maxs = dists.l2_apo_max
    ixs_order_maxs = np.argsort(maxs)
    data = [maxs[ix] for ix in ixs_order_maxs]
    fig.add_trace(
        go.Table(
            header=dict(
                values=["Order", "Index", "Distance"],
                #font=dict(size=10),
                align="left"
            ),
            cells=dict(
                values=[list(range(1,len(data)+1)), ixs_order_maxs, data],                
                align = "left")
        ),
        row=row, col=col
    )
    return row+1

def plot_distance_relative_to_trj1(dists: PocketDistances, fig, ix_row, ix_col):
    trace1 = go.Scatter(y=dists.rmsd_trj_relative_to_start_all_atoms, mode='lines', name='RMSD (all atoms)')
    

    fig.add_trace(
        trace1,
        row=ix_row, col=ix_col
    )
    ix_row+=1
    
    fig.add_trace(
        go.Scatter(y=dists.sasa_trj_relative_to_start, mode='lines', name='SASA'),
        secondary_y=False,
        row=ix_row, col=ix_col
    )
    fig.add_trace(
        go.Scatter(y=dists.sasa_pct_trj_relative_to_start, mode='lines', name='SASA %'),
        secondary_y=True,
        row=ix_row, col=ix_col
    )
    fig.update_yaxes(title_text="Absolute difference", row=ix_row, col=ix_col, secondary_y=False)
    fig.update_yaxes(title_text="% difference", row=ix_row, col=ix_col, secondary_y=True)

    return ix_row + 1

def plot_pocket_descriptors(dists: PocketDistances, pocket_trj:list[structure.Structure], fig, row, col):
    sasa:list[float] = [p.get_sasa(pymol.cmd) for p in pocket_trj]
    fig.add_trace(
        go.Table(
            header=dict(
                values=["Index", 
                        "RMSD (all TRJ pocket atoms) to APO",
                        "RMSD (CA-CoM TRJ pocket atoms) to APO",
                        "Maximum L2 distance of TRJ pocket atoms to APO atoms",
                        "RMSD (all pocket atoms) of TRJ relative to TRJ1",
                        "SASA",
                        "SASA (difference) of TRJ relative to TRJ1"
                        ],
                #font=dict(size=10),
                align="left"
            ),
            cells=dict(
                values=[list(range(1,len(dists.rmsd_all_atoms)+1)), dists.rmsd_all_atoms, dists.rmsd_ca_com, dists.l2_apo_max, dists.rmsd_trj_relative_to_start_all_atoms, sasa, dists.sasa_trj_relative_to_start],
                align = "left")
        ),
        row=row, col=col
    )
    return row+1


def plot_distances(dists: PocketDistances, pocket_trj:list[structure.Structure], in_dir:str, out_dir:str):

    titles = ["Trajectories descriptors", "RMSDs", "Histogram", "Ordered distances", "RMSD realative to the first trajectory", "SASA realative to the first trajectory", f"Minimum distance under which all trajectory pocket atoms are from apo pocket"]    
    cnt_rows = len(titles)
    ix_row = 1
    fig = make_subplots(
        rows=cnt_rows, cols=1,
        # shared_xaxes=True,
        vertical_spacing=1.0 / float(cnt_rows - 1)/10.0,
        specs=[
            [{"type": "table"}],
            [{"type": "scatter"}],
            [{"type": "bar"}],            
            [{"type": "table"}],
            [{"type": "scatter"}],
            [{"type": "xy", "secondary_y": True}],
            [{"type": "table"}],            
            ],
        subplot_titles=titles
    )

    ix_row = plot_pocket_descriptors(dists, pocket_trj, fig, ix_row, 1)
    
    trace1 = go.Scatter(y=dists.rmsd_all_atoms, mode='lines', name='RMSD (all atoms)')
    trace2 = go.Scatter(y=dists.rmsd_ca_com, mode='lines', name='RMSD (CA atoms and center of mass)')    

    fig.add_trace(
        trace1,
        row=ix_row, col=1
    )
    fig.add_trace(
        trace2,
        row=ix_row, col=1
    )
    ix_row+=1

    fig.add_trace(
        go.Histogram(x=dists.rmsd_all_atoms, texttemplate="%{y}", name="Hist RMSD (all atoms)"),
        row=ix_row, col=1
    )
    fig.add_trace(
        go.Histogram(x=dists.rmsd_ca_com, texttemplate="%{y}", name="Hist RMSD (CA+CoM)"),
        row=ix_row, col=1
    )
    # fig.update_layout(barmode='overlay', row=ix_row, col=1)    
    fig.update_traces(opacity=0.75, row=ix_row, col=1)

    ix_row+=1
    
    ix_row = plot_distances_table(dists, fig, ix_row, 1)
    ix_row = plot_distance_relative_to_trj1(dists, fig, ix_row, 1)
    ix_row = plot_min_distances_table(dists, fig, ix_row, 1)
    
    

    fig.update_xaxes(title_text="Trajectory", row=1, col=1)
    fig.update_yaxes(title_text="RMSD", row=1, col=1)

    fig.update_layout(         
        title_text=f"Report for {in_dir}",
        title_font_size=30,  
        height=cnt_rows*600              
    )

    fig.write_html(f"{out_dir}/plot.html")

    #sns.lineplot(data=pd.DataFrame({'rmsd': dists.rmsd_all_atoms, 'ix': list(range(len(dists.rmsd_all_atoms)))}), x="ix", y="rmsd").set(title='RMSD (all atoms)')
    #sns.lineplot(data=pd.DataFrame({'rmsd': dists.rmsd_ca_com, 'ix': list(range(len(dists.rmsd_ca_com)))}), x="ix", y="rmsd").set(title='RMSD (CA atoms and center of mass)')
    #plt.savefig(f"{out_dir}/dists_rmsd_all.svg", format='svg')

    # sns.lineplot(data=df, dashes=False)
    # plt.savefig(f"{out_dir}/dists.svg", format='svg')


def check_structures():

    apo = structure.Structure(pymol.cmd.get_model('apo'), 'apo')
    # holo = structure.Structure(pymol.cmd.get_model('holo'),'holo')
    trj = structure.Structure(pymol.cmd.get_model('trj', 1), 'trj')

    # apo.sort_atoms()
    # holo.sort_atoms()
    # trj.sort_atoms()

    a = [(a.resn, a.resi, a.name) for a in apo.atoms]
    # h = [(a.resn, a.resi, a.name) for a in holo.atoms]
    t = [(a.resn, a.resi, a.name) for a in trj.atoms]

    if len(a) != len(t):        
        logging.error(f'Atoms in apo not in trj: {sorted(set(a)-set(t), key=lambda x: (x[1], x[0], x[2]) )}')
        logging.error(f'Atoms in trj not in apo: {sorted(set(t)-set(a), key=lambda x: (x[1], x[0], x[2]) )}')
        raise Exception('Different number of atoms in apo vs trj')
    

def trim_apo():
    apo = structure.Structure(pymol.cmd.get_model('apo'), 'apo')    
    trj = structure.Structure(pymol.cmd.get_model('trj', 1), 'trj')

    a = [(a.resn, a.resi, a.name) for a in apo.atoms]
    t = [(a.resn, a.resi, a.name) for a in trj.atoms]

    if len(set(a)-set(t)):
        pymol.cmd.remove(f"apo and ({' or '.join([f'(resn {x[0]} and resi {x[1]} and name {x[2]})' for x in set(set(a)-set(t))])})")
    if len(set(t)-set(a)):
        pymol.cmd.remove(f"trj and ({' or '.join([f'(resn {x[0]} and resi {x[1]} and name {x[2]})' for x in set(set(t)-set(a))])})")


def create_pymol_script(fn):
    NotImplemented
    # with f as open(fn, 'w'):


def analyse(input_dir:str, output_dir:str, pocket_threshold:int,):    
    global CHAIN_APO, CHAIN_HOLO, LIGAND

    fn_settings = f'{input_dir}/vis.json'
    try:
        with open(fn_settings) as f:
        
            settings = json.load(f)
            
            PATH.GRO = f'{input_dir}/{settings["path"]["gro"]}'
            PATH.TRJ = f'{input_dir}/{settings["path"]["trj"]}'
            PATH.APO = f'{input_dir}/{settings["path"]["apo"]}'
            PATH.HOLO = f'{input_dir}/{settings["path"]["holo"]}'

            CHAIN_APO = settings["chain_apo"]
            CHAIN_HOLO = settings["chain_holo"]
            LIGAND = settings["ligand"]

    except FileNotFoundError:
        logging.error(f'File {fn_settings} not found')
        exit()
    
    except IOError:
        logging.error(f'Error opening the file {fn_settings}')
        exit()

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    _stdouterr = sys.stdout, sys.stderr #PyMOL captures sys.stdout and sys.stderr, to control it with it's own feedback mechanism. To prevent that, save and restore both streams
    pymol.finish_launching(['pymol', '-c', '-s', f'{output_dir}/pymol.log'])
    sys.stdout, sys.stderr = _stdouterr

    # pymol.cmd.pwd()
    pymol.cmd.load(PATH.GRO, 'trj')
    # pymol.cmd.load_traj(PATH.TRJ, 'trj', max=10)
    pymol.cmd.load_traj(PATH.TRJ, 'trj')
    pymol.cmd.load(PATH.APO, 'apo')
    pymol.cmd.load(PATH.HOLO, 'holo')
    pymol.cmd.remove(f'holo and not(holo and chain {CHAIN_HOLO} and (polymer or resn {LIGAND}))')
    pymol.cmd.remove(f'apo and not(apo and chain {CHAIN_APO})')
    #pymol.cmd.remove(f'apo and (solvent or not polymer)')    
    pymol.cmd.remove(f'apo and not polymer')    
    pymol.cmd.remove("apo and not alt ''+A") #Removing alterantive coordinates based on https://pymol.org/dokuwiki/doku.php?id=concept:alt
    pymol.cmd.alter("apo", "alt=''")
    pymol.cmd.remove(f'holo and solvent')
    pymol.cmd.remove(f'trj and hydrogens')
    pymol.cmd.remove(f'apo and hydrogens')
    pymol.cmd.remove(f'holo and hydrogens')

    trim_apo()
    check_structures()
    
    align()
    create_ligand()
    color()
    def_apo, def_holo, def_trj = create_pockets(pocket_threshold=pocket_threshold)    

    pymol.cmd.save(f'{output_dir}/session.pse')
    create_pymol_script(f'{output_dir}/script.py')

    logging.info("Structures loaded, pockets extracted, session file stored")
    
    logging.info("Converting pockets into internal representation")
    pocket_apo, pocket_holo, pocket_trj = get_pockets(def_apo, def_holo, def_trj)
    logging.info("Computing SASA")
    compute_sasa(pocket_apo, pocket_holo, pocket_trj)
    logging.info("Computing pocket distances")
    dists = compute_pocket_distances(pocket_apo=pocket_apo, pocket_trj=pocket_trj)

    logging.info("Plotting distances")
    plot_distances(dists, pocket_trj, input_dir, output_dir)
    

if __name__ == '__main__':

    init_logging()

    parser = argparse.ArgumentParser(description='Analysis of gromacs results')
    
    parser.add_argument('-i', help='Path to the MD results')
    parser.add_argument('-o', help='Output directory')
    parser.add_argument('-t', default=POCKET_THRESHOLD, nargs='?', help='Pocket definition threshold in Angstroms')
    
    args = parser.parse_args()

    start = timer()
    analyse(args.i, args.o, args.t)
    end = timer()
    logging.info(end - start)



