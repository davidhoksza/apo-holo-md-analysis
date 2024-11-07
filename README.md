# Analysing MD simulations of holo structures of apo-holo pairs

To play with the individual results, check the `pymol.ipynb`. It loads the MD data, starts PyMOL, does superposition, extracts pockets, colors them, computes RMSD between apo pocket and each snapshot of the simulation, ...

`run_analyses.py` analyses multiple MD runs by running `analyse.py` for each of them.

# Analysing apo-holo simulations

The `apo-holo-trj-comp` contains code to obtain all-to-all pocket RMSD distances from MD simulations of apo and holo version of a protein. The pocket definition is expected to be the same in both apo and holo structures (simulations).

 - `apo-holo-trj-comp.ipynb` allows for interactive exploration of the two simulations and RMSD computations - good for prototypng and debugging.
 - `apo-holo-trj-comp.py` takes two structures (defined by gromacs, trajectory and index files) and two ranges and computes all-to-all pocket RMSD distances of snapshots in the defined ranges
 - `job.sh` - runs the previous pythons script in [Metacentrum](https://metavo.metacentrum.cz/) ev
 - `prepare_jobs.py` - creates a shell script with with index ranges covering all the snapshots where the window size is passed as an argument. You can modify the parameters of qsub by editing the script (for example, if window size is selected too large [like 5000, for instance], then the job might require more RAM)

 Test data is available [here](https://cunicz-my.sharepoint.com/:f:/g/personal/51137390_cuni_cz/Em9WJpFIlJ9Fk0EWAL1CC5MBMiP1dHf_0nRzO4-QwogpSQ?e=te8TJK).

The following scripts are here just as an example to c&p when prototyping.
```
python .\apo-holo-trj-comp.py `
    --apo-gro c:/git/apo-holo-md-analysis/data/apo-holo-trj-comp/G01/apo_2FJY/conf_wh20.gro `
    --apo-trj c:/git/apo-holo-md-analysis/data/apo-holo-trj-comp/G01/apo_2FJY/traj_350ns_w_protein.xtc `
    --apo-ndx c:/git/apo-holo-md-analysis/data/apo-holo-trj-comp/G01/apo_2FJY/index.ndx `
    --apo-ix-range 0-20 `
    --holo-gro c:/git/apo-holo-md-analysis/data/apo-holo-trj-comp/G01/holo_2P70/conf_wh20.gro `
    --holo-trj c:/git/apo-holo-md-analysis/data/apo-holo-trj-comp/G01/holo_2P70/traj_350ns_w_protein.xtc `
    --holo-ndx c:/git/apo-holo-md-analysis/data/apo-holo-trj-comp/G01/holo_2P70/index.ndx `
    --holo-ix-range 0-20 `
    --output rmsd_0_20.npy
```

```
python apo-holo-trj-comp.py \
    --apo-gro data/G01/apo_2FJY/conf_wh20.gro \
    --apo-trj data/G01/apo_2FJY/traj_350ns_w_protein.xtc \
    --apo-ndx data/G01/apo_2FJY/index.ndx \
    --apo-ix-range 0-20 \
    --holo-gro data/G01/holo_2P70/conf_wh20.gro \
    --holo-trj data/G01/holo_2P70/traj_350ns_w_protein.xtc \
    --holo-ndx data/G01/holo_2P70/index.ndx \
    --holo-ix-range 0-20 \
    --output rmsd_0_20.npy
```