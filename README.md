# Analysing MD simulations of holo structures of apo-holo pairs

To play with the individual results, check the `pymol.ipynb`. It loads the MD data, starts PyMOL, does superposition, extracts pockets, colors them, computes RMSD between apo pocket and each snapshot of the simulation, ...

`run_analyses.py` analyses multiple MD runs by running `analyse.py` for each of them.