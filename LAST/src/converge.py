#!/usr/bin/env python

"""Convergence check
"""

import sys
import os
import pickle
import numpy as np


PDB = sys.argv[1]
PATIENCE = int(sys.argv[2])

if os.path.isfile(f"../results/{PDB}_rmsd.pkl"):
    with open(f"../results/{PDB}_rmsd.pkl", "rb") as rmsd_file:
        rmsds = pickle.load(rmsd_file)
        max_val_idx = np.argmax(rmsds)
        if len(rmsds) - 1 - max_val_idx >= PATIENCE:
            sys.exit("converged!")
