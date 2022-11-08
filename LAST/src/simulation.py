#!/usr/bin/env python

"""Conduct MD simulations
"""

import sys
import glob
import MDAnalysis as mda
from MDAnalysis.analysis import align
from amber import Amber


PDB = sys.argv[1]
ITER_ROUND = int(sys.argv[2])

amber = Amber()
amber.set_top_file(f"../inputs/{PDB}.prmtop")

for i in range(1, 11):
    if ITER_ROUND == 0:
        amber.set_cor_file(f"../rst/{PDB}_100ps_npt_rts")
        amber.set_up_simulation(sim_type='nvt', minimize_energy=True)
    else:
        cor_file = f"../rst/{PDB}_r{ITER_ROUND}_rst." + "{0:02d}".format(i)
        amber.set_cor_file(cor_file)
        amber.set_up_simulation(sim_type='nvt', minimize_energy=False)

    amber.set_simulation_time(total_steps=50000, report_interval=500)
    dcd_file_name = f"../trajs/{PDB}_r{ITER_ROUND}_s" + "{0:02d}".format(i)
    amber.do_simulation(file_name=dcd_file_name)

    # prompt progress info
    print(f"round {ITER_ROUND} seed {i} finished")


# align trajs
trajs_dir = sorted(
    glob.glob(f"../trajs/{PDB}_r*.dcd"),
    key=lambda x: (int(x.split("_")[1][1:]), int(x.split("_")[2][1:-4])))

trajs = mda.Universe(f"../inputs/{PDB}.prmtop", trajs_dir)
ref = mda.Universe(f"../inputs/{PDB}.prmtop", trajs_dir[0])

align.AlignTraj(
    trajs,
    ref,
    select='protein and not name H*',
    filename=f'../trajs/{PDB}_aligned.dcd',
    match_atoms=True
).run()
