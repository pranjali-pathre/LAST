#!/usr/bin/env python

"""Preliminary MD simulations
"""

import sys
import mdtraj as md
from amber import Amber


PDB = sys.argv[1]

# do 100ps nvt
sim_amber = Amber()
sim_amber.set_top_file(f"../inputs/{PDB}.prmtop")
sim_amber.set_cor_file(f"../inputs/{PDB}.inpcrd")
sim_amber.set_simulation_time(total_steps=50000, report_interval=1000)
sim_amber.set_up_simulation(sim_type='nvt', minimize_energy=True)
sim_amber.do_simulation(file_name=f'../trajs/{PDB}_100ps_nvt')

trajs = md.load(f'../trajs/{PDB}_100ps_nvt.dcd', top=f'../inputs/{PDB}.prmtop')
trajs[-1].save_amberrst7(f'../rst/{PDB}_100ps_nvt_rts')

# do 100ps npt
sim_amber.set_cor_file(f"../rst/{PDB}_100ps_nvt_rts")
sim_amber.set_simulation_time(total_steps=50000, report_interval=1000)
sim_amber.set_up_simulation(sim_type='npt', minimize_energy=True)
sim_amber.do_simulation(file_name=f'../trajs/{PDB}_100ps_npt')

trajs = md.load(f'../trajs/{PDB}_100ps_npt.dcd', top=f'../inputs/{PDB}.prmtop')
trajs[-1].save_amberrst7(f'../rst/{PDB}_100ps_npt_rts')
