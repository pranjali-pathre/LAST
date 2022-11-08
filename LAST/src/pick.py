#!/usr/bin/env python
# pylint: disable=consider-using-f-string

"""Seed structure selection
"""

import sys
import os
import glob
import pickle
import numpy as np
from tensorflow import keras
from sklearn.preprocessing import MinMaxScaler
import mdtraj as md
import statsmodels.api as sm

PDB = sys.argv[1]
ITER_ROUND = int(sys.argv[2])

# heavy atom indices
trajs = md.load(
    f"../trajs/{PDB}_aligned.dcd", top=f"../inputs/{PDB}.prmtop"
)
trajs = trajs.atom_slice(
    trajs.topology.select_atom_indices("heavy")
)

frames = trajs.n_frames // (ITER_ROUND + 1)
coors = trajs.xyz
coors = coors.reshape(trajs.n_frames, -1)

scaler = MinMaxScaler()
coors_scaled = scaler.fit_transform(coors)

# load latents and decoder
with open(f"../results/{PDB}_r{ITER_ROUND}_latents.pkl", "rb") as latent_file:
    latents = pickle.load(latent_file)
latent_dim = latents.shape[1]
decoder = keras.models.load_model(
    f"../models/{PDB}_r{ITER_ROUND}_decoder"
)

decoded_structure = decoder(latents).numpy()
decoded_structure = scaler.inverse_transform(decoded_structure) * 10
real_structure = scaler.inverse_transform(coors_scaled) * 10


def rmsd(p_1, p_2, val=None):
    """RMSD calculation
    """
    if val is None:
        val = trajs.n_atoms
    p_1 = p_1.reshape(-1, 3)
    p_2 = p_2.reshape(-1, 3)

    assert p_1.shape  == p_2.shape
    assert val != 0
    return np.sqrt(np.sum(np.square(p_1 - p_2)) / val)

# nonparametric fit
dens_u = sm.nonparametric.KDEMultivariate(
    data=latents, var_type='c' * latent_dim, bw='normal_reference'
)
cdf = dens_u.cdf()

# sort based on CDF
idxs = np.arange(0, latents.shape[0])
cdf_idx = np.stack((cdf, idxs), axis=1)
cdf_idx = sorted(cdf_idx, key=lambda x: x[0])

# detect outliers
# traverse from beginning

if os.path.isfile(f"../results/{PDB}_structure.pkl"):
    outliers = []
    with open(f"../results/{PDB}_structure.pkl", "rb") as structure:
        outliers_structure = pickle.load(structure)
    CURR_IDX = 0
else:
    outliers = [int(cdf_idx[0][1])]
    outliers_structure = [real_structure[outliers[-1]]]
    CURR_IDX = 1

while CURR_IDX < len(cdf_idx) and len(outliers) < 5:
    current_structure = real_structure[int(cdf_idx[CURR_IDX][1])]
    current_latent = latents[int(cdf_idx[CURR_IDX][1])]
    # whether this point is NEAR within 1A
    NEAR = False

    for out_structure in outliers_structure:
        if rmsd(current_structure, out_structure) < 1:
            NEAR = True
            break

    if not NEAR:
        outliers.append(int(cdf_idx[CURR_IDX][1]))
        outliers_structure.append(real_structure[outliers[-1]])

    CURR_IDX += 1


# traverse from the other side
if os.path.isfile(f"../results/{PDB}_structure.pkl"):
    CURR_IDX = len(cdf_idx) - 1
else:
    outliers.append(int(cdf_idx[-1][1]))
    outliers_structure.append(real_structure[outliers[-1]])
    CURR_IDX = len(cdf_idx) - 2

while CURR_IDX >= 0 and len(outliers) < 10:
    current_structure = real_structure[int(cdf_idx[CURR_IDX][1])]
    current_latent = latents[int(cdf_idx[CURR_IDX][1])]
    NEAR = False

    for out_structure in outliers_structure:
        if rmsd(current_structure, out_structure) < 1:
            NEAR = True
            break

    if not NEAR:
        outliers.append(int(cdf_idx[CURR_IDX][1]))
        outliers_structure.append(real_structure[outliers[-1]])

    CURR_IDX -= 1


# delete data to spare some space
del trajs
del coors

for i, num in enumerate(outliers):
    round_num = num // frames
    seed_num = (num - round_num * frames) // 100 + 1
    STR1 = '../trajs/' + str(PDB) + '_r' + str(round_num) + '_'
    STR2 = 's{0:02d}.dcd'.format(seed_num)
    STR3 = "../inputs/" + str(PDB) + ".prmtop"
    trajs = md.load(STR1 + STR2, top=STR3)

    out_idx = num - round_num * frames - 100 * (seed_num - 1)
    STR4 = '../rst/' + str(PDB) + '_r' + str(ITER_ROUND + 1) + '_rst.'
    STR5 = '{0:02d}'.format(i + 1)
    trajs[out_idx].save_amberrst7(STR4 + STR5)

# save structure
with open(f"../results/{PDB}_structure.pkl", "wb") as structure_file:
    pickle.dump(outliers_structure,structure_file)

# save outliers index
with open(f"../results/{PDB}_out_{ITER_ROUND}.pkl", "wb") as outlier_file:
    pickle.dump(outliers,outlier_file)

# save rmsd
trajs_dir = sorted(
    glob.glob(f"../trajs/{PDB}_r*.dcd"),
    key=lambda x: (int(x.split("_")[1][1:]), int(x.split("_")[2][1:-4])))

# get the latest iteration
trajs_dir = trajs_dir[-10:]

ref = md.load(f"../inputs/{PDB}.PDB")
ref = ref.atom_slice(
    ref.topology.select_atom_indices("alpha")
)

trajs = md.load(trajs_dir, top=f"../inputs/{PDB}.prmtop")
trajs = trajs.atom_slice(
    trajs.topology.select_atom_indices("alpha")
)

rmsd = md.rmsd(trajs, ref) * 10
rmsd_mean = np.mean(rmsd)

# save mean rmsd
if os.path.isfile(f"../results/{PDB}_rmsd.pkl"):
    with open(f"../results/{PDB}_rmsd.pkl", "rb") as rmsd_file:
        rmsds = pickle.load(rmsd_file)
else:
    rmsds = []

rmsds.append(rmsd_mean)
with open(f"../results/{PDB}_rmsd.pkl", "wb") as rmsd_file:
    pickle.dump(rmsds, rmsd_file)
