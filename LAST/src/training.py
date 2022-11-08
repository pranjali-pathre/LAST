#!/usr/bin/env python

"""Train VAE model
"""

import sys
import pickle
import mdtraj as md
from sklearn.preprocessing import MinMaxScaler
from vae import build_vae

PDB = sys.argv[1]
ITER_ROUND = sys.argv[2]
ITER_ROUND = int(ITER_ROUND)
STR1 = "../trajs/" + str(PDB) + "_aligned.dcd"
STR2 = "../inputs/" + str(PDB) + ".prmtop"
# heavy atom indices
trajs = md.load(STR1, top=STR2)
trajs = trajs.atom_slice(
    trajs.topology.select_atom_indices("heavy")
)

frames = trajs.n_frames
frames = frames // (ITER_ROUND + 1)
coors = trajs
coors = coors.xyz
coors = coors.reshape(trajs.n_frames, -1)

# scale
scaler = MinMaxScaler()
coors_scaled = scaler.fit_transform(coors)

encoder, decoder, vae = build_vae(
    input_dim=coors_scaled.shape[1:],
    encoder_neuron_nums=[512, 128, 32],
    latent_dim=2
)

HISTORY = vae.fit(
    x=coors_scaled, y=coors_scaled,
    shuffle=True, epochs=400, batch_size=16
)

# get latents
latents = encoder.predict(coors_scaled)[0]

# save latents, models
with open(f"../results/{PDB}_r{ITER_ROUND}_latents.pkl", "wb") as latent_file:
    pickle.dump(latents, latent_file)
encoder.save(f"../models/{PDB}_r{ITER_ROUND}_encoder")
decoder.save(f"../models/{PDB}_r{ITER_ROUND}_decoder")
