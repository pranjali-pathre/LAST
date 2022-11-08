#!/usr/bin/env python

"""Variational autoencoder model Testing
"""
# from tensorflow import keras
from keras import backend as K
from keras.layers import Input, Dense, Lambda
from vae import vae_encoder
from vae import sampling

def test_vae_encoder():
    """unit test for VAE Encoder
    """
    input_dim = (256,256)
    neuron_nums = [512, 128, 32]
    latent_dim = 2

    _, z_mean, z_log_var, encoder_input = vae_encoder(input_dim, neuron_nums, latent_dim)

    assert K.int_shape(z_mean) == (None, 256, 2)
    assert K.int_shape(z_log_var) == (None, 256, 2)
    assert K.int_shape(encoder_input) == (None, 256, 256)

def test_sampling():
    """unit test for sampling function
    """
    input_dim = (256,256)
    neuron_nums = [512, 128, 32]
    latent_dim = 2

    layer = Input(shape=input_dim)

    for neuron_num in neuron_nums:
        layer = Dense(neuron_num, activation='relu')
        layer = layer(layer)

    z_mean = Dense(latent_dim)
    layer = layer(layer)
    z_log_var = Dense(latent_dim)
    z_log_var = z_log_var(layer)

    sampled_z = Lambda(sampling)([z_mean, z_log_var, latent_dim])

    assert K.int_shape(sampled_z) == (None, 256, 256)

test_vae_encoder()
test_sampling()
