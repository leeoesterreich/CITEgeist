"""Tests for VAE architecture."""
import sys
import os
import importlib.util

import torch
import pytest


def _import_vae():
    """Import vae module directly without triggering model __init__.py."""
    vae_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "model",
        "vae.py"
    )
    spec = importlib.util.spec_from_file_location("vae", vae_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_encoder_output_shape():
    """Encoder should output mu and logvar of correct shape."""
    vae_module = _import_vae()
    VAEEncoder = vae_module.VAEEncoder

    encoder = VAEEncoder(in_channels=3, latent_dim=128)
    x = torch.randn(4, 3, 96, 96)  # batch of 4 patches

    mu, logvar = encoder(x)

    assert mu.shape == (4, 128)
    assert logvar.shape == (4, 128)


def test_decoder_output_shape():
    """Decoder should reconstruct correct shape."""
    vae_module = _import_vae()
    VAEDecoder = vae_module.VAEDecoder

    decoder = VAEDecoder(out_channels=3, latent_dim=128)
    z = torch.randn(4, 128)

    x_recon = decoder(z)

    assert x_recon.shape == (4, 3, 96, 96)


def test_vae_forward_and_loss():
    """VAE forward pass and loss should work."""
    vae_module = _import_vae()
    VAE = vae_module.VAE

    vae = VAE(in_channels=3, latent_dim=128)
    x = torch.randn(4, 3, 96, 96)

    x_recon, mu, logvar = vae(x)

    assert x_recon.shape == x.shape
    assert mu.shape == (4, 128)

    loss, recon_loss, kl_loss = VAE.loss_function(x, x_recon, mu, logvar)

    assert loss.ndim == 0  # scalar
    assert loss > 0


def test_vae_encode():
    """Encode method should return deterministic latents."""
    vae_module = _import_vae()
    VAE = vae_module.VAE

    vae = VAE(in_channels=3, latent_dim=128)
    vae.eval()
    x = torch.randn(4, 3, 96, 96)

    z1 = vae.encode(x)
    z2 = vae.encode(x)

    # Should be deterministic (no sampling)
    assert torch.allclose(z1, z2)
