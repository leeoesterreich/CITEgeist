"""Variational Autoencoder for nucleus patch representation learning.

This module implements a convolutional VAE for learning compressed
representations of nucleus image patches. The architecture is designed
for 96x96 pixel patches extracted around detected nuclei.

The VAE consists of:
- Encoder: 4 conv blocks (96->48->24->12->6) + adaptive pooling + FC
- Decoder: FC + 4 transposed conv blocks (6->12->24->48->96)
- Reparameterization trick for differentiable sampling

Usage:
    vae = VAE(in_channels=3, latent_dim=128)
    x_recon, mu, logvar = vae(x)
    loss, recon_loss, kl_loss = VAE.loss_function(x, x_recon, mu, logvar)

    # For inference (deterministic encoding)
    z = vae.encode(x)
"""
import torch
import torch.nn as nn
import torch.nn.functional as F


class VAEEncoder(nn.Module):
    """CNN encoder for nucleus patches.

    Encodes input patches to parameters (mu, logvar) of a latent
    Gaussian distribution using a series of strided convolutions.

    Architecture:
        Conv blocks with BatchNorm and ReLU, each halving spatial dims:
        96 -> 48 -> 24 -> 12 -> 6
        Then adaptive pooling and FC layers for mu/logvar.

    Attributes:
        latent_dim: Dimensionality of the latent space.
    """

    def __init__(self, in_channels: int = 3, latent_dim: int = 128):
        """Initialize encoder.

        Args:
            in_channels: Number of input channels (default 3 for RGB).
            latent_dim: Dimensionality of latent space.
        """
        super().__init__()
        self.latent_dim = latent_dim

        # Conv blocks: 96 -> 48 -> 24 -> 12 -> 6
        self.conv1 = nn.Sequential(
            nn.Conv2d(in_channels, 32, 3, stride=2, padding=1),
            nn.BatchNorm2d(32),
            nn.ReLU(inplace=True),
        )
        self.conv2 = nn.Sequential(
            nn.Conv2d(32, 64, 3, stride=2, padding=1),
            nn.BatchNorm2d(64),
            nn.ReLU(inplace=True),
        )
        self.conv3 = nn.Sequential(
            nn.Conv2d(64, 128, 3, stride=2, padding=1),
            nn.BatchNorm2d(128),
            nn.ReLU(inplace=True),
        )
        self.conv4 = nn.Sequential(
            nn.Conv2d(128, 256, 3, stride=2, padding=1),
            nn.BatchNorm2d(256),
            nn.ReLU(inplace=True),
        )

        # Adaptive pooling to handle variable input sizes
        self.pool = nn.AdaptiveAvgPool2d((1, 1))

        # FC layers for mu and logvar
        self.fc_mu = nn.Linear(256, latent_dim)
        self.fc_logvar = nn.Linear(256, latent_dim)

    def forward(self, x: torch.Tensor) -> tuple:
        """Encode patches to latent distribution parameters.

        Args:
            x: (B, C, H, W) input patches

        Returns:
            mu: (B, latent_dim) mean of latent distribution
            logvar: (B, latent_dim) log variance of latent distribution
        """
        x = self.conv1(x)
        x = self.conv2(x)
        x = self.conv3(x)
        x = self.conv4(x)
        x = self.pool(x)
        x = x.view(x.size(0), -1)

        mu = self.fc_mu(x)
        logvar = self.fc_logvar(x)

        return mu, logvar


class VAEDecoder(nn.Module):
    """CNN decoder for nucleus patches.

    Decodes latent vectors back to image patches using a series of
    transposed convolutions.

    Architecture:
        FC layer to 256*6*6, then transposed conv blocks:
        6 -> 12 -> 24 -> 48 -> 96
    """

    def __init__(self, out_channels: int = 3, latent_dim: int = 128):
        """Initialize decoder.

        Args:
            out_channels: Number of output channels (default 3 for RGB).
            latent_dim: Dimensionality of latent space.
        """
        super().__init__()

        self.fc = nn.Linear(latent_dim, 256 * 6 * 6)

        # Transposed conv blocks: 6 -> 12 -> 24 -> 48 -> 96
        self.deconv1 = nn.Sequential(
            nn.ConvTranspose2d(256, 128, 4, stride=2, padding=1),
            nn.BatchNorm2d(128),
            nn.ReLU(inplace=True),
        )
        self.deconv2 = nn.Sequential(
            nn.ConvTranspose2d(128, 64, 4, stride=2, padding=1),
            nn.BatchNorm2d(64),
            nn.ReLU(inplace=True),
        )
        self.deconv3 = nn.Sequential(
            nn.ConvTranspose2d(64, 32, 4, stride=2, padding=1),
            nn.BatchNorm2d(32),
            nn.ReLU(inplace=True),
        )
        self.deconv4 = nn.Sequential(
            nn.ConvTranspose2d(32, out_channels, 4, stride=2, padding=1),
        )

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        """Decode latent vector to patch.

        Args:
            z: (B, latent_dim) latent vectors

        Returns:
            x_recon: (B, C, 96, 96) reconstructed patches
        """
        x = self.fc(z)
        x = x.view(x.size(0), 256, 6, 6)
        x = self.deconv1(x)
        x = self.deconv2(x)
        x = self.deconv3(x)
        x = self.deconv4(x)
        return x


class VAE(nn.Module):
    """Full VAE for nucleus patch representation learning.

    Combines encoder and decoder with reparameterization trick for
    end-to-end training. Provides both stochastic forward pass
    (for training) and deterministic encoding (for inference).

    Attributes:
        encoder: VAEEncoder instance
        decoder: VAEDecoder instance
        latent_dim: Dimensionality of latent space
    """

    def __init__(self, in_channels: int = 3, latent_dim: int = 128):
        """Initialize VAE.

        Args:
            in_channels: Number of input/output channels.
            latent_dim: Dimensionality of latent space.
        """
        super().__init__()
        self.encoder = VAEEncoder(in_channels, latent_dim)
        self.decoder = VAEDecoder(in_channels, latent_dim)
        self.latent_dim = latent_dim

    def reparameterize(self, mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
        """Reparameterization trick: z = mu + std * epsilon.

        Enables backpropagation through the stochastic sampling step
        by expressing z as a deterministic function of mu, logvar,
        and random noise epsilon.

        Args:
            mu: (B, latent_dim) mean of latent distribution
            logvar: (B, latent_dim) log variance of latent distribution

        Returns:
            z: (B, latent_dim) sampled latent vectors
        """
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x: torch.Tensor) -> tuple:
        """Forward pass through VAE.

        Args:
            x: (B, C, H, W) input patches

        Returns:
            x_recon: (B, C, H, W) reconstructed patches
            mu: (B, latent_dim) latent means
            logvar: (B, latent_dim) latent log variances
        """
        mu, logvar = self.encoder(x)
        z = self.reparameterize(mu, logvar)
        x_recon = self.decoder(z)
        return x_recon, mu, logvar

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        """Encode to latent space (use mean, no sampling).

        For inference/embedding extraction, returns the mean of the
        latent distribution without stochastic sampling.

        Args:
            x: (B, C, H, W) input patches

        Returns:
            z: (B, latent_dim) latent vectors (deterministic)
        """
        mu, _ = self.encoder(x)
        return mu

    @staticmethod
    def loss_function(x: torch.Tensor, x_recon: torch.Tensor,
                      mu: torch.Tensor, logvar: torch.Tensor,
                      beta: float = 0.5) -> tuple:
        """Compute VAE loss.

        The loss is a weighted combination of reconstruction loss
        (MSE) and KL divergence to standard normal prior.

        Args:
            x: Original input
            x_recon: Reconstructed input
            mu: Latent means
            logvar: Latent log variances
            beta: KL weight (default 0.5 for beta-VAE behavior)

        Returns:
            loss: Total loss
            recon_loss: Reconstruction loss (MSE)
            kl_loss: KL divergence to N(0,1)
        """
        # Reconstruction loss (MSE)
        recon_loss = F.mse_loss(x_recon, x, reduction='mean')

        # KL divergence: -0.5 * sum(1 + log(sigma^2) - mu^2 - sigma^2)
        kl_loss = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())

        loss = recon_loss + beta * kl_loss

        return loss, recon_loss, kl_loss
