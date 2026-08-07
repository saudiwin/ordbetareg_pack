import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions.beta import Beta

class BetaDiffusionModel(nn.Module):
    def __init__(self, in_channels, hidden_channels):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Conv2d(in_channels, hidden_channels, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.Conv2d(hidden_channels, hidden_channels, kernel_size=3, padding=1),
            nn.ReLU(),
        )
        self.decoder = nn.Sequential(
            nn.Conv2d(hidden_channels, hidden_channels, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.Conv2d(hidden_channels, in_channels * 2, kernel_size=3, padding=1),
        )
    
    def forward(self, x, t):
        features = self.encoder(x)
        alpha_beta = self.decoder(features)  # Predict alpha and beta
        alpha, beta = torch.chunk(alpha_beta, 2, dim=1)
        alpha = F.softplus(alpha) + 1e-6  # Ensure positivity
        beta = F.softplus(beta) + 1e-6
        return alpha, beta

class BetaDiffusion:
    def __init__(self, model, timesteps):
        self.model = model
        self.timesteps = timesteps
        self.beta_schedule = torch.linspace(0.01, 0.1, timesteps)
    
    def forward_diffusion(self, x):
        noisy_images = []
        for t in range(self.timesteps):
            beta_t = self.beta_schedule[t]
            alpha_t = 1 - beta_t
            x_noisy = Beta(alpha_t, beta_t).sample(x.shape)
            noisy_images.append(x_noisy)
        return noisy_images

    def loss(self, x, alpha, beta):
        beta_dist = Beta(alpha, beta)
        log_prob = beta_dist.log_prob(x)
        return -log_prob.mean()

    def train_step(self, x, t):
        alpha, beta = self.model(x, t)
        loss = self.loss(x, alpha, beta)
        return loss

# Example Usage
in_channels = 3  # RGB images
hidden_channels = 64
timesteps = 100

model = BetaDiffusionModel(in_channels, hidden_channels)
diffusion = BetaDiffusion(model, timesteps)

# Dummy data
x = torch.rand(8, 3, 32, 32)  # Batch of 8 images (32x32)

# Training step
t = torch.randint(0, timesteps, (x.size(0),))  # Random timesteps
loss = diffusion.train_step(x, t)

print("Loss:", loss.item())
