import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions.beta import Beta
from torchvision import datasets, transforms
from torch.utils.data import DataLoader

# Define the BetaDiffusionModel
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

# Define the BetaDiffusion process
class BetaDiffusion:
    def __init__(self, model, timesteps):
        self.model = model
        self.timesteps = timesteps
        self.beta_schedule = torch.linspace(0.01, 0.1, timesteps).to(device)
    
    def loss(self, x, alpha, beta):
        beta_dist = Beta(alpha, beta)
        log_prob = beta_dist.log_prob(x)
        return -log_prob.mean()

    def train_step(self, x, t):
        alpha, beta = self.model(x, t)
        loss = self.loss(x, alpha, beta)
        return loss

# Load CIFAR-10 dataset
transform = transforms.Compose([
    transforms.ToTensor(),
    transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))  # Normalize to [-1, 1]
])

train_dataset = datasets.CIFAR10(root="./data", train=True, transform=transform, download=True)
train_loader = DataLoader(train_dataset, batch_size=64, shuffle=True)

# Model parameters
in_channels = 3  # RGB images
hidden_channels = 64
timesteps = 100
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Initialize model and diffusion process
model = BetaDiffusionModel(in_channels, hidden_channels).to(device)
diffusion = BetaDiffusion(model, timesteps)
optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)

# Training loop
num_epochs = 1  # For demonstration purposes
model.train()
for epoch in range(num_epochs):
    for batch_idx, (x, _) in enumerate(train_loader):
        x = x.to(device)  # Move to GPU if available
        x = (x + 1) / 2  # Normalize to [0, 1] for Beta distribution
        t = torch.randint(0, timesteps, (x.size(0),), device=device)  # Random timesteps

        optimizer.zero_grad()
        loss = diffusion.train_step(x, t)
        loss.backward()
        optimizer.step()

        if batch_idx % 100 == 0:
            print(f"Epoch: {epoch+1}/{num_epochs}, Batch: {batch_idx}, Loss: {loss.item()}")
