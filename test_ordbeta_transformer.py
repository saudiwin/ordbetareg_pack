import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.datasets import fetch_california_housing
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_squared_error

# Load and preprocess dataset
def load_data():
    data = fetch_california_housing()
    X, y = data['data'], data['target']
    # Normalize X and y to [0, 1]
    scaler_X = MinMaxScaler()
    scaler_y = MinMaxScaler()
    X = scaler_X.fit_transform(X)
    y = scaler_y.fit_transform(y.reshape(-1, 1)).flatten()
    return train_test_split(X, y, test_size=0.2, random_state=42)

def dordbeta(x, mu, phi, cutpoints=torch.tensor([-1.0, 1.0]), log=False):
    # Clamp inputs to avoid invalid values
    x = torch.clamp(x, min=1e-6, max=1 - 1e-6)
    mu = torch.clamp(mu, min=1e-6, max=1 - 1e-6)
    phi = torch.clamp(phi, min=1e-6, max=1e6)

    # Transform mu to logit space
    mu_ql = torch.logit(mu)

    # Ensure compatibility of cutpoints
    if mu_ql.numel() == 1:
        mu_ql = mu_ql.expand(x.shape)

    # Compute probabilities for three regions
    low = 1 - torch.sigmoid(mu_ql - cutpoints[0])
    middle = torch.sigmoid(mu_ql - cutpoints[0]) - torch.sigmoid(mu_ql - cutpoints[1])
    high = torch.sigmoid(mu_ql - cutpoints[1])

    # Beta distribution for the middle region
    alpha = mu * phi
    beta = (1 - mu) * phi
    beta_dist = torch.distributions.Beta(alpha, beta)

    # Compute probabilities safely
    try:
        out_beta = torch.exp(beta_dist.log_prob(x)) if not log else beta_dist.log_prob(x)
    except ValueError:
        print("Invalid alpha/beta values:", alpha, beta)
        out_beta = torch.zeros_like(x)

    # Assign probabilities based on region
    outcomes = torch.zeros_like(x)
    outcomes[x == 0] = low[x == 0]
    outcomes[x == 1] = high[x == 1]
    outcomes[(x > 0) & (x < 1)] = middle[(x > 0) & (x < 1)] * out_beta[(x > 0) & (x < 1)]

    return torch.log(outcomes) if log else outcomes



# Define Transformer model
class TransformerRegressor(nn.Module):
    def __init__(self, input_dim, hidden_dim=64, nhead=4, num_layers=3):
        super(TransformerRegressor, self).__init__()
        self.input_proj = nn.Linear(input_dim, hidden_dim)
        self.pos_encoder = nn.Linear(hidden_dim, hidden_dim)
        self.transformer = nn.Transformer(hidden_dim, nhead, num_layers)
        self.fc_out = nn.Linear(hidden_dim, 1)

    def forward(self, x):
        x = self.input_proj(x).unsqueeze(1)
        x = self.pos_encoder(x)
        x = self.transformer(x, x)
        return torch.sigmoid(self.fc_out(x[:, -1, :]))  # Ensure outputs in [0, 1]

# Training function
def train_model(model, criterion, optimizer, X_train, y_train, epochs=50, use_dordbeta=False):
    for epoch in range(epochs):
        optimizer.zero_grad()
        inputs = torch.tensor(X_train, dtype=torch.float32)
        targets = torch.tensor(y_train, dtype=torch.float32)
        outputs = model(inputs).squeeze()

        if use_dordbeta:
            phi = torch.exp(torch.tensor(1.0, requires_grad=True))  # Learnable precision
            loss = -torch.sum(dordbeta(targets, outputs, phi))
        else:
            loss = criterion(outputs, targets)

        loss.backward()
        optimizer.step()

        if epoch % 10 == 0:
            print(f"Epoch {epoch}, Loss: {loss.item():.4f}")
            
class TransformerRegressor(nn.Module):
    def __init__(self, input_dim, hidden_dim=64, nhead=4, num_layers=3):
        super(TransformerRegressor, self).__init__()
        self.input_proj = nn.Linear(input_dim, hidden_dim)
        self.transformer = nn.Transformer(
            hidden_dim, nhead, num_layers, batch_first=True  # Enable batch_first
        )
        self.fc_out = nn.Linear(hidden_dim, 1)

    def forward(self, x):
        # Input shape: (batch_size, input_dim)
        x = self.input_proj(x)  # Project input to hidden dimension
        x = x.unsqueeze(1)  # Add sequence dimension: (batch_size, seq_len=1, hidden_dim)
        x = self.transformer(x, x)  # Pass through transformer
        return torch.sigmoid(self.fc_out(x[:, -1, :]))  # Ensure outputs are in [0, 1]


# Evaluate model
def evaluate_model(model, X_test, y_test):
    with torch.no_grad():
        inputs = torch.tensor(X_test, dtype=torch.float32)
        targets = torch.tensor(y_test, dtype=torch.float32)
        outputs = model(inputs).squeeze()
        mse = mean_squared_error(targets.numpy(), outputs.numpy())
        print(f"Test MSE: {mse:.4f}")
    return mse

# Main script
if __name__ == "__main__":
    # Load data
    X_train, X_test, y_train, y_test = load_data()

    # Define models
    input_dim = X_train.shape[1]
    model_baseline = TransformerRegressor(input_dim)
    model_custom = TransformerRegressor(input_dim)

    # Define optimizers
    optimizer_baseline = optim.Adam(model_baseline.parameters(), lr=0.001)
    optimizer_custom = optim.Adam(model_custom.parameters(), lr=0.001)

    # Train baseline model (MSE Loss)
    print("Training baseline model...")
    train_model(model_baseline, nn.MSELoss(), optimizer_baseline, X_train, y_train, epochs=50)

    # Train custom model (dordbeta Loss)
    print("Training custom model...")
    train_model(model_custom, None, optimizer_custom, X_train, y_train, epochs=50, use_dordbeta=True)

    # Evaluate models
    print("Evaluating baseline model...")
    mse_baseline = evaluate_model(model_baseline, X_test, y_test)

    print("Evaluating custom model...")
    mse_custom = evaluate_model(model_custom, X_test, y_test)

    print(f"Baseline MSE: {mse_baseline:.4f}, Custom MSE: {mse_custom:.4f}")
