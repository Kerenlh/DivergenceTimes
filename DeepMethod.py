import torch
import torch.nn as nn
import torch.optim as optim
import torch.distributions as dist
import random
import numpy as np
from torch.nn import init
import math
import matplotlib.pyplot as plt

# Define the DeepSet model
class DeepSetModel(nn.Module):
    def __init__(self, input_dim, hidden_dim):
        super(DeepSetModel, self).__init__()
        self.encoder = nn.Sequential(
            EquivLayer(input_dim, hidden_dim),  # Updated input dimension
            nn.ReLU(),
            EquivLayer(hidden_dim, hidden_dim),
            nn.ReLU(),
            EquivLayer(hidden_dim, hidden_dim),
            nn.ReLU(),
            EquivLayer(hidden_dim, 1),
        )

    def forward(self, X, Z):
        # Aggregate information from sets of (X_i, Z_i) pairs
        combined = torch.stack((X, Z), dim=2)  # Concatenate along the feature dimension
        #print(combined.shape)
        # Feed the aggregated information through the encoder
        encoded = self.encoder(combined)
        return torch.mean(encoded.squeeze(), dim=1)

# Function to generate the dataset
def generate_dataset(n_samples, n_features):
    X_list = []
    Z_list = []
    p_list = []

    for _ in range(n_samples):
        shape = 0.2305047*torch.ones(n_features) # alpha estimated from XmtDNA
        rate =  0.1643353*torch.ones(n_features) # beta estimated from XmtDNA

        gamma_dist = dist.Gamma(shape, rate)
        lambdas = gamma_dist.sample()  # Generate lambdas from gamma distribution

        #lambdas = torch.rand(n_features) * 10
        # Simulate X ~ Poisson(lambdas)
        X = torch.poisson(lambdas)

        # Sample p
        p = torch.rand(1) * 2  # Random p between 0 and 2

        # Simulate Y ~ Poisson(lambdas * p)
        Y = torch.poisson(lambdas * p)

        # Determine Z as the parity of Y
        Z = torch.remainder(Y, 2)

        X_list.append(X)
        Z_list.append(Z)
        p_list.append(p)

    return torch.stack(X_list), torch.stack(Z_list), torch.stack(p_list)



# Function to generate the dataset
def generate_dataset2(n_samples, n_features):
    X_list = []
    Z_list = []
    p_list = []

    for _ in range(n_samples):
        # Sample p
        p = torch.rand(1) * 2  # Random p between 0 and 2
        X = torch.normal(p*torch.ones(n_features))
        Z = X*0

        X_list.append(X)
        Z_list.append(Z)
        p_list.append(p)

    return torch.stack(X_list), torch.stack(Z_list), torch.stack(p_list)

# Generate training and testing datasets
n_train_samples = 10000
n_test_samples = 100
n_features = 16569  # Length of X and Lambdas vectors

X_train, Z_train, p_train = generate_dataset(n_train_samples, n_features)
print(X_train.shape, Z_train.shape, p_train.shape)
X_test, Z_test, p_test = generate_dataset(n_test_samples, n_features)


# Define DeepSet model
input_dim = n_features  # Dimensionality of X and Z vectors
hidden_dim = 4
model = DeepSetModel(input_dim=2, hidden_dim=hidden_dim)


# Define loss function and optimizer
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.0001)
#optimizer = optim.SGD(model.parameters(), lr=0.0001, momentum=0.9)


# Define batch size
batch_size = 32

res = []
# Training loop
num_epochs = 200000
for epoch in range(num_epochs):
    epoch_loss = 0.0
    # for i in range(0, 10):
    indices = np.array(random.sample(range(n_train_samples), batch_size))
    # indices = torch.randperm(n_train_samples)[:batch_size]

    optimizer.zero_grad()

    batch_X = X_train[indices].float()
    batch_Z = Z_train[indices].float()
    batch_p = p_train[indices].float()

    output = model(batch_X, batch_Z).squeeze()
    #loss = criterion(output, batch_p)
    loss = torch.mean((output - batch_p.squeeze())**2)
    loss.backward()
    optimizer.step()

    epoch_loss += loss.item() * batch_X.size(0)

    epoch_loss /= len(X_train)
    res.append(epoch_loss)

    if epoch % 1000 == 0:
#         plt.plot(res[-1000:])
#         plt.show()
        print(f"Epoch [{epoch}/{num_epochs}], Loss: {np.mean(res[-200:]):.8f}")

    if epoch % 10000 == 9999:
        # Evaluate on testing data and print real and estimated p values
        with torch.no_grad():
            estimated_p = model(X_test.float(), Z_test.float()).squeeze()
            print("Real p vs Estimated p:")
            for real_p, est_p in zip(p_test.squeeze(), estimated_p):
                print(f"Real p: {real_p.item():.3f}, Estimated p: {est_p.item():.3f}")

        # Calculate and print total MSE
        total_mse = torch.mean((estimated_p.squeeze() - p_test.float().squeeze()) ** 2).item()
        print("Total MSE:", total_mse)
