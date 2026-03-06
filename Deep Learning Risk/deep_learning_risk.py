# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 05:17:47 2026

@author: nburgessx
"""

import torch
import torch.nn as nn
import torch.optim as optim

# -----------------------------
# 1. Create training data
# -----------------------------
# X_train: each row is a market scenario, columns are input features
# For demo: 10000 scenarios, 12 features (5 yield curve points, 5 credit spreads, 2 correlation params)
X_train = torch.rand(10000, 12)  # random numbers for illustration

# y_train: corresponding portfolio prices (Nth-to-Default PVs)
# In reality, these would come from Monte Carlo simulations
y_train = torch.rand(10000, 1)

# -----------------------------
# 2. Define the neural network
# -----------------------------
# The network learns the pricing kernel: mapping market inputs -> price
class PricingKernelNN(nn.Module):
    def __init__(self, input_size):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_size, 32),  # input layer -> 32 neurons
            nn.ReLU(),                  # activation introduces nonlinearity
            nn.Linear(32, 16),          # hidden layer -> 16 neurons
            nn.ReLU(),
            nn.Linear(16, 1)            # output layer -> predicted price
        )
    def forward(self, x):
        return self.net(x)

# Instantiate the model with 12 input features
model = PricingKernelNN(X_train.shape[1])

# -----------------------------
# 3. Training setup
# -----------------------------
criterion = nn.MSELoss()          # mean squared error loss (price differences)
optimizer = optim.Adam(model.parameters(), lr=0.01)  # Adam optimizer

# -----------------------------
# 4. Training loop
# -----------------------------
for epoch in range(20):           # loop over epochs
    optimizer.zero_grad()          # reset gradients
    y_pred = model(X_train)        # forward pass: predict prices
    loss = criterion(y_pred, y_train)  # compute loss
    loss.backward()                # backward pass: compute gradients
    optimizer.step()               # update weights
    if epoch % 5 == 0:             # print loss every 5 epochs
        print(f"Epoch {epoch}, Loss: {loss.item():.6f}")

# -----------------------------
# 5. Real-time pricing & risk
# -----------------------------
# X_new: new market scenarios for pricing
# requires_grad=True allows us to compute sensitivities (gradients) later
X_new = torch.rand(5, 12, requires_grad=True)

# Predict prices for new scenarios
predicted_prices = model(X_new)

# Compute gradients to get sensitivities
# Automatic differentiation computes how price changes with each input
predicted_prices.sum().backward()  # sum ensures scalar output for grad

# CS01: sensitivity of price to credit spreads (columns 5-9)
CS01 = X_new.grad[:, 5:10]

# Correlation Delta: sensitivity to correlation parameters (columns 10-11)
Correlation_Delta = X_new.grad[:, 10:]

# Output results
print("Predicted Prices:", predicted_prices.detach().numpy())
print("CS01 (credit spread sensitivity):", CS01.detach().numpy())
print("Correlation Delta:", Correlation_Delta.detach().numpy())