import numpy as np

# Define constants
s = 1.0  # Example value for s
t = 0.5  # Example value for t

# Initial values for the currents
i_an_init = 0.01  # Initial anodic current
i_ca_init = 0.02  # Initial cathodic current

# Increment in currents
delta_i_an = 0.001  # Small increment for anodic current
delta_i_ca = 0.001  # Small increment for cathodic current

# Function to calculate incremental resistances
def incremental_resistance(i_an, i_ca, s, t):
    # Anodic incremental resistance
    R_an = s / (t * i_an + t)
    # Cathodic incremental resistance
    R_ca = s / (t * i_ca + t)
    return R_an, R_ca

# Calculate initial incremental resistances
R_an_init, R_ca_init = incremental_resistance(i_an_init, i_ca_init, s, t)

print(f"Initial R_an: {R_an_init:.4f}, Initial R_ca: {R_ca_init:.4f}")

# Increment the currents
i_an_new = i_an_init + delta_i_an
i_ca_new = i_ca_init + delta_i_ca

# Calculate new incremental resistances after the increment
R_an_new, R_ca_new = incremental_resistance(i_an_new, i_ca_new, s, t)

print(f"New R_an after increment: {R_an_new:.4f}, New R_ca after increment: {R_ca_new:.4f}")

# Calculate change in potentials (u_an, u_ca)
def potential_change(i_an, i_ca, s, t):
    # Calculate the potential based on the equation
    u_an = s * np.log((1 / t) * i_an + 1)
    u_ca = s * np.log((1 / t) * i_ca + 1)
    return u_an, u_ca

# Calculate initial potentials
u_an_init, u_ca_init = potential_change(i_an_init, i_ca_init, s, t)
print(f"Initial u_an: {u_an_init:.4f}, Initial u_ca: {u_ca_init:.4f}")

# Calculate new potentials after increment
u_an_new, u_ca_new = potential_change(i_an_new, i_ca_new, s, t)
print(f"New u_an after increment: {u_an_new:.4f}, New u_ca after increment: {u_ca_new:.4f}")

# Calculate changes in potentials
delta_u_an = u_an_new - u_an_init
delta_u_ca = u_ca_new - u_ca_init

print(f"Change in u_an: {delta_u_an:.4f}, Change in u_ca: {delta_u_ca:.4f}")
