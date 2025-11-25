import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler

# Load the CSV file
file_path = 'dataThyVC.csv'  # Replace with your file path
data = pd.read_csv(file_path)

# %%
# Set the font to Times and enable LaTeX rendering with pdflatex
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times']
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{newtxmath}'
plt.rcParams['text.latex.preamble'] = r'\usepackage{gensymb}'
plt.rcParams['pgf.rcfonts'] = False
plt.rcParams['pgf.texsystem'] = 'pdflatex'

plt.rcParams["axes.prop_cycle"] = cycler(color="brgcmyk")

# Extract relevant data
time = data['Time / s']
voltage = data['Vm3:Measured voltage']
current = data['Am1:Measured current']
# current = data['I_dc:Source current']

# Create a figure with two subplots
fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Subplot 1: Voltage vs Time
ax[0].plot(time, voltage, color='b', label='Voltage',  lw=1.5)
ax[0].set_ylabel('Voltage (V)')
ax[0].set_xlim([0, 0.1])
#ax[0].set_ylim([30, 35])
ax[0].grid(True)

# Subplot 2: Current vs Time
ax[1].plot(time, current, color='r', label='Current',  lw=1.5)
ax[1].set_xlabel('Time (s)')
ax[1].set_ylabel('Current (A)')
ax[1].set_xlim([0, 0.1])
# ax[1].set_ylim([58, 62])
#ax[1].set_yticks(range(58, 62, 1))
ax[1].grid(True)

# Adjust layout and show the plot
plt.tight_layout()
plt.show()
