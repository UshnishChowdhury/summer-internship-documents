import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler

# Load the CSV file
file_path = 'datathd.csv'  # Replace with your file path
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

# Data for the plots
frequency = data['Frequency / Hz']
voltage = data['Vm3:Measured voltage']
current = data['Am1:Measured current']

# Create a figure with subplots
fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Subplot 1: Voltage vs Frequency
axes[0].bar(frequency, voltage, color='b', alpha=0.7, width=20)
#axes[0].set_title('Voltage')
axes[0].set_ylabel('Voltage (V)')
axes[0].grid(True, alpha=0.5)
axes[0].set_xticks(range(0, 600, 50))
axes[0].set_xlim([0, 600])

# Subplot 2: Current vs Frequency
axes[1].bar(frequency, current, color='r', alpha=0.7, width=20)
#axes[1].set_title('Current')
axes[1].set_xlabel('Frequency (Hz)')
axes[1].set_ylabel('Current (A)')
axes[1].grid(True, alpha=0.5)
axes[1].set_xticks(range(0, 600, 50))
axes[1].set_xlim([0, 600])
# Adjust layout and show the plot
plt.tight_layout()
plt.show()
