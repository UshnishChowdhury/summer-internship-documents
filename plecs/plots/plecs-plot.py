import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler

# Load the CSV file
file_path_15 = 'data15.csv'  # Replace with the actual file path if different
data_15 = pd.read_csv(file_path_15)

file_path_35 = 'data35.csv'  # Replace with the actual file path if different
data_35 = pd.read_csv(file_path_35)

file_path_65 = 'data65.csv'  # Replace with the actual file path if different
data_65 = pd.read_csv(file_path_65)

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

# Extract Voltage and Current data
voltage_15 = data_15['Vm3:Measured voltage']
current_15 = data_15['I:Source current']

voltage_35 = data_35['Vm3:Measured voltage']
current_35 = data_35['I:Source current']

voltage_65 = data_65['Vm3:Measured voltage']
current_65 = data_65['I:Source current']

# Create the plot
#plt.figure(figsize=(10, 6))
plt.plot(current_15, voltage_15,  label=r'15$^{\circ}$ C', lw=1.5)
plt.plot(current_35, voltage_35, color='red', label=r'35$^{\circ}$ C', lw=1.5)
plt.plot(current_65, voltage_65, color='green', label=r'65$^{\circ}$ C', lw=1.5)

# Add labels and title
plt.xlabel(f'Current (A)')
plt.ylabel(f'Voltage (V)')
plt.xlim([0, 120])
plt.ylim([26, 42])
plt.title(f'Constant pressure: 25 bar')
plt.legend(loc='lower right')

# Add grid and legend
plt.grid()
# Show the plot
plt.tight_layout()
plt.show()