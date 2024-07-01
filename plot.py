import matplotlib.pyplot as plt

# Another set of nonuniform fractions and corresponding dt values
nonuniform_fractions_log = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
dt_values_log = [1.32848e-12, 1.30447e-12, 1.30904e-12, 1.41098e-12, 1.268e-12, 1.25313e-12,
                 1.22048e-12, 1.20112e-12, 1.19602e-12, 1.21243e-12, 3.70968e-16]

# Plotting the new data on a logarithmic scale
plt.figure(figsize=(10, 6))
plt.plot(nonuniform_fractions_log, dt_values_log, marker='o', linestyle='-', color='g')
plt.title('Nonuniform Fraction vs. dt (Another Data Set) - Log Scale')
plt.xlabel('Nonuniform Fraction')
plt.ylabel('dt (s)')
# plt.yscale('log')
plt.grid(True, which="both", ls="--")
plt.show()
