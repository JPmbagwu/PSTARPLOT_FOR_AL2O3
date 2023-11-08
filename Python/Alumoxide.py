#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 16:54:39 2023

@author: johnpaulmbagwu
"""

import pandas as pd
import matplotlib.pyplot as plt
import math

# Read the Excel data file
df = pd.read_csv('alumoxide.csv', delimiter=' ')

# Get the stopping power and energy columns
stopping_power_pstar = df['TotalStp.Pow']
energy = df['KineticEnergy']

def bethe_bloch(E, Z, A, I, Z_proj, m_e):
    # Constants
    K = 0.307075  # MeV*cm^2/g
    me_c_squared = 0.511  # MeV

    # Calculate the beta value (velocity of the particle)
    beta = math.sqrt(E * (E + 2 * m_e)) / (E + m_e)

    # Initialize density correction (delta), gamma, and T_max
    delta = 0.0
    gamma = 1.0  # Initialize gamma with a default value
    T_max = 0.0

    if beta >= 0.1:
        # Calculate the gamma value (Lorentz factor)
        gamma = 1 / math.sqrt(1 - beta ** 2)

        # Calculate the maximum energy transfer (T_max)
        T_max = (2 * m_e * beta ** 2 * gamma ** 2) / (1 + 2 * gamma * (m_e / (2 * me_c_squared)) + (m_e / me_c_squared) ** 2)

        # Calculate the density correction (delta)
        X = math.log10(beta * gamma)
        C = -3.80
        m = 2
        X0 = 0.2

        if X >= X0:
            delta = 2 * math.log10(10) * X - C + m * math.log10(X / X0)

    # Handle the case where beta is less than 0.1
    if beta < 0.1:
        return 0.0  # Return 0 for stopping power

    # Calculate the stopping power (dE/dx)
    dEdx = K * Z / A * (1 / beta ** 2) * (0.5 * math.log(2 * m_e * beta ** 2 * gamma ** 2 * T_max / (I ** 2)) - beta ** 2 - delta)

    return dEdx

# Input parameters for aluminum oxide (Al2O3)
E = 1.0  # Energy of the incident particle in MeV
Z_al2o3 = 13 + 8  # Atomic number of aluminum (Z=13) + oxygen (Z=8)
A_al2o3 = 26.98 + 16  # Atomic mass of aluminum (A=26.98) + oxygen (A=16) in g/mol
I_al2o3 = 166  # Mean excitation energy of aluminum oxide in eV
Z_proj = 1  # Atomic number of the incident particle (e.g., 1 for electrons)
m_e = 0.511  # Rest mass of the electron in MeV/c^2

# Calculate the stopping power for aluminum oxide
stopping_power_al2o3 = [bethe_bloch(e, Z_al2o3, A_al2o3, I_al2o3, Z_proj, m_e) for e in energy]

# Create a figure and axis
fig, ax1 = plt.subplots()

# Plot PSTAR data
ax1.plot(energy, stopping_power_pstar, label='PSTAR', color='blue')

# Plot your calculations
ax1.plot(energy, stopping_power_al2o3, label='My Calculations', color='red')

# Set the x and y-axis scales to logarithmic
ax1.set_xscale('log')
ax1.set_yscale('log')

# Set the x-axis and y-axis limits
ax1.set_xlim(1e-2, 1e4)

# Set the y-axis label
ax1.set_ylabel('Stopping Power (MeV cm^2/g)')

# Set the x-axis label
ax1.set_xlabel('Energy (MeV)')

# Add a legend with labels for both graphs
ax1.legend()

# Show the plot
plt.show()
