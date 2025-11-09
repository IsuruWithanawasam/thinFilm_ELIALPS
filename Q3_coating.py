# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 08:39:50 2025

@author: Asus
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import find_peaks
# import math # Not needed when using numpy

# --- Load Data ---
T_1mm_data = pd.read_csv('1mm_T.Sample.Raw.csv')
R_1mm_data = pd.read_csv('1mm_thick_R.Sample.Raw.csv')

T_6mm_data = pd.read_csv('6mm_T.Sample.Raw.csv')
R_6mm_data = pd.read_csv('6mm_thick_R.Sample.Raw.csv')

data_T_TL = pd.read_csv('Nb2O4_T_TL.csv', sep=';')
std_n_tl = pd.read_csv('Nb2O5_std.csv')


# --- Prepare Arrays ---
wl = np.array(T_1mm_data['nm'])
l_1mm = 1.0 # thickness in mm
l_6mm = 6.0 # thickness in mm

# =============================================================================
# FIX 1: Convert T and R from percentage (0-100) to fraction (0-1.0)
# =============================================================================
T_1mm = np.array(T_1mm_data[' %T']) / 100.0
R_1mm = np.array(R_1mm_data[' %R']) / 100.0
T_6mm = np.array(T_6mm_data[' %T']) / 100.0
R_6mm = np.array(R_6mm_data[' %R']) / 100.0

########################################## Functions ####################################

def calculate_k(lambda_vals, l_mm, T_vals, R_vals):
  l_nm = l_mm * 1e6 # convert mm to nm
  log_argument = T_vals / ((1 - R_vals)**2)
  k = -(lambda_vals / (4 * np.pi * l_nm)) * np.log(log_argument)
  
  return k

def calculate_n(R_vals, k_vals):
  # Define A, B, and C for the quadratic equation
  A = (1 - R_vals)
  B = -2 * (1 + R_vals)
  C = (1 - R_vals) * (1 + k_vals**2)
  
  # Use numpy's sqrt function
  discriminant = B**2 - 4 * A * C
  
  # Handle any remaining numerical instability
  # Set any small negative values in discriminant to 0 before sqrt
  discriminant[discriminant < 0] = 0
  
  # Solve for the two roots
  n1 = (-B + np.sqrt(discriminant)) / (2 * A)
  n2 = (-B - np.sqrt(discriminant)) / (2 * A)
  
  # Choose the physically valid (positive) root
  n_values = np.where(n1 > 0, n1, n2)
  
  return n_values

################################################################################

# --- Plot Input Data (using 100* for display) ---

plt.figure()
plt.plot(wl, T_1mm * 100, label='T_1mm')
plt.plot(wl, T_6mm * 100, label='T_6mm')
plt.grid(True)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Transmission (%)')
plt.legend()
plt.show()

plt.figure()
plt.plot(wl, R_1mm * 100, label='One side reflectance')
plt.plot(wl, T_1mm * 100, label='Transmission')
plt.grid(True)
plt.xlabel('Wavelength (nm)')
plt.ylabel('T / R (%)')
plt.title('Transmission & Reflectance of Fused Silica')
plt.legend()
plt.show()

plt.figure()
plt.plot(wl, R_6mm * 100, label='One side reflectance')
plt.plot(wl, T_6mm * 100, label='Transmission')
plt.grid(True)
plt.xlabel('Wavelength (nm)')
plt.ylabel('T / R (%)')
plt.title('Transmission & Reflectance of BK7 Glass')
plt.legend()
plt.show()

# --- Calculate n and k ---
# (Calculating for both 1mm and 6mm)

k_1mm = calculate_k(wl, l_1mm, T_1mm, R_1mm)
n_1mm = calculate_n(R_1mm, k_1mm)

k_6mm = calculate_k(wl, l_6mm, T_6mm, R_6mm) 
n_6mm = calculate_n(R_6mm, k_6mm) 

# =============================================================================
# for 6mm sample
# =============================================================================

fig, ax1 = plt.subplots()

# Plot Refractive Index (n) on the first axis (ax1)
color = 'tab:blue'
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('Refractive Index (n)', color=color)
ax1.plot(wl, n_6mm, color=color, label='n')
ax1.tick_params(axis='y', labelcolor=color)
ax1.grid(True, linestyle='--')

# Create a second y-axis that shares the same x-axis
ax2 = ax1.twinx()  

# Plot Extinction Coefficient (k) on the second axis (ax2)
color = 'tab:red'
ax2.set_ylabel('Extinction Coeff (k)', color=color)
ax2.plot(wl, k_6mm, color=color, label='k')
ax2.tick_params(axis='y', labelcolor=color)
# Optional: Use scientific notation for the k-axis
# ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# Add a combined legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc='upper right')
plt.title('Refractive Index (n) and Extinction (k) of BK7 Glass')
plt.xlim(250, 700)
plt.show()


# =============================================================================
# for 1mm sample
# =============================================================================


fig, ax1 = plt.subplots()

# Plot Refractive Index (n) on the first axis (ax1)
color = 'tab:blue'
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('Refractive Index (n)', color=color)
ax1.plot(wl, n_1mm, color=color, label='n (1mm sample)')
ax1.tick_params(axis='y', labelcolor=color)
ax1.grid(True, linestyle='--')

# Create a second y-axis that shares the same x-axis
ax2 = ax1.twinx()  

# Plot Extinction Coefficient (k) on the second axis (ax2)
color = 'tab:red'
ax2.set_ylabel('Extinction Coeff (k)', color=color)
ax2.plot(wl, k_1mm, color=color, label='k (1mm sample)')
ax2.tick_params(axis='y', labelcolor=color)
# Optional: Use scientific notation for the k-axis
# ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# Add a combined legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc='upper right')
plt.title('Refractive Index (n) and Extinction (k) of Fused Silica')
plt.xlim(250, 700)
plt.show()

############################### Swanepole ########################################


T_swan = np.array(data_T_TL['T'])/100.0
TL_swan = np.array(data_T_TL['Lower Envelope'])/100.0

n_sub= n_6mm


def calculate_N(n_sub, TL_vals):
    N = (2 * n_sub / TL_vals) - ((n_sub**2 + 1) / 2)
    return N

def calculate_n_layer(N_vals, n_sub):
    # Calculate the discriminant
    discriminant = N_vals**2 - n_sub**2
    
    # FIX: Set any small negative values (from noise) to 0 to prevent sqrt error
    discriminant[discriminant < 0] = 0
    
    # Calculate n_l
    n_l = np.sqrt(N_vals + np.sqrt(discriminant))
    
    # Clean up any remaining NaNs that might come from n_sub
    n_l = np.nan_to_num(n_l)
    
    return n_l


peaks, properties = find_peaks(T_swan, height=0.8, distance=5)

# 5. Plot the input layer data
plt.figure()
plt.plot(wl, T_swan, label='T') # Plot as %
plt.plot(wl, TL_swan, 'r--', label='TL') # Plot as %
#plt.plot(wl[peaks], T_swan[peaks], "rx", label="Peaks")
plt.grid(True)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Transmission (%)')
plt.title('Layer Transmission and Envelope')
plt.legend(loc='lower right')
plt.xlim(300, max(wl))
plt.show()


indices = np.where(wl > 350) # Using 350nm is a safer choice

wl_calc = wl[indices]
n_sub_calc = n_sub[indices]
TL_calc = TL_swan[indices]


N_values = calculate_N(n_sub_calc, TL_calc)
n_layer = calculate_n_layer(N_values, n_sub_calc)

wl_std = (np.array(std_n_tl['wl']))*1000
n_std = np.array(std_n_tl['n'])

plt.figure()
plt.plot(wl_calc, n_layer, label='$n_l$ (Layer)')
plt.plot(wl_std, n_std, label='from lit')
plt.grid(True)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Refractive Index ($n$)')
plt.title('Layer Refractive Index (Swanepoel Method)')
plt.xlim(350, 650) # Start plot from the transparent region
plt.legend()
plt.show()


wl_max_ar=wl[peaks]
n_ar=n_layer[peaks]

print(wl_max_ar)
print(n_ar)

wl1=wl_max_ar[2]
wl2=wl_max_ar[1]
n1=n_ar[2]
n2=n_ar[1]

d=(wl1*wl2)/(2*(n1*wl2-n2*wl1))
print('d= ',d,' nm')





