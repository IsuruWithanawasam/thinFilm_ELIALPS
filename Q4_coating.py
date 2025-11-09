# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 11:56:34 2025

@author: Asus
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


data_50 = pd.read_csv('50deg.csv', sep=';')
data_55 = pd.read_csv('55deg.csv', sep=';')
data_60 = pd.read_csv('60deg.csv', sep=';')

# --- Prepare Arrays ---
wl_ev = np.array(data_50['wavelength (eV)'])
wl_nm = 1240 / wl_ev

p_50 = np.array(data_50['psi'])
tan_p_50= np.array(data_50['tan(Psi)'])
dt_50 = np.array(data_50['delta'])
cos_dt_50= np.array(data_50['cos(Delta)'])

p_55 = np.array(data_55['psi'])
tan_p_55= np.array(data_55['tan(psi)'])
dt_55 = np.array(data_55['delta'])
cos_dt_55= np.array(data_55['cos(delta)'])

p_60 = np.array(data_60['psi'])
tan_p_60= np.array(data_60['tan(psi)'])
dt_60 = np.array(data_60['delta'])
cos_dt_60= np.array(data_60['cos(delta)'])

###########################################func


def calculate_rho(psi, delta):
    
  psi_rad=np.radians(psi)
  delta_rad=np.radians(delta)
  rho = np.tan(psi_rad) * np.exp(1j * delta_rad)
  return rho

def calculate_n_k_from_rho(rho, phi_deg):
  # Convert angle of incidence to radians
  phi_rad = np.radians(phi_deg)
  sin_phi = np.sin(phi_rad)
  tan_phi = np.tan(phi_rad)
  
  # Calculate complex dielectric constant <epsilon> 
  term_in_brackets = (1 - rho) / (1 + rho)
  epsilon_tilde = (sin_phi**2) * (1 + (tan_phi**2) * (term_in_brackets**2))
  
  # Extract epsilon_1 and epsilon_2
  # Syllabus uses <epsilon> = epsilon_1 - i*epsilon_2 [cite: 204]
  epsilon_1 = np.real(epsilon_tilde)
  epsilon_2 = -np.imag(epsilon_tilde)
  
  # Invert to find n and k
  # Based on epsilon_1 = n^2 - k^2 [cite: 207] and epsilon_2 = 2*n*k [cite: 211]
  mag_epsilon = np.sqrt(epsilon_1**2 + epsilon_2**2)
  
  n = np.sqrt((mag_epsilon + epsilon_1) / 2)
  k = np.sqrt((mag_epsilon - epsilon_1) / 2)
  
  return n, k


plt.figure()
plt.plot(wl_nm, p_50, label='50 Degree')
plt.plot(wl_nm, p_55, label='55 Degree')
plt.plot(wl_nm, p_60, label='60 Degree')
plt.grid(True)
plt.xlabel('Wavelength (nm)')
plt.ylabel('$\Psi$')
plt.legend()
plt.show()

plt.figure()
plt.plot(wl_nm, dt_50, label='50 Degree')
plt.plot(wl_nm, dt_55, label='55 Degree')
plt.plot(wl_nm, dt_60, label='60 Degree')
plt.grid(True)
plt.xlabel('Wavelength (nm)')
plt.ylabel('relative phase difference ($\Delta$)')
plt.legend()
plt.show()


#########50
rho_50 = calculate_rho(p_50, dt_50)
n_50, k_50 = calculate_n_k_from_rho(rho_50, 50)

#########50
rho_55 = calculate_rho(p_55, dt_55)
n_55, k_55 = calculate_n_k_from_rho(rho_55, 55)

#########50
rho_60 = calculate_rho(p_60, dt_60)
n_60, k_60 = calculate_n_k_from_rho(rho_60, 60)

plt.figure()
plt.plot(wl_nm, n_50, label='50 Degree')
plt.plot(wl_nm, n_55, label='55 Degree')
plt.plot(wl_nm, n_60, label='60 Degree')
plt.grid(True)
plt.xlabel('Wavelength (nm)')
plt.ylabel('n')
plt.legend()
plt.show()

plt.figure()
plt.plot(wl_nm, k_50, label='50 Degree')
plt.plot(wl_nm, k_55, label='55 Degree')
plt.plot(wl_nm, k_60, label='60 Degree')
plt.grid(True)
plt.xlabel('Wavelength (nm)')
plt.ylabel('k')
plt.legend()
plt.show()

