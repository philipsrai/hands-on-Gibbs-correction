###################################################################              
#     CODE : Thermal Correction to Gibbs Free Energy Calculation  #
#   Method : Using molecular vibrational frequencies              #
#        By: Philips Kumar Rai                                    #
#    Malaviya National Technology Jaipur                          #
#    Email : philipsrai786@gmail.com                              #
###################################################################               

import numpy as np

#######################################
h = 6.62607015e-34         # J.s      #
kB = 1.380649e-23          # J/K      #
c = 2.99792458e10          # cm/s     #
R = 8.314462618            # J/mol·K  #
Na = 6.02214076e23         # mol-1    #
T = 298.15                 # K        #
P = 101325                 # Pa       #
#######################################
mass_amu = 129.018780
symmetry_number = 1
mass_kg = mass_amu * 1.66053906660e-27
# --- Translational ---
q_trans = ((2 * np.pi * mass_kg * kB * T) / h**2)**(3/2) * (kB * T / P)
U_trans = 1.5 * R * T
S_trans = R * (1 + 1.5 + np.log(q_trans))
print(f"S_trans : {S_trans:.6f} J/mol·K")
# --- Rotational ---
moments_amu_bohr2 = [641.28469, 1712.41424, 2129.86756]
bohr2_to_m2 = (5.29177210903e-11)**2
moments_kg_m2 = [I * 1.66053906660e-27 * bohr2_to_m2 for I in moments_amu_bohr2]
I_prod = np.prod(moments_kg_m2)
q_rot = np.sqrt(np.pi) * ((8 * np.pi**2 * kB * T) / h**2)**(1.5) * np.sqrt(I_prod) / symmetry_number
U_rot = 1.5 * R * T
S_rot = R * (1.5 + np.log(q_rot))
print(f"S_rot : {S_rot:.6f} J/mol·K")
# --- Vibrational ---
vib_freqs_cm = [11.8, 54.3, 66.9, 111.8, 184.8, 193.6, 253.2, 262.5, 509.2, 602.1, 612.0, 654.6, 771.6, 804.1, 890.9, 906.5, 945.2, 982.9, 1042.3, 
1050.6, 1141.0, 1204.6, 1269.3, 1281.0, 1293.3, 1423.0, 1472.4, 1517.8, 1545.8, 1654.2, 1791.3, 3020.3, 3285.1, 3295.2, 3310.3, 3361.2]  
U_vib = 0.0
S_vib = 0.0
for freq in vib_freqs_cm:
    theta = h * c * freq / kB  
    x = theta / T
    U_i = R * theta / (np.exp(x) - 1) + 0.5 * R * theta
    S_i = R * (x / (np.exp(x) - 1) - np.log(1 - np.exp(-x)))
    U_vib += U_i
    S_vib += S_i
print(f"S_vib : {S_vib:.6f} J/mol·K")
# --- Electronic ---
g_elec = 2  
U_elec = 0.0
S_elec = R * np.log(g_elec)
# --- Total ---
U_total = U_trans + U_rot + U_vib + U_elec
S_total = S_trans + S_rot + S_vib + S_elec
H_total = U_total + R * T
G_total = H_total - T * S_total
G_hartree = G_total / 1000 / 2625.5
print(f"Thermal correction to Gibbs Free Energy: {G_hartree:.6f} Hartree")


