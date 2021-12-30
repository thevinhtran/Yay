""" Module that defines particle classes and provides initial list of baryons, leptons, mesons 
    Hopefully would allow for ease of use of declaring particles
"""

import numpy as np
import pandas as pd
import sympy as sym
import matplotlib.pyplot as plt


pi = np.pi
Pi = sym.symbols('pi')
hc = 197.33 # MeV fm
n_sat = 0.16 # fm-3

U = sym.Function('U')

class eos:
    # equation of state class
    # can use to declare an EOS object that contains the coupling constants for the model 

    def __init__(self, g_sigma_N = 0.0, g_omega_N = 0.0, g_rho_N = 0.0, g_phi_N = 0.0, b = 0.0, c = 0.0,\
                    g_sigma_H = 0.0, g_omega_H = 0.0, g_rho_H = 0.0, g_phi_H = 0.0):
        
        self.g_sigma_N = g_sigma_N
        self.g_omega_N = g_omega_N
        self.g_rho_N = g_rho_N
        self.g_phi_N = g_phi_N
        
        self.g_sigma_H = g_sigma_H
        self.g_omega_H = g_omega_H
        self.g_rho_H = g_rho_H
        self.g_phi_H = g_phi_H
        
        self.b = b
        self.c = c

        # mixed coupling
        # self.g_sigma_omega... 
        # self.g_sigma_rho... 
        # self.g_omega_rho... etc, etc 


class baryon:
    # baryon class
    # use to declare a baryon object that holds attributes of that particle such as mass, charge, etc 
    def __init__(self, mass, isospin, charge, kind, var_type, g_sigma, mass_eff = 0.0, num_density = 0.0,\
                 frac = 0.0, kf = 0.0, ef = 0.0, chem_pot = 0.0):
    
        # variables to be established at baryon declaration
        self.mass = mass
        self.isospin = isospin
        self.charge = charge
        self.kind = kind
        self.var_type = var_type
        self.g_sigma = g_sigma
    
        # variables to be stored later
        self.mass_eff = mass_eff
        self.num_density = num_density
        self.frac = frac
        self.kf = kf
        self.ef = ef
        self.chem_pot = chem_pot


class lepton: 
    def __init__(self, mass, charge, num_density = 0.0, frac = 0.0, var_type = 0.0, kf = 0.0, chem_pot = 0.0):
        self.mass = mass
        self.charge = charge

        self.num_density = num_density
        self.frac = frac
        self.var_type = var_type
        self.kf = kf
        self.chem_pot = chem_pot 


class meson:
    # meson class: initialize with mass and coupling constants 
    def __init__(self, mass, g_nucleon = 0.0, g_hyperon = 0.0, field = 0.0):
        self.mass = mass # in MeV
        self.g_nucleon = g_nucleon
        self.g_hyperon = g_hyperon
        self.field = 0.0  


class independent_var:
    # independent variable class
    def __init__(self, var, func, total_deriv = 0.0, tilde_chem_pot = 0.0, tilde_chem_pot_val = 0.0, num_val = 0.0, total_deriv_num = 0.0):
        self.var = var # stores symbol for independent varible
        self.func = func # not sure why I defined this one

        self.total_deriv = total_deriv # total derivative wrt to nB, used in the system of equations part
        
        self.tilde_chem_pot = tilde_chem_pot
        self.tilde_chem_pot_val = tilde_chem_pot_val
        
        self.num_val = num_val
        
        # total derivative of fraction wrt to nB
        self.total_deriv_num = total_deriv_num




# would be nice if in the future we could already have these defined here so to initialize all a user would need to do would be to 
# call a function and just pass a baryon name... 
Neutron = [sym.symbols('m_n'), sym.symbols('n_n'), sym.symbols('n_B')*sym.symbols('x_n'), sym.symbols('x_n'), sym.symbols('mu_n'), 0.0]
Proton = [sym.symbols('m_p'), sym.symbols('n_p'), sym.symbols('n_B')*sym.symbols('x_p'), sym.symbols('x_p'), sym.symbols('mu_p'), 1.0]
Lambda_0 = [sym.symbols('n_Lambda_0'), sym.symbols('n_B')*sym.symbols('x_Lambda_0'), sym.symbols('x_Lambda_0'),\
             sym.symbols('mu_Lambda_0'),  0.0]
Sigma_0 = [sym.symbols('n_Sigma_0'), sym.symbols('n_B')*sym.symbols('x_Sigma_0'), sym.symbols('x_Sigma_0'),\
           sym.symbols('mu_Sigma_0'), 0.0]
Sigma_neg = [sym.symbols('n_Sigma_-'), sym.symbols('n_B')*sym.symbols('x_Sigma_-'), sym.symbols('x_Sigma_-'),\
           sym.symbols('mu_Sigma_-'), 0.0]
Cascade_neg = [sym.symbols('n_\Xi_-'), sym.symbols('n_B')*sym.symbols('x_\Xi-'), sym.symbols('x_\Xi-'),\
               sym.symbols('mu_\Xi_-'), -1.0]




# MESONS 

#Sigma = 
#Omega = 
#Rho = 
#phi = 


