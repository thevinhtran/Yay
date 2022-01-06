""" Sound Speeds Helper Files """ 

import numpy as np 
import pandas as pd
import sympy as sym 



""" 
    Establishing the relevant classes to our system
"""

class eos:
    # equation of state coupling constants 
    def __init__(self,\
                g_sigma_N = 0.0, g_omega_N = 0.0, g_rho_N = 0.0, g_phi_N = 0.0, b = 0.0, c = 0.0,\
                g_sigma_H = 0.0, g_omega_H = 0.0, g_rho_H = 0.0, g_phi_H = 0.0):
        
        self.g_sigma_N = g_sigma_N
        self.g_omega_N = g_omega_N
        self.g_rho_N = g_rho_N
        self.g_phi_N = g_phi_N

        self.b = b
        self.c = c

        self.g_sigma_H = g_sigma_H
        self.g_omega_H = g_omega_H
        self.g_rho_H = g_rho_H
        self.g_phi_H = g_phi_H


class baryon:
    # baryon class, stores both symbolic and numeric values in a single class 
    # baryon particle 
    def __init__(self, name, spin, isospin, charge, kind, var_type,\
                sym_mass, sym_mass_eff, sym_density, sym_frac, sym_kf, sym_ef, sym_chem_pot,\
                num_mass = 0.0, num_mass_eff = 0.0, num_density = 0.0, num_frac = 0.0, num_kf = 0.0, num_ef = 0.0, num_chem_pot = 0.0):

        # variables common to both classes
        self.name = name
        self.kind = kind
        self.var_type = var_type 
        self.charge = charge 
        self.spin = spin 
        self.isospin = isospin 

        # variables to be established at baryon declaration
        self.sym_mass = sym_mass
        self.num_mass = num_mass

        # variables to be stored later
        self.sym_mass_eff = sym_mass_eff
        self.sym_num_density = sym_density
        self.sym_frac = sym_frac
        self.sym_kf = sym_kf
        self.sym_ef = sym_ef
        self.sym_chem_pot = sym_chem_pot

        self.num_mass_eff = num_mass_eff
        self.num_num_density = num_density
        self.num_frac = num_frac
        self.num_kf = num_kf
        self.num_ef = num_ef
        self.num_chem_pot = num_chem_pot
        

        # coupling constants 
        self.sym_g_sigma = 0.0
        self.sym_g_omega = 0.0 
        self.sym_g_rho = 0.0 
        self.sym_g_phi = 0.0 

        self.num_g_sigma = 0.0
        self.num_g_omega = 0.0 
        self.num_g_rho = 0.0 
        self.num_g_phi = 0.0



class lepton: 



class meson: 