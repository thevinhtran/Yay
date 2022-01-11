import numpy as np
import pandas as pd
import sympy as sym
import matplotlib.pyplot as plt 


# coupling constant symbols 
Pi = sym.symbols('pi')
hc = 197.33
g_sigma_N_sym, g_omega_N_sym, g_rho_N_sym, g_phi_N_sym = sym.symbols('g_sigma_N, g_omega_N, g_rho_N, g_phi_N')
b_sym, c_sym = sym.symbols('b, c')
g_sigma_H_sym, g_omega_H_sym, g_rho_H_sym, g_phi_H_sym = sym.symbols('g_sigma_H, g_omega_H, g_rho_H, g_phi_H')




# minimal classes 

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
    # lepton particle 
    def __init__(self, name, charge, var_type,\
                sym_mass, sym_density, sym_frac, sym_kf, sym_chem_pot,\
                num_mass = 0.0, num_density = 0.0, num_frac = 0.0, num_kf = 0.0, num_chem_pot = 0.0):
        self.name = name
        self.charge = charge
        self.var_type = var_type


        self.sym_mass = sym_mass
        self.sym_density = sym_density
        self.sym_frac = sym_frac
        self.sym_kf = sym_kf
        self.sym_chem_pot = sym_chem_pot

        self.num_mass = num_mass
        self.num_density = num_density
        self.num_frac = num_frac
        self.num_kf = num_kf
        self.num_chem_pot = num_chem_pot

def partial_omega(nb, baryon_list, ind_var):
    result = 0
    for baryon in baryon_list:
        result += baryon.sym_g_omega * nb * baryon.sym_frac 
    result = result/omega.sym_mass**2
    return result.diff(ind_var) 

def partial_rho(nb, baryon_list, ind_var):
    result = 0 
    for baryon in baryon_list:
        result += baryon.sym_g_rho * baryon.sym_frac * nb * baryon.isospin 
    result = result/rho.sym_mass**2 
    return result.diff(ind_var) 

def partial_phi(nb, baryon_list, ind_var):
    result = 0 
    for baryon in baryon_list: 
        result += baryon.sym_g_phi * nb * baryon.sym_frac 
    result = result/omega.sym_mass**2
    return result.diff(ind_var)



class meson:
    def __init__(self, name, sym_mass, sym_field, num_mass, num_field = 0.0, sym_g_N = 0.0, sym_g_H = 0.0, num_g_N = 0.0, num_g_H = 0.0, partial = partial_omega):
        self.name = name 
        self.sym_mass = sym_mass # in MeV
        self.sym_field = sym_field

        self.num_mass = num_mass
        self.num_field = num_field

        # coupling constants 
        # could update this in the future to store coupling constants here
        # but for now only makes sense for the sigma meson to store self coupling 

        self.sym_b = sym.symbols('b') 
        self.sym_c = sym.symbols('c')
        self.num_b = 0.0 
        self.num_c = 0.0 

        # eos coupling
        self.sym_g_N = sym_g_N 
        self.sym_g_H = sym_g_H
        self.num_g_N = num_g_N
        self.num_g_H = num_g_H 

        #
        self.partial = partial 



# symbolic baryon objects
Proton = baryon(name = 'Proton', spin = 1/2, isospin = 1/2, charge = 1, kind = 'Nucleon', var_type = 'Dependent',\
                sym_mass = sym.symbols('m_p'), sym_mass_eff = sym.symbols('m_p^*'), sym_density = sym.symbols('n_p'),\
                sym_frac = sym.symbols('x_p'), sym_kf = sym.symbols('k_p'), sym_ef = sym.symbols('E^*_p'), sym_chem_pot = sym.symbols('mu_p'),\
                num_mass = 939.0)

Neutron = baryon(name = 'Neutron', spin = 1/2, isospin = -1/2, charge = 0.0, kind = 'Nucleon', var_type = 'Dependent',\
                sym_mass = sym.symbols('m_n'), sym_mass_eff = sym.symbols('m_n^*'), sym_density = sym.symbols('n_n'),\
                sym_frac = sym.symbols('x_n'), sym_kf = sym.symbols('k_n'), sym_ef = sym.symbols('E^*_n'), sym_chem_pot = sym.symbols('mu_n'),\
                num_mass = 939.0)

Lambda = baryon(name = 'Lambda', spin = 1/2, isospin = 0, charge = 0, kind = 'Hyperon', var_type = 'Indepdent',\
                sym_mass = sym.symbols('m_Lambda'), sym_mass_eff =  sym.symbols('m_Lambda^*'), sym_density = sym.symbols('n_Lambda'),\
                sym_frac =  sym.symbols('x_Lambda'), sym_kf = sym.symbols('k_Lambda'), sym_ef = sym.symbols('E^*_Lambda'), sym_chem_pot = sym.symbols('mu_Lambda'),\
                num_mass = 1116.0)

Sigma_neu = baryon(name = 'Sigma_0', spin = 1/2, isospin = 1.0, charge = 0.0, kind = 'Hyperon', var_type = '',\
                    sym_mass = sym.symbols('m_Sigma_0'), sym_mass_eff = sym.symbols('m_Sigma_0^*'),\
                    sym_density = sym.symbols('n_Sigma_0'), sym_frac = sym.symbols('x_Sigma_0'), sym_kf = sym.symbols('k_Sigma'),\
                    sym_ef = sym.symbols('E^*_Sigma'), sym_chem_pot = sym.symbols('mu_Sigma'),\
                    num_mass = 1192.642)

Sigma_plus = baryon(name = 'Sigma_+', spin = 1/2, isospin = 1.0, charge = 1.0, kind = 'Hyperon', var_type = '',\
                    sym_mass = sym.symbols('m_Sigma_+'), sym_mass_eff = sym.symbols('m_Sigma_+^*'),\
                    sym_density = sym.symbols('n_Sigma_+'), sym_frac = sym.symbols('x_Sigma_+'), sym_kf = sym.symbols('k_Sigma_+'),\
                    sym_ef = sym.symbols('E^*_Sigma_+'), sym_chem_pot = sym.symbols('mu_Sigma_+'),\
                    num_mass = 1189.37)

Sigma_min = baryon(name = 'Sigma_-', spin = 1/2, isospin = 1.0, charge = -1.0, kind = 'Hyperon', var_type = '',\
                    sym_mass = sym.symbols('m_Sigma_-'), sym_mass_eff = sym.symbols('m_Sigma_-^*'),\
                    sym_density = sym.symbols('n_Sigma_-'), sym_frac = sym.symbols('x_Sigma_-'), sym_kf = sym.symbols('k_Sigma_-'),\
                    sym_ef = sym.symbols('E^*_Sigma_-'), sym_chem_pot = sym.symbols('mu_Sigma_-'),\
                    num_mass = 1197.5)

xi_neu_sym = baryon(name = 'Xi_0', spin = 1/2, isospin = 1/2, charge = 0.0, kind = 'Hyperon', var_type = '',\
                    sym_mass = sym.symbols('m_Xi_0'), sym_mass_eff =  sym.symbols('m_Xi_0^*'),\
                    sym_density = sym.symbols('n_Xi_0'), sym_frac = sym.symbols('x_Xi_0'), sym_kf =  sym.symbols('k_Xi_0'),\
                    sym_ef = sym.symbols('E^*_Xi_0'), sym_chem_pot = sym.symbols('mu_Xi_0'),\
                    num_mass = 1314.86)

xi_min_sym = baryon(name = 'Xi_-', spin = 1/2, isospin = 1/2, charge = -1.0, kind = 'Hyperon', var_type = '',\
                    sym_mass = sym.symbols('m_Xi_-'), sym_mass_eff = sym.symbols('m_Xi_-^*'),\
                    sym_density = sym.symbols('n_Xi_-'), sym_frac = sym.symbols('x_Xi_-'), sym_kf = sym.symbols('k_Xi_-'),\
                    sym_ef = sym.symbols('E^*_Xi_-'), sym_chem_pot = sym.symbols('mu_Xi_-'),\
                    num_mass = 1321.72)


""" Leptons """

# symbolic lepton objects 
electron = lepton(name = 'electron', charge = -1.0, var_type = 'Independent',\
                sym_mass = sym.symbols('m_e'), sym_density = sym.symbols('n_e'), sym_frac = sym.symbols('x_e'),\
                sym_kf = sym.symbols('k_e'), sym_chem_pot = sym.symbols('\mu_e'),\
                num_mass = 0.510)

muon = lepton(name = 'muon', charge = -1.0, var_type = 'Dependent',\
            sym_mass = sym.symbols('m_mu'), sym_density = sym.symbols('n_mu'), sym_frac = sym.symbols('x_\mu'),\
            sym_kf = sym.symbols('k_mu'), sym_chem_pot =  sym.symbols('\mu_mu'),\
            num_mass = 105.65) 




""" Mesons """

# symbolic meson objects 
sigma = meson(name = 'sigma', sym_mass = sym.symbols('m_sigma'), sym_field = sym.symbols('sigma'), num_mass = 550.0,\
            sym_g_N = sym.symbols('g_sigma_N'), sym_g_H = sym.symbols('g_sigma_H'))
omega = meson(name = 'omega', sym_mass = sym.symbols('m_omega'), sym_field = sym.symbols('omega'), num_mass = 783.0,\
            sym_g_N = sym.symbols('g_omega_N'), sym_g_H = sym.symbols('g_omega_H'), partial = partial_omega)
rho = meson(name = 'rho', sym_mass = sym.symbols('m_rho'), sym_field = sym.symbols('rho'), num_mass = 770.0,\
            sym_g_N = sym.symbols('g_rho_N'), sym_g_H = sym.symbols('g_rho_H'), partial = partial_rho)
phi = meson(name = 'phi', sym_mass = sym.symbols('m_phi'), sym_field = sym.symbols('phi'), num_mass = 1020.0,\
            sym_g_N = sym.symbols('g_phi_N'), sym_g_H = sym.symbols('g_phi_H'), partial = partial_phi)





def baryon_coupling(baryon, eos ):
    if (baryon.kind == 'Nucleon'):
        baryon.sym_g_sigma = g_sigma_N_sym
        baryon.sym_g_omega = g_omega_N_sym
        baryon.sym_g_rho = g_rho_N_sym
        baryon.sym_g_phi = g_phi_N_sym

        baryon.num_g_sigma = eos.g_sigma_N
        baryon.num_g_omega = eos.g_omega_N
        baryon.num_g_rho = eos.g_rho_N
        baryon.num_g_phi = eos.g_phi_N

    elif (baryon.kind == 'Hyperon'):
        baryon.sym_g_sigma = g_sigma_H_sym
        baryon.sym_g_omega = g_omega_H_sym
        baryon.sym_g_rho = g_rho_H_sym
        baryon.sym_g_phi = g_phi_H_sym

        baryon.num_g_sigma = eos.g_sigma_H
        baryon.num_g_omega = eos.g_omega_H
        baryon.num_g_rho = eos.g_rho_H
        baryon.num_g_phi = eos.g_phi_H


# load in sigma field self coupling 
def sigma_coupling(eos):
    sigma.num_b = eos.b
    sigma.num_c = eos.c 

# meson load 
def meson_coupling(eos, meson_list):
    for meson in meson_list:
        if (meson == sigma):
            meson.num_g_N = eos.g_sigma_N
            meson.num_g_N = eos.g_sigma_H
        elif (meson == omega):
            meson.num_g_N = eos.g_omega_N
            meson.num_g_N = eos.g_omega_H
        elif (meson == rho):
            meson.num_g_N = eos.g_rho_N
            meson.num_g_N = eos.g_rho_H
        elif (meson == phi): 
            meson.num_g_N = eos.g_phi_N
            meson.num_g_N = eos.g_phi_H

def init_system(eos, meson_list, baryon_list):
    for baryon in baryon_list:
        baryon_coupling(baryon, eos)

    # re-write dependent variables in term of independent variables
    Proton.sym_frac = electron.sym_frac
    Neutron.sym_frac = 1 - Proton.sym_frac 
    
    meson_coupling(eos, meson_list)

    sigma_coupling(eos)