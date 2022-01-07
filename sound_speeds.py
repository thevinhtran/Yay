""" Sound Speeds Helper Files """ 

import numpy as np 
import pandas as pd
import sympy as sym 

# coupling constant symbols 
Pi = sym.symbols('pi')
hc = 197.33
g_sigma_N_sym, g_omega_N_sym, g_rho_N_sym, g_phi_N_sym = sym.symbols('g_sigma_N, g_omega_N, g_rho_N, g_phi_N')
b_sym, c_sym = sym.symbols('b, c')
g_sigma_H_sym, g_omega_H_sym, g_rho_H_sym, g_phi_H_sym = sym.symbols('g_sigma_H, g_omega_H, g_rho_H, g_phi_H')

U = sym.Function('U')


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


# meson partial derivatives

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


class independent_var:
    def __init__(self, var, func, tilde_chem_pot = 0.0, tilde_chem_pot_val = 0.0, num_val = 0.0, total_deriv = 0.0):
        self.var = var # var is the symbol for the variable
        self.func = func 
        
        self.tilde_chem_pot = tilde_chem_pot
        self.tilde_chem_pot_val = tilde_chem_pot_val
        
        self.num_val = num_val
        
        # total derivative of fraction wrt to nB
        self.total_deriv = total_deriv


# Pre-initialize our objects so that at run time all we need to do is
# call the lists!! 
""" Baryons """

# symbolic baryon objects
Proton = baryon(name = 'Proton', spin = 1/2, isospin = 1/2, charge = 1, kind = 'Nucleon', var_type = 'Dependent',\
                sym_mass = sym.symbols('m_p'), sym_mass_eff = sym.symbols('m_p^*'), sym_density = sym.symbols('n_p'),\
                sym_frac = sym.symbols('x_p'), sym_kf = sym.symbols('k_p'), sym_ef = sym.symbols('E^*_p'), sym_chem_pot = sym.symbols('mu_p'),\
                num_mass = 939.0)

Neutron = baryon(name = 'Neutron', spin = 1/2, isospin = -1/2, charge = 0.0, kind = 'Nucleon', var_type = 'Dependent',\
                sym_mass = sym.symbols('m_n'), sym_mass_eff = sym.symbols('m_n^*'), sym_density = sym.symbols('n_n'),\
                sym_frac = sym.symbols('x_n'), sym_kf = sym.symbols('k_n'), sym_ef = sym.symbols('E^*_n'), sym_chem_pot = sym.symbols('mu_n'),\
                num_mass = 939.0)

Lambda = baryon(name = 'Lambda', spin = 1/2, isospin = 0, charge = 0, kind = 'Nucleon', var_type = 'Indepdent',\
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


# independent variable object
nb = independent_var(sym.symbols('n_B'), sym.Function('n_B'))
xe = independent_var(sym.symbols('x_e'), sym.Function('x_e'), sym.symbols('mu tilde_x_e'))
xl = independent_var(sym.symbols('x_Lambda'), sym.Function('x_Lambda'), sym.symbols('mu tilde_x_Lambda'))



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


# helper functions 
def baryon_chemical_potential(baryon, meson_list):
    meson_contribution = 0.0
    if (baryon.kind == 'Nucleon'):
        for meson in meson_list:
            if (meson == rho):
                meson_contribution += meson.num_g_N * baryon.isospin * meson.num_field 
            else: 
                meson_contribution += meson.num_g_N * meson.num_field 
    elif (baryon.kind == 'Hyperon'):
        for meson in meson_list:
            if (meson == rho):
                meson_contribution += meson.num_g_H * baryon.isospin * meson.num_field 
            else: 
                meson_contribution += meson.num_g_H * meson.num_field 
    
    return meson_contribution + np.sqrt(baryon.num_kf**2 + (baryon.num_mass - baryon.num_g_sigma * sigma.num_field)**2)


# chemical potential electron derivative 
def chem_pot_electron(ind_var):
    mu_e = sym.sqrt(electron.sym_kf**2 + electron.sym_mass**2)
    mu_e = mu_e.subs(electron.sym_kf, (3*Pi**2 * sym.symbols('n_B')* electron.sym_frac)**(1/3))
    return mu_e.diff(ind_var)

def chem_pot_lepton(lepton, ind_var):
    mu_lep = sym.sqrt(lepton.sym_kf**2 + lepton.sym_mass**2)
    mu_lep = mu_lep.subs(lepton.sym_kf, (3*Pi**2 * sym.symbols('n_B')* lepton.sym_frac)**(1/3))
    return mu_lep.diff(ind_var)

# partial baryon chemical potential
def partial_mu_R(baryon_sym, baryon_list, meson_list, ind_var):
    
    result = 0 
    if (baryon_sym.kind == 'Nucleon'):
        for meson in meson_list:
            result += meson.sym_g_N * meson.partial(sym.symbols('n_B'), baryon_list, ind_var)
    elif (baryon_sym.kind == 'Hyperon'): 
        for meson in meson_list:
            result += meson.sym_g_H * meson.partial(ind_var)
    
    return result 


# derivative of fermi momentum
def partial_fermi(baryon_sym, ind_var):
    kFi = (Pi**2 * sym.symbols('n_B') * baryon_sym.sym_frac)/baryon_sym.sym_kf**2
    return kFi.diff(ind_var)


# 
def alpha(baryon_sym):
    term1 = (sym.S(3)/2)*(1/Pi**2)*baryon_sym.sym_g_sigma * baryon_sym.sym_mass_eff**2 * sym.log((baryon_sym.sym_kf + baryon_sym.sym_ef)/baryon_sym.sym_mass_eff)
    term2 = (sym.S(1)/2)*baryon_sym.sym_kf * baryon_sym.sym_ef
    term3 = baryon_sym.sym_mass_eff**2 * baryon_sym.sym_kf/baryon_sym.sym_ef

    return 2*(term1 - baryon_sym.sym_g_sigma / Pi**2 * (term2 + term3))

def beta(baryon_sym):
    return 2*baryon_sym.sym_mass_eff * baryon_sym.sym_kf**2 / Pi**2 / baryon_sym.sym_ef

def partial_sigma(baryon_list, ind_var):
    numerator = 0
    denominator = sigma.sym_mass**2 + 2*sym.symbols('b')*Neutron.sym_g_sigma**3*sigma.sym_field + 3*sym.symbols('c')*Neutron.sym_g_sigma**4*sigma.sym_field**2
    for baryon in baryon_list:
        numerator += baryon.sym_g_sigma * beta(baryon) * partial_fermi(baryon, ind_var)
        denominator -= baryon.sym_g_sigma * alpha(baryon) 
    return numerator/denominator 

def partial_mu_prime(baryon, baryon_list, ind_var):
    denominator = sym.sqrt(baryon.sym_kf**2 + baryon.sym_mass_eff**2)
    numerator = baryon.sym_kf * partial_fermi(baryon, ind_var) - baryon.sym_g_sigma * baryon.sym_mass_eff * partial_sigma(baryon_list, ind_var)
    return numerator/denominator 

def chem_pot_part_deriv(baryon, baryon_list, meson_list, ind_var):
    return partial_mu_prime(baryon, baryon_list, ind_var) + partial_mu_R(baryon, baryon_list, meson_list, ind_var)


def chem_pot_electron_num(ind_var):
    # takes symbolic electron chemical potential partial derivative
    # and returns numeric expression using information stored in electron_sym and 
    # electron_num objects. Idea is that electron_sym contains symbols like m_e and electron_num contains
    # numerical values for those symbols like: m_e = 0.510 MeV
    
    # load in symbolic expression 
    
    symbolic = chem_pot_electron(ind_var)

    # replace symbolic variables using sympy subs method 
    symbolic = symbolic.subs([(Pi, np.pi), (electron.sym_mass, electron.num_mass),\
        (electron.sym_frac, electron.num_frac), (nb.var, nb.num_val)])
    
    return sym.simplify(symbolic)


def chem_pot_part_deriv_num(baryon, baryon_list, meson_list, lepton_list, x_j):
    # substitute in numerical values and get a numerical 
    # result for partial derivative of baryon wrt to independent variable 
    
    # call the symbolic expression 
    symbolic_part_deriv = chem_pot_part_deriv(baryon, baryon_list, meson_list, x_j)
    

    # replace Pi using subs 
    symbolic_part_deriv = symbolic_part_deriv.subs(Pi, np.pi)
    

    # baryon substitution 
    for baryon in baryon_list: 
        symbolic_part_deriv = symbolic_part_deriv.subs(baryon.sym_mass, baryon.num_mass)
        symbolic_part_deriv = symbolic_part_deriv.subs([(baryon.sym_g_sigma, baryon.num_g_sigma),\
                                (baryon.sym_g_omega, baryon.num_g_omega), (baryon.sym_g_phi, baryon.num_g_phi),\
                                (baryon.sym_g_rho, baryon.num_g_rho)])
        symbolic_part_deriv = symbolic_part_deriv.subs(baryon.sym_mass_eff, baryon.num_mass - baryon.num_g_omega * sigma.num_field)
        # fermi momentum 
        symbolic_part_deriv = symbolic_part_deriv.subs(baryon.sym_kf, baryon.num_kf)
        # effective energy
        symbolic_part_deriv = symbolic_part_deriv.subs(baryon.sym_ef, baryon.num_ef)
        # particle fraction
        symbolic_part_deriv = symbolic_part_deriv.subs(baryon.sym_frac, baryon.num_frac)


    # meson substitution 
    for meson in meson_list:
        symbolic_part_deriv = symbolic_part_deriv.subs(meson.sym_mass, meson.num_mass) 
        symbolic_part_deriv = symbolic_part_deriv.subs(meson.sym_field, meson.num_field)


    # lepton substitution 
    for lepton in lepton_list:
        symbolic_part_deriv = symbolic_part_deriv.subs(lepton.sym_frac, lepton.num_frac)
    
    # replace partial derivative of U self energy
    symbolic_part_deriv = symbolic_part_deriv.subs(sym.diff(U(sym.symbols('sigma')),sym.symbols('sigma'),sym.symbols('sigma')),\
                                                  2*sigma.num_b*Neutron.num_g_sigma**3*sigma.num_field + 3*sigma.num_c*Neutron.num_g_sigma**4*sigma.num_field**2)
    symbolic_part_deriv = symbolic_part_deriv.subs(nb.var, nb.num_val)
    
    return symbolic_part_deriv


def partial_mu_x_e_tilde(x_j, neutron_sym, proton_sym):
    # returns numerical values for the tilde chemical potentials
    result = -chem_pot_part_deriv_num(neutron_sym, x_j) + chem_pot_part_deriv_num(proton_sym, x_j)\
        + chem_pot_electron_num(x_j)
    return sym.N(result)


def partial_mu_x_l_tilde(x_j, neutron_sym, lambda_sym):
    # returns numerical values for the tilde chemical potentials
    result = - chem_pot_part_deriv_num(neutron_sym, x_j) + chem_pot_part_deriv_num(lambda_sym, x_j)
    return sym.N(result)


