import numpy as np
import pandas as pd
import sympy as sym


# important constants
hc = 197.33

# common symbols
Pi = sym.symbols('pi')

# coupling constant symbols 
g_sigma_N_sym, g_omega_N_sym, g_rho_N_sym, g_phi_N_sym = sym.symbols('g_sigma_N, g_omega_N, g_rho_N, g_phi_N')
b_sym, c_sym = sym.symbols('b, c')
g_sigma_H_sym, g_omega_H_sym, g_rho_H_sym, g_phi_H_sym = sym.symbols('g_sigma_H, g_omega_H, g_rho_H, g_phi_H')

# Class declaration 
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
    def __init__(self, spin, isospin, charge, kind, var_type,\
                sym_mass, sym_mass_eff, sym_density, sym_frac, sym_kf, sym_ef, sym_chem_pot,\
                num_mass = 0.0, num_mass_eff = 0.0, num_density = 0.0, num_frac = 0.0, num_kf = 0.0, num_ef = 0.0, num_chem_pot = 0.0):

        # variables common to both classes
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
    def __init__(self, charge, var_type,\
                sym_mass, sym_density, sym_frac, sym_kf, sym_chem_pot,\
                num_mass = 0.0, num_density = 0.0, num_frac = 0.0, num_kf = 0.0, num_chem_pot = 0.0):
        
        self.charge = charge
        self.var_type = var_type


        self.sym_mass = sym_mass
        self.sym_density = sym_density
        self.sym_frac = sym_frac
        self.sym_kf = sym_kf
        self.sym_chem_pot = sym_chem_pot

        self.num_mass = num_mass
        self.num_density = num_density
        self.frac = num_frac
        self.kf = num_kf
        self.chem_pot = num_chem_pot


class meson:
    def __init__(self, sym_mass, sym_field, num_mass, num_field = 0.0):
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

# Pre-initialize our objects so that at run time all we need to do is
# call the lists!! 
""" Baryons """

# symbolic baryon objects
Proton = baryon(spin = 1/2, isospin = 1/2, charge = 1, kind = 'Nucleon', var_type = 'Dependent',\
                sym_mass = sym.symbols('m_p'), sym_mass_eff = sym.symbols('m_p^*'), sym_density = sym.symbols('n_p'),\
                sym_frac = sym.symbols('x_p'), sym_kf = sym.symbols('k_p'), sym_ef = sym.symbols('E^*_p'), sym_chem_pot = sym.symbols('mu_p'),\
                num_mass = 939.0)

Neutron = baryon(spin = 1/2, isospin = -1/2, charge = 0.0, kind = 'Nucleon', var_type = 'Dependent',\
                sym_mass = sym.symbols('m_n'), sym_mass_eff = sym.symbols('m_n^*'), sym_density = sym.symbols('n_n'),\
                sym_frac = sym.symbols('x_n'), sym_kf = sym.symbols('k_n'), sym_ef = sym.symbols('E^*_n'), sym_chem_pot = sym.symbols('mu_n'),\
                num_mass = 939.0)

Lambda = baryon(spin = 1/2, isospin = 0, charge = 0, kind = 'Nucleon', var_type = 'Indepdent',\
                sym_mass = sym.symbols('m_Lambda'), sym_mass_eff =  sym.symbols('m_Lambda^*'), sym_density = sym.symbols('n_Lambda'),\
                sym_frac =  sym.symbols('x_Lambda'), sym_kf = sym.symbols('k_Lambda'), sym_ef = sym.symbols('E^*_Lambda'), sym_chem_pot = sym.symbols('mu_Lambda'),\
                num_mass = 1116.0)

Sigma_neu = baryon(spin = 1/2, isospin = 1.0, charge = 0.0, kind = 'Hyperon', var_type = '',\
                    sym_mass = sym.symbols('m_Sigma_0'), sym_mass_eff = sym.symbols('m_Sigma_0^*'),\
                    sym_density = sym.symbols('n_Sigma_0'), sym_frac = sym.symbols('x_Sigma_0'), sym_kf = sym.symbols('k_Sigma'),\
                    sym_ef = sym.symbols('E^*_Sigma'), sym_chem_pot = sym.symbols('mu_Sigma'),\
                    num_mass = 1192.642)

Sigma_plus = baryon(spin = 1/2, isospin = 1.0, charge = 1.0, kind = 'Hyperon', var_type = '',\
                    sym_mass = sym.symbols('m_Sigma_+'), sym_mass_eff = sym.symbols('m_Sigma_+^*'),\
                    sym_density = sym.symbols('n_Sigma_+'), sym_frac = sym.symbols('x_Sigma_+'), sym_kf = sym.symbols('k_Sigma_+'),\
                    sym_ef = sym.symbols('E^*_Sigma_+'), sym_chem_pot = sym.symbols('mu_Sigma_+'),\
                    num_mass = 1189.37)

Sigma_min = baryon(spin = 1/2, isospin = 1.0, charge = -1.0, kind = 'Hyperon', var_type = '',\
                    sym_mass = sym.symbols('m_Sigma_-'), sym_mass_eff = sym.symbols('m_Sigma_-^*'),\
                    sym_density = sym.symbols('n_Sigma_-'), sym_frac = sym.symbols('x_Sigma_-'), sym_kf = sym.symbols('k_Sigma_-'),\
                    sym_ef = sym.symbols('E^*_Sigma_-'), sym_chem_pot = sym.symbols('mu_Sigma_-'),\
                    num_mass = 1197.5)

xi_neu_sym = baryon(spin = 1/2, isospin = 1/2, charge = 0.0, kind = 'Hyperon', var_type = '',\
                    sym_mass = sym.symbols('m_Xi_0'), sym_mass_eff =  sym.symbols('m_Xi_0^*'),\
                    sym_density = sym.symbols('n_Xi_0'), sym_frac = sym.symbols('x_Xi_0'), sym_kf =  sym.symbols('k_Xi_0'),\
                    sym_ef = sym.symbols('E^*_Xi_0'), sym_chem_pot = sym.symbols('mu_Xi_0'),\
                    num_mass = 1314.86)

xi_min_sym = baryon(spin = 1/2, isospin = 1/2, charge = -1.0, kind = 'Hyperon', var_type = '',\
                    sym_mass = sym.symbols('m_Xi_-'), sym_mass_eff = sym.symbols('m_Xi_-^*'),\
                    sym_density = sym.symbols('n_Xi_-'), sym_frac = sym.symbols('x_Xi_-'), sym_kf = sym.symbols('k_Xi_-'),\
                    sym_ef = sym.symbols('E^*_Xi_-'), sym_chem_pot = sym.symbols('mu_Xi_-'),\
                    num_mass = 1321.72)


""" Leptons """

# symbolic lepton objects 
electron = lepton(charge = -1.0, var_type = 'Independent',\
                sym_mass = sym.symbols('m_e'), sym_density = sym.symbols('n_e'), sym_frac = sym.symbols('x_e'),\
                sym_kf = sym.symbols('k_e'), sym_chem_pot = sym.symbols('\mu_e'),\
                num_mass = 0.510)

muon = lepton(charge = -1.0, var_type = 'Dependent',\
            sym_mass = sym.symbols('m_mu'), sym_density = sym.symbols('n_mu'), sym_frac = sym.symbols('x_\mu'),\
            sym_kf = sym.symbols('k_mu'), sym_chem_pot =  sym.symbols('\mu_mu'),\
            num_mass = 105.65) 




""" Mesons """

# symbolic meson objects 
sigma = meson(sym_mass = sym.symbols('m_sigma'), sym_field = sym.symbols('sigma'), num_mass = 550.0)
omega = meson(sym_mass = sym.symbols('m_omega'), sym_field = sym.symbols('omega'), num_mass = 783.0)
rho = meson(sym_mass = sym.symbols('m_rho'), sym_field = sym.symbols('rho'), num_mass = 770.0)
phi = meson(sym_mass = sym.symbols('m_phi'), sym_field = sym.symbols('phi'), num_mass = 1020.0)



#  Load in baryon objects their coupling constants
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

# write an initialization function 

def init(eos, baryon_list):
    # load in baryons

    # load coupling constants into baryon objects
    for baryon in baryon_list:
        baryon_coupling(baryon, eos)

    # load in sigma self coupling 
    sigma_coupling(eos)



""" When initializing a system, would just need to
1. Initialize the independent variables. 
2. Declare the baryon,lepton, meson, and independent variable lists 
"""


""" Need to then make the functions to generate the equations of motion """

# First: generating the sigma equation of motion 

def scalar_density(baryon):
    # returns scalar density n_s for a given baryon
    # note: we modify the argument of the natural log as to be positive definite to avoid complex numbers 
    # this is a symbolic expression 

    coeff_1 = (2*baryon.spin + 1)/(2*Pi**2)
    coeff_2 = baryon.sym_g_sigma*baryon.sym_mass_eff
    term_2 = baryon.sym_ef*baryon.sym_kf 
    term_2_2 = sym.sqrt(((baryon.sym_ef + baryon.sym_kf)/baryon.sym_mass_eff)**2)
    term_3 = sym.log(term_2_2)
    
    return coeff_1*coeff_2*(term_2 - baryon.sym_mass_eff**2*term_3)


def sigma_eom_init(baryon_list):
    # returns symbolic sigma equation of motion in a mostly simplified form 
    
    term_1 = sigma.sym_mass**2*sigma.sym_field 
    term_2 = sigma.sym_b * Neutron.sym_mass * Neutron.sym_g_sigma**3 * sigma.sym_field**2 
    term_3 = sigma.sym_c * Neutron.sym_g_sigma**4 * sigma.sym_field**3
    
    tot = 0
    
    for baryon in baryon_list:
        tot += baryon.sym_g_sigma * scalar_density(baryon)

    return term_1 + term_2 + term_3 - tot

def sigma_eom(baryon_list):
    # substitutes into fully symbolic sigma EOM for
    # effective mass and effective energy to get a function in terms of 
    # fermi momenta and sigma field 

    result = sigma_eom_init(baryon_list)

    for baryon in baryon_list:
        result = result.subs(baryon.sym_ef, sym.sqrt(baryon.sym_kf**2 + baryon.sym_mass_eff**2))
        result = result.subs(baryon.sym_mass_eff, baryon.sym_mass - baryon.sym_g_sigma * sigma.sym_field)
    
    return result 

# omega meson EOM

def omega_eom(baryon_list):
    # returns expression for omega equation of motion in terms of baryon fermi momentum
    # note this expression is equal to zero 
    result = 0 
    for baryon in baryon_list:
        result += baryon.sym_g_omega * baryon.sym_kf**3 
    return omega.sym_mass**2 * omega.sym_field * (3*Pi**2) - result 


# rho meson EOM

def rho_eom(baryon_list):
    result = 0 
    for baryon in baryon_list:
        result += baryon.sym_g_rho * baryon.sym_kf**3 * baryon.isospin
    return rho.sym_mass**2 * rho.sym_field * (3*Pi**2) - result 



# phi meson EOM 

def phi_eom(baryon_list):
    # returns eom for phi meson in terms of baryon momenta rather than
    # number density

    result = 0
    for baryon in baryon_list:
        result += baryon.sym_g_phi * baryon.sym_kf**3 
    return phi.sym_mass**2 * (3*Pi**2) * phi.sym_field - result 



# beta equilibrium 

def baryon_chem_pot(baryon):
    # returns baryon chemical potential in terms of expanded effective mass, 
    # ie, in terms of m - gsigma*sigma 
    expr = sym.sqrt(baryon.sym_kf**2 + baryon.sym_mass_eff**2) + baryon.sym_g_omega*omega.sym_field\
            + baryon.sym_g_rho*rho.sym_field*baryon.isospin + baryon.sym_g_phi*phi.sym_field 
    return expr.subs(baryon.sym_mass_eff, baryon.sym_mass - baryon.sym_g_sigma * sigma.sym_field)


def lepton_chem_pot(lepton):
    # returns chemical potential in terms of momenta symbol 
    # and mass symbol
    return sym.sqrt(lepton.sym_kf**2 + lepton.sym_mass**2)


def beta_equilibrium(baryon_list):
    # generates a list of all beta equilibrium conditions 
    neutron_chem_pot = baryon_chem_pot(Neutron)
    electron_chem_pot = lepton_chem_pot(electron)

    equation_array =  [] 
    for baryon in baryon_list:
        if (baryon != Neutron):
            equation_array.append(baryon_chem_pot(baryon) - neutron_chem_pot + baryon.charge*electron_chem_pot)
                
    return equation_array


# charge conservation

def charge_conservation(baryon_list, lepton_list):
    # gives charge conservation equation condition on fermi momentum 
    particles = baryon_list + lepton_list 
    expression = 0 

    for particle in particles:
        if (particle.charge > 0):
            expression += particle.sym_kf
        elif (particle.charge < 0):
            expression -= particle.sym_kf 
    
    return expression 


# baryon number conservation

def baryon_num_conservation(baryon_list):
    # gives baryon number conservation equation condition on fermi momentum
    # of individual baryons rather than the individual species' number density 
    result = 0 
    for baryon in baryon_list:
        result += baryon.sym_kf**3
    return 3*Pi**2*sym.symbols('n_B') 


""" System of Equations Generator 
        - Takes the above and generates our system of equations 
"""
def sys_eqn_gen(baryon_list, lepton_list):
    # function to generate all our equations and to store in an array
    # called sys_eqn_gen 
    func_gen = [sigma_eom, omega_eom, rho_eom, phi_eom, beta_equilibrium,\
                    charge_conservation, baryon_num_conservation]
    sys_eqn = []
    
    for function in func_gen:
        if (function == charge_conservation):
            # since charge conservation depends on both baryons and leptons, we need to 
            # pass to it both baryon and lepton lists 
            sys_eqn.append(function(baryon_list, lepton_list))
        elif (function == beta_equilibrium):
            # beta condition function returns an array with (possibly) multiple equations
            # we unload those functions and append to array 
            beta_conditions = function(baryon_list)
            for equation in beta_conditions:
                sys_eqn.append(equation)
        else:
            sys_eqn.append(function(baryon_list))
            
    return sys_eqn


""" Substitution Function 
    - Up to this point all of our expressions have been symbolic which makes it easy
    - to verify their legitimateness. Now we want to perform numerical calculations
    - so we want to substitute in for the symbolic expressions for the masses and stuff 
"""
#arg_list = [baryon_list, baryon_list, meson_sym_list, meson_num_list, lepton_sym_list, lepton_num_list]

def substitution(equation, baryon_list, meson_list, lepton_list):

    # loops through baryons, mesons, leptons to replace masses and stuff with numeric values

    # baryons 
    for i in range(len(baryon_list)):
        equation = equation.subs([(baryon_list[i].sym_mass, baryon_list[i].num_mass),\
                                 (sigma.sym_b, sigma.num_b), (sigma.sym_c, sigma.num_c),\
                                 (baryon_list[i].sym_g_sigma, baryon_list[i].num_g_sigma),\
                                 (baryon_list[i].sym_g_omega, baryon_list[i].num_g_omega),\
                                 (baryon_list[i].sym_g_rho, baryon_list[i].num_g_rho),\
                                 (baryon_list[i].sym_g_phi, baryon_list[i].num_g_phi)])
    
    # mesons 
    for i in range(len(meson_list)):
        equation = equation.subs(meson_list[i].sym_mass, meson_list[i].num_mass)
    
    # leptons 
    for i in range(len(lepton_list)):
        equation = equation.subs(lepton_list[i].sym_mass, lepton_list[i].num_mass)
    
    equation = equation.subs(Pi, np.pi)

    return equation 


""" Fraction finder """

def fraction(fermi, nB):
    # after having solved for fermi momentum, can 
    # get the corresponding particle fraction 
    return fermi**3/3/np.pi**2/nB



""" Would be good to then include things here that do the solving
    for us 
"""