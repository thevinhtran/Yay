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
    def __init__(self, g_sigma_N = 0.0, g_omega_N = 0.0, g_rho_N = 0.0, g_phi_N = 0.0, b = 0.0, c = 0.0,\
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
    # baryon particle 
    def __init__(self, mass, spin, isospin, charge, kind, var_type, mass_eff = 0.0, num_density = 0.0,\
                 frac = 0.0, kf = 0.0, ef = 0.0, chem_pot = 0.0):
    
        # variables to be established at baryon declaration
        self.mass = mass
        self.spin = spin 
        self.isospin = isospin
        self.charge = charge
        self.kind = kind
        self.var_type = var_type
    
        # variables to be stored later
        self.mass_eff = mass_eff
        self.num_density = num_density
        self.frac = frac
        self.kf = kf
        self.ef = ef
        self.chem_pot = chem_pot
        
        # coupling constants 
        self.g_sigma = 0.0
        self.g_omega = 0.0 
        self.g_rho = 0.0 
        self.g_phi = 0.0 


class lepton:
    # lepton particle 
    def __init__(self, mass, charge, num_density = 0.0, frac = 0.0, var_type = 0.0, kf = 0.0, chem_pot = 0.0):
        self.mass = mass
        self.charge = charge

        self.num_density = num_density
        self.frac = frac
        self.var_type = var_type
        self.kf = kf
        self.chem_pot = chem_pot


class meson:
    def __init__(self, mass, field = 0.0):
        self.mass = mass # in MeV
        self.field = field

        # coupling constants 
        # could update this in the future to store coupling constants here
        # but for now only makes sense for the sigma meson to store self coupling 
        self.b = 0.0 
        self.c = 0.0 

# Pre-initialize our objects so that at run time all we need to do is
# call the lists!! 
""" Baryons """

# symbolic baryon objects
proton_sym = baryon(sym.symbols('m_p'), 1/2, 1/2, 1, 'Nucleon', 'Dependent', sym.symbols('m_p^*'),\
                    sym.symbols('n_p'), sym.symbols('x_p'), sym.symbols('k_p'),\
                    sym.symbols('E^*_p'), sym.symbols('mu_p'))
neutron_sym = baryon(sym.symbols('m_n'), 1/2, -1/2, 0, 'Nucleon', 'Dependent', sym.symbols('m_n^*'),\
                    sym.symbols('n_n'), sym.symbols('x_n'), sym.symbols('k_n'),\
                    sym.symbols('E^*_n'), sym.symbols('mu_n'))

lambda_sym = baryon(sym.symbols('m_Lambda'), 1/2, 0, 0, 'Hyperon', 'Independent', sym.symbols('m_Lambda^*'),\
                    sym.symbols('n_Lambda'), sym.symbols('x_Lambda'), sym.symbols('k_Lambda'),\
                    sym.symbols('E^*_Lambda'), sym.symbols('mu_Lambda'))

sigma_neu_sym = baryon(sym.symbols('m_Sigma_0'), 1/2, 1.0, 0.0, 'Hyperon', '', sym.symbols('m_Sigma_0^*'),\
                    sym.symbols('n_Sigma_0'), sym.symbols('x_Sigma_0'), sym.symbols('k_Sigma'),\
                    sym.symbols('E^*_Sigma'), sym.symbols('mu_Sigma'))

sigma_plus_sym = baryon(sym.symbols('m_Sigma_+'), 1/2, 1.0, 1.0, 'Hyperon', '', sym.symbols('m_Sigma_+^*'),\
                    sym.symbols('n_Sigma_+'), sym.symbols('x_Sigma_+'), sym.symbols('k_Sigma_+'),\
                    sym.symbols('E^*_Sigma_+'), sym.symbols('mu_Sigma_+'))

sigma_min_sym = baryon(sym.symbols('m_Sigma_-'), 1/2, 1.0, -1.0, 'Hyperon', '', sym.symbols('m_Sigma_-^*'),\
                    sym.symbols('n_Sigma_-'), sym.symbols('x_Sigma_-'), sym.symbols('k_Sigma_-'),\
                    sym.symbols('E^*_Sigma_-'), sym.symbols('mu_Sigma_-'))

xi_neu_sym = baryon(sym.symbols('m_Xi_0'), 1/2, 1/2, 0, 'Hyperon', '', sym.symbols('m_Xi_0^*'),\
                    sym.symbols('n_Xi_0'), sym.symbols('x_Xi_0'), sym.symbols('k_Xi_0'),\
                    sym.symbols('E^*_Xi_0'), sym.symbols('mu_Xi_0'))

xi_min_sym = baryon(sym.symbols('m_Xi_-'), 1/2, 1/2, -1.0, 'Hyperon', '', sym.symbols('m_Xi_-^*'),\
                    sym.symbols('n_Xi_-'), sym.symbols('x_Xi_-'), sym.symbols('k_Xi_-'),\
                    sym.symbols('E^*_Xi_-'), sym.symbols('mu_Xi_-'))




# numeric baryon objects 
proton_num = baryon(939.0, 1/2, 1/2, 1, 'Nucleon', 'Dependent')
neutron_num = baryon(939.0, 1/2, -1/2, 0, 'Nucleon', 'Dependent')
lambda__num = baryon(1116.0, 1/2, 0, 0, 'Hyperon', 'Independent')

sigma_neu_num = baryon(1192.642, 1/2, 1, 0, 'Hyperon', 'Dependent')
sigma_min_num = baryon(1197.5, 1/2, 1, -1, 'Hyperon', 'Dependent')
sigma_plus_num = baryon(1189.37, 1/2, 1, 1, 'Hyperon', 'Dependent')

xi_neu_num = baryon(1314.86, 1/2, 1/2, 0, 'Hyperon', 'Dependent')
xi_min_num = baryon(1321.72, 1/2, 1/2, -1, 'Hyperon', 'Dependent')

""" Leptons """

# symbolic lepton objects 
electron_sym = lepton(sym.symbols('m_e'), -1, sym.symbols('n_e'), sym.symbols('x_e'), 'Independent',\
                      sym.symbols('k_F_e'), sym.symbols('\mu_e'))
muon_sym = lepton(sym.symbols('m_mu'), -1.0, sym.symbols('n_mu'), sym.symbols('x_\mu'), 'Independent',\
                    sym.symbols('k_F_mu'), sym.symbols('\mu_mu'))

# numeric lepton objects 
electron_num = lepton(0.510, -1.0)
muon_num = lepton(105.65, -1.0)




""" Mesons """

# symbolic meson objects 
sigma_sym = meson(sym.symbols('m_sigma'), sym.symbols('sigma'))
omega_sym = meson(sym.symbols('m_omega'), sym.symbols('omega'))
rho_sym = meson(sym.symbols('m_rho'), sym.symbols('rho'))
phi_sym = meson(sym.symbols('m_phi'), sym.symbols('phi'))

# numeric meson objects
sigma_num = meson(550.0)
omega_num = meson(783.0)
rho_num = meson(770.0)
phi_num = meson(1020.0)



#  Load in baryon objects their coupling constants
def baryon_coupling(baryon, eos):
    if (baryon.kind == 'Nucleon'):
        baryon.g_sigma, baryon.g_omega, baryon.g_rho, baryon.g_phi = eos.g_sigma_N, eos.g_omega_N, eos.g_rho_N, eos.g_phi_N
    elif (baryon.kind == 'Hyperon'):
        baryon.g_sigma = eos.g_sigma_H
        baryon.g_omega = eos.g_omega_H
        baryon.g_rho = eos.g_rho_H
        baryon.g_phi = eos.g_phi_H


# load in sigma field self coupling 
def sigma_coupling(eos):
    sigma_sym.b = eos.b
    sigma_sym.c = eos.c 

# write an initialization function 

# def initialization(eos, baryon_list):
    # load in baryons

    # for baryon in baryon_list:
    #    baryon_coupling(baryon, eos)

    # sigma_coupling(eos)

    # 


""" When initializing a system, would just need to
1. Initialize the independent variables. 
2. Declare the baryon,lepton, meson, and independent variable lists 
"""


""" Need to then make the functions to generate the equations of motion """

# First: generating the sigma equation of motion 

def scalar_density(baryon):
    # returns scalar density n_s for a given baryon
    # note: we modify the argument of the natural log as to be positive definite to avoid complex numbers 

    coeff_1 = (2*baryon.spin + 1)/(2*Pi**2)
    coeff_2 = baryon.g_sigma*baryon.mass_eff
    term_2 = baryon.ef*baryon.kf 
    term_2_2 = sym.sqrt(((baryon.ef + baryon.kf)/baryon.mass_eff)**2)
    term_3 = sym.log(term_2_2)
    
    return coeff_1*coeff_2*(term_2 - baryon.mass_eff**2*term_3)


def sigma_eom_init(baryon_list):
    # returns symbolic sigma equation of motion in a mostly simplified form 
    
    term_1 = sigma_sym.mass**2*sigma_sym.field 
    term_2 = sigma_sym.b * neutron_sym.mass * neutron_sym.g_sigma**3 * sigma_sym.field**2 
    term_3 = sigma_sym.c * neutron_sym.g_sigma**4 * sigma_sym.field**3
    
    tot = 0
    
    for baryon in baryon_list:
        tot += baryon.g_sigma * scalar_density(baryon)

    return term_1 + term_2 + term_3 - tot

def sigma_eom(baryon_list):
    # substitutes into fully symbolic sigma EOM for
    # effective mass and effective energy to get a function in terms of 
    # fermi momenta and sigma field 

    result = sigma_eom_init(baryon_list)

    for baryon in baryon_list:
        result = result.subs(baryon.ef, sym.sqrt(baryon.kf**2 + baryon.mass_eff**2))
        result = result.subs(baryon.mass_eff, baryon.mass - baryon.g_sigma * sigma_sym.field)
    
    return result 

# omega meson EOM

#def omega_eom(baryon_list):
#    result = 0 
#    for baryon in baryon_list:
#        result += baryon.g_omega*baryon.num_density
#    return sym.symbols('m_omega')**2*sym.symbols('omega') - result

#def omega_eom_expanded(baryon_list):
#    result = omega_eom(baryon_list) 
#    for baryon in baryon_list:
#        result = result.subs(baryon.num_density, sym.symbols('n_B')*baryon.frac)
#    return result

def omega_eom(baryon_list):
    # returns expression for omega equation of motion in terms of baryon fermi momentum
    # note this expression is equal to zero 
    result = 0 
    for baryon in baryon_list:
        result += baryon.g_omega * baryon.kf**3 
    return omega_sym.mass**2 * omega_sym.field * (3*Pi**2) - result 



# rho meson EOM 

#def rho_eom(baryon_list):
#    result = 0
#    for baryon in baryon_list:
#        result += baryon.g_rho * baryon.num_density * baryon.isospin 
#    return sym.symbols('m_rho')**2 * sym.symbols('rho') - result 

#def rho_eom_expanded(baryon_list):
#    result = rho_eom(baryon_list) 
#    for baryon in baryon_list:
#        result = result.subs(baryon.num_density, sym.symbols('n_B')*baryon.frac)
#    return result 

def rho_eom(baryon_list):
    result = 0 
    for baryon in baryon_list:
        result += baryon.g_rho * baryon.kf**3 * baryon.isospin
    return rho_sym.mass**2 * rho_sym.field * (3*Pi**2) - result 



# phi meson EOM 

#def phi_eom(baryon_list):
#    result = 0 
#    for baryon in baryon_list:
#        result += baryon.g_phi*baryon.num_density
#    return sym.symbols('m_phi')**2*sym.symbols('phi') - result

#def phi_eom_expanded(baryon_list):
#    result = phi_eom(baryon_list)
#    for baryon in baryon_list:
#        result = result.subs(baryon.num_density, sym.symbols('n_B')*baryon.frac)
#    return result 

def phi_eom(baryon_list):
    # returns eom for phi meson in terms of baryon momenta rather than
    # number density

    result = 0
    for baryon in baryon_list:
        result += baryon.g_phi * baryon.kf**3 
    return phi_sym.mass**2 * (3*Pi**2) * phi_sym.field - result 



# beta equilibrium 

def baryon_chem_pot(baryon):
    # returns baryon chemical potential in terms of expanded effective mass, 
    # ie, in terms of m - gsigma*sigma 
    expr = sym.sqrt(baryon.kf**2 + baryon.mass_eff**2) + baryon.g_omega*omega_sym.field\
            + baryon.g_rho*rho_sym.field*baryon.isospin + baryon.g_phi*phi_sym.field 
    return expr.subs(baryon.mass_eff, baryon.mass - baryon.g_sigma * sigma_sym.field)


def lepton_chem_pot(lepton):
    # returns chemical potential in terms of momenta symbol 
    # and mass symbol
    return sym.sqrt(lepton.kf**2 + lepton.mass**2)


def beta_equilibrium(baryon_list):
    # generates a list of all beta equilibrium conditions 
    neutron_chem_pot = baryon_chem_pot(neutron_sym)
    electron_chem_pot = lepton_chem_pot(electron_sym)

    equation_array =  [] 
    for baryon in baryon_list:
        if (baryon != neutron_sym):
            equation_array.append(baryon_chem_pot(baryon) - neutron_chem_pot + baryon.charge*electron_chem_pot)
                
    return equation_array


# charge conservation

def charge_conservation(baryon_list, lepton_list):
    # gives charge conservation equation condition on fermi momentum 
    particles = baryon_list + lepton_list 
    expression = 0 

    for particle in particles:
        if (particle.charge > 0):
            expression += particle.kf
        elif (particle.charge < 0):
            expression -= particle.kf 
    
    return expression 


# baryon number conservation

def baryon_num_conservation(baryon_list):
    # gives baryon number conservation equation condition on fermi momentum
    # of individual baryons rather than the individual species' number density 
    result = 0 
    for baryon in baryon_list:
        result += baryon.kf**3
    return 3*Pi**2*sym.symbols('n_B') 


""" System of Equations Generator 
        - Takes the above and generates our system of equations 
"""
def sys_eqn_gen(baryon_sym_list, lepton_sym_list):
    # function to generate all our equations and to store in an array
    # called sys_eqn_gen 
    func_gen = [sigma_eom, omega_eom, rho_eom, phi_eom, beta_equilibrium,\
                    charge_conservation, baryon_num_conservation]
    sys_eqn = []
    
    for function in func_gen:
        if (function == charge_conservation):
            # since charge conservation depends on both baryons and leptons
            sys_eqn.append(function(baryon_sym_list, lepton_sym_list))
        elif (function == beta_equilibrium):
             # beta condition function returns an array with multiple equations
            # we unpack that here
            beta_conditions = function(baryon_sym_list)
            for equation in beta_conditions:
                sys_eqn.append(equation)
        else:
            sys_eqn.append(function(baryon_sym_list))
            
    return sys_eqn


""" Substitution Function 
    - Up to this point all of our expressions have been symbolic which makes it easy
    - to verify their legitimateness. Now we want to perform numerical calculations
    - so we want to substitute in for the symbolic expressions for the masses and stuff 
"""
#arg_list = [baryon_sym_list, baryon_num_list, meson_sym_list, meson_num_list, lepton_sym_list, lepton_num_list]

def substitution(equation, baryon_sym_list, baryon_num_list,\
    meson_sym_list, meson_num_list, lepton_sym_list, lepton_num_list):

    # loops through baryons, mesons, leptons to replace masses and stuff 

    # baryons 
    for i in range(len(baryon_sym_list)):
        equation = equation.subs([(baryon_sym_list[i].mass, baryon_num_list[i].mass),\
                                 (sym.symbols('b'), sigma_num.b), (sym.symbols('c'), sigma_num.c),\
                                 (baryon_sym_list[i].g_sigma, baryon_num_list[i].g_sigma),\
                                 (baryon_sym_list[i].g_omega, baryon_num_list[i].g_omega),\
                                 (baryon_sym_list[i].g_rho, baryon_num_list[i].g_rho),\
                                 (baryon_sym_list[i].g_phi, baryon_num_list[i].g_phi)])
    
    # mesons 
    for i in range(len(meson_sym_list)):
        equation = equation.subs(meson_sym_list[i].mass, meson_num_list[i].mass)
    
    # leptons 
    for i in range(len(lepton_sym_list)):
        equation = equation.subs(lepton_sym_list[i].mass, lepton_num_list[i].mass)
    
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