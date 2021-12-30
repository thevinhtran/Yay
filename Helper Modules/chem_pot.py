""" Collection of symbolic functions that together return a symbolic chemical potential
This uses sympy and the idea is to arrive first at a symbolic expression for the chemical potential after
which we can substitute in numerical values. 
"""

import numpy as np
import sympy as sym

Pi = sym.symbols('pi')

""" Electron/Lepton Chemical Potential Derivative """
def chem_pot_electron(ind_var):
    # calculates partial derivative of electron chemical potential
    # with respect to independent variable ind_var
    mu_e = sym.sqrt(sym.symbols('k_F_e')**2 + sym.symbols('m_e**2'))
    mu_e = mu_e.subs(sym.symbols('k_F_e')**2, (3*Pi**2*sym.symbols('n_B')*sym.symbols('x_e'))**(1/3))
    return mu_e.diff(ind_var)

""" "Trivial" Meson Field Partial Derivatives """

def partial_omega(eos, omega_mass, baryon_list, ind_var):
    # calculates partial derivative of omega or phi field wrt to independent variable
    # takes as argument equation of state, omega meson mass, list of baryons, 
    # and independent variable to which the partial derivative is taken 
    omega = 0
    for i in range(len(baryon_list)):
        if (baryon_list[i].kind == 'Nucleon'):
            omega = eos.g_omega_N*baryon_list[i].frac*sym.symbols('n_B') + omega
        elif (baryon_list[i].kind == 'Hyperon'):
            omega = eos.g_omega_H*baryon_list[i].frac*sym.symbols('n_B') + omega
    omega = 1/omega_mass**2*omega
    
    return sym.simplify(omega.diff(ind_var))


def partial_omega_phi(eos, phi_mass, baryon_list, ind_var):
    # calculates partial derivative of phi field wrt to independent variable
    # takes as argument equation of state, omega meson mass, list of baryons, 
    # and independent variable to which the partial derivative is taken 
    phi = 0
    for i in range(len(baryon_list)):
        if (baryon_list[i].kind == 'Nucleon'):
            phi = eos.g_phi_N*baryon_list[i].frac*sym.symbols('n_B') + phi
        elif (baryon_list[i].kind == 'Hyperon'):
            phi = eos.g_phi_H*baryon_list[i].frac*sym.symbols('n_B') + phi
    phi = 1/phi_mass**2*phi
    
    return sym.simplify(phi.diff(ind_var))


def partial_rho(eos, rho_mass, baryon_list, ind_var):
    # calculates partial derivative of omega or phi field wrt to independent variable
    # takes as argument equation of state, omega meson mass, list of baryons, 
    # and independent variable to which the partial derivative is taken 
    rho = 0
    for i in range(len(baryon_list)):
        if (baryon_list[i].kind == 'Nucleon'):
            rho = eos.g_rho*baryon_list[i].frac*sym.symbols('n_B')*baryon_list[i].isospin + rho
        elif (baryon_list[i].kind == 'Hyperon'):
            rho = eos.g_omega_H*baryon_list[i].frac*sym.symbols('n_B')*baryon_list[i].isospin + rho
    rho = 1/rho_mass**2*rho
    
    return sym.simplify(rho(ind_var))


def partial_mu_R(eos, baryon, baryon_list, omega_mass, phi_mass, rho_mass, ind_var):
    # returns dmu_i^R/dx_j
    
    if (baryon.kind == 'Nucleon'):
        return eos.g_omega_N*partial_omega_phi(eos, omega_mass, baryon_list, ind_var) + eos.g_phi_N*partial_omega_phi(eos, phi_mass, baryon_list, ind_var)\
            + baryon.isospin*eos.g_rho_N*partial_rho(eos, rho_mass, baryon_list, ind_var)
    
    elif (baryon.kind == 'Hyperon'):
        return eos.g_omega_H*partial_omega_phi(eos, omega_mass, baryon_list, ind_var) + eos.g_phi_H*partial_omega_phi(eos, phi_mass, baryon_list, ind_var)\
            + baryon.isospin*eos.g_rho_H*partial_rho(eos, rho_mass, baryon_list, ind_var)


''' Now, getting to functions that calculate partial derivative of the baryon chemical potentials...'''


def partial_fermi(baryon, ind_var):
    # calculates partial derivative of fermi momentum of input baryon with respect to ind_var
    # assumes input baryon number density has already been re-written in terms of independent variables nB, xe, xL

    k_fermi = (3*Pi**2*baryon.num_density)**(1/3)
    return k_fermi.diff(ind_var) 


def alpha(eos, baryon):
    # calculates alpha for input symbolic baryon. This alpha then is used to find the partial derivative
    # of sigma field with respect to ind_var 

    if (baryon.kind == 'Nucleon'):
        g_sigma = eos.g_sigma_N
    elif (baryon.kind == 'Hyperon'):
        g_sigma = eos.g_sigma_H
    
    term_1_1 = (3/2/Pi**2)*eos.g_sigma*baryon.mass_eff**2
    term_1_2 = sym.log((baryon.kf + baryon.ef)/baryon.mass_eff)

    term_2 = (1/2)*baryon.kf*baryon.ef
    term_3 = (baryon.mass_eff**2)*baryon.kf/baryon.ef

    return term_1_1*term_1_2 - eos.g_sigma/Pi**2*(term_2 + term_3)


def beta(baryon):
    # calculates beta for input symbolic baryon. This beta then is used to find the partial derivative 
    # of sigma field with respect to ind_var

    return baryon.mass_eff*baryon.kf**2/Pi**2/baryon.ef


def partial_sigma(eos, sigma_mass, baryon_list, ind_var):
    # calculates symbolic partial derivative of sigma field with respect to ind_var

    # second partial derivative of U(sigma) wrt sigma 
    U = sym.Function('U')
    sec_deriv_U = sym.diff(U(sym.symbols('sigma')),sym.symbols('sigma'),sym.symbols('sigma'))

    numerator = 0
    denominator = sigma_mass**2 + sec_deriv_U

    for i in range(len(baryon_list)):
        if (baryon_list[i].kind == 'Nucleon'):
            numerator = numerator + eos.g_sigma_N*beta(baryon_list[i])*partial_fermi(baryon_list[i], ind_var)
        elif (baryon_list[i].kind == 'Hyperon'):
            numerator = numerator + eos.g_sigma_H*beta(baryon_list[i])*partial_fermi(baryon_list[i], ind_var) 
    
    for i in range(len(baryon_list)):
        if (baryon_list[i].kind == 'Nucleon'):
            denominator = denominator - eos.g_sigma_N*alpha(eos, baryon_list)
        elif (baryon_list[i].kind == 'Hyperon'):
            denominator = denominator - eos.g_sigma_H*alpha(eos, baryon_list)
    
    return numerator/denominator 


def partial_mu_prime(eos, sigma_mass, baryon_list, baryon, ind_var):
    # calculates partial mu 

    if (baryon.kind == 'Nucleon'):
        g_sigma = eos.g_sigma_N 
    elif (baryon.kind == 'Hyperon'):
        g_sigma = eos.g_sigma_H 
    
    numerator = baryon.kf*partial_fermi(baryon, ind_var) - eos.g_sigma*baryon.mass_eff*partial_sigma(eos, sigma_mass, baryon_list, ind_var)
    denominator = sym.sqrt(baryon.kf**2 + baryon.mass_eff**2)
    return numerator/denominator 


def chem_pot_part_deriv(baryon, ind_var,\
    eos, baryon_list, meson_list):
    # adds partial_mu_prime and partial_mu_R together to get the full partial derivative of the chemical potential 
    # this is the function that returns the full chemical potential!!
    sigma_mass = meson_list[0].mass
    omega_mass = meson_list[1].mass
    rho_mass = meson_list[2].mass
    phi_mass = meson_list[3].mass
    return partial_mu_prime(eos, sigma_mass, baryon_list, baryon, ind_var) + partial_mu_R(eos, baryon, baryon_list, omega_mass, phi_mass, rho_mass, ind_var)


""" Baryon Chemical Potential """
def baryon_chemical_potential(eos, baryon, omega, phi, rho):
    # returns chemical potential for a baryon 
    # really only need this for the neutron chemical potential when calculating sound speeds 

    if (baryon.kind == 'Nucleon'):
        baryon.chem_pot = np.sqrt(baryon.kf**2 + baryon.ef**2) + eos.g_omega_N*omega.field\
                    + eos.g_phi_N*phi.field + baryon.isospin*eos.g_rho_N*rho.field
    
    elif (baryon.kind == 'Hyperon'):
        baryon.chem_pot = np.sqrt(baryon.kf**2 + baryon.ef**2) + eos.g_omega_H*omega.field\
                    + eos.g_phi_H*phi.field + baryon.isospin*eos.g_rho_H*rho.field