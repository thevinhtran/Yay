""" Functions to perform substitutions """

from chem_pot import * 

def chem_pot_electron_num(nb, electron_sym, electron_num, ind_var):
    # takes symbolic electron chemical potential partial derivative
    # and returns numeric expression using information stored in electron_sym and 
    # electron_num objects. Idea is that electron_sym contains symbols like m_e and electron_num contains
    # numerical values for those symbols like: m_e = 0.510 MeV
    
    # load in symbolic expression 
    symbolic = chem_pot_electron(ind_var)

    # replace symbolic variables using sympy subs method 
    symbolic = symbolic.subs([(Pi, np.pi), (electron_sym.mass, electron_num.mass),\
        (electron_sym.frac, electron_num.frac), (nb.var, nb.num_val)])
    
    return sym.simplify(symbolic)


def chem_pot_part_deriv_num(baryon, ind_var,\
    eos_sym, eos_num, baryon_list, baryon_num_list, meson_list, meson_num_list, lepton_list, lepton_num_list, nb):
    # calculates numerical partial derivative of input baryon chemical potential wrt to ind_var 

    # first calculate symbolic expression for partial derivative 
    symbolic = chem_pot_part_deriv(baryon, ind_var,\
        eos_sym, baryon_list, meson_list)

    # in this section: we replace symbols with their numerical values

    # replace Pi
    symbolic = symbolic.subs(Pi, np.pi)

    # replace coupling constants
    symbolic = symbolic.subs([(eos_sym.g_sigma_N, eos_num.g_sigma_N),\
        (eos_sym.g_sigma_H, eos_num.g_sigma_H), (eos_sym.g_omega_N, eos_num.g_omega_N),\
            (eos_sym.g_omega_H, eos_num.g_omega_H), (eos_sym.g_phi_N, eos_num.g_phi_N),\
                (eos_sym.g_phi_H, eos_num.g_phi_H), (eos_sym.g_rho_N, eos_num.g_rho_N),\
                    (eos_sym.g_rho_H, eos_num.g_rho_H)])

     # replace meson field masses and coupling constants
    for i in range(len(meson_list)):
        symbolic = symbolic.subs(meson_list[i].mass, meson_num_list[i].mass)

    # replace Baryon masses, baryon number density, effective energy, fermi momentum, fractions 
    for i in range(len(baryon_list)):
        symbolic = symbolic.subs(baryon_list[i].mass, baryon_num_list[i].mass)
        symbolic = symbolic.subs(baryon_list[i].kf, baryon_num_list[i].kf)
        symbolic = symbolic.subs(baryon_list[i].ef, baryon_num_list[i].ef)
        symbolic = symbolic.subs(baryon_list[i].num_density, baryon_num_list[i].num_density)
        symbolic = symbolic.subs(baryon_list[i].frac, baryon_num_list[i].frac)

    # replace lepton fractions
    for i in range(len(lepton_list)):
        symbolic  = symbolic.subs(lepton_list[i].frac, lepton_num_list[i].frac)
    
    # replace partial derivative of U self energy
    symbolic = symbolic.subs(sym.diff(U(sym.symbols('sigma')),sym.symbols('sigma'),sym.symbols('sigma')),\
                2*eos_num.b*eos_num.g_sigma_N**3*meson_num_list[0].field + 3*eos_num.c*eos_num.g_sigma_N**4*meson_num_list[0].field**2)
    
    # replace nB
    symbolic = symbolic.subs(nb.var, nb.num_val)

    return symbolic
                

""" Tilde Chemical Potentials """

def mu_xe_tilde_part_deriv_num(ind_var,\
    eos_sym, eos_num, baryon_list, baryon_num_list, meson_list, meson_num_list, lepton_list, lepton_num_list, nb, electron_sym, electron_num):
    # calculates partial derivative of mu_xe_tilde
    return chem_pot_part_deriv_num(baryon_list[1], ind_var, eos_sym, eos_num, baryon_list, baryon_num_list, meson_list, meson_num_list, lepton_list, lepton_num_list, nb)\
        - chem_pot_part_deriv_num(baryon_list[0], ind_var, eos_sym, eos_num, baryon_list, baryon_num_list, meson_list, meson_num_list, lepton_list, lepton_num_list, nb)\
             - chem_pot_electron_num(nb, electron_sym, electron_num, ind_var)


def mu_xl_tilde_part_deriv_num(ind_var,\
    eos_sym, eos_num, baryon_list, baryon_num_list, meson_list, meson_num_list, lepton_list, lepton_num_list, nb):
    # calculates partial derivative of mu_xe_tilde
    return chem_pot_part_deriv_num(baryon_list[1], ind_var, eos_sym, eos_num, baryon_list, baryon_num_list, meson_list, meson_num_list, lepton_list, lepton_num_list, nb)\
        - chem_pot_part_deriv_num(baryon_list[2], ind_var, eos_sym, eos_num, baryon_list, baryon_num_list, meson_list, meson_num_list, lepton_list, lepton_num_list, nb)


""" Systems of Equation Solver """

#def total_derivatives():
#    eqn_1 = equation_1


""" Sound Speed Difference Calculator """

#def sound_speed_diff():
    