# -*- coding: utf-8 -*-

'''
1. Améliorer optimisation eps2, l
2. Changer distrib initiale
3. Changer le vecteur d'onde
4. Différents Omega
'''


# Python packages
import matplotlib.pyplot
import numpy
import os

# MRG packages
import _env
import preprocessing
import processing
import postprocessing
import alpha_compute
import pdb
#import solutions

import random


def your_optimization_procedure(domain_omega, spacestep, wavenumber, f, f_dir, f_neu, f_rob,
                                beta_pde, alpha_pde, alpha_dir, beta_neu, beta_rob, alpha_rob,
                                Alpha, mu, chi, V_obj, mu1, V_0):
    eps1 = 0.01
    eps2_0 = 30*spacestep*40
    eps2 = eps2_0*spacestep*40
    eps0 = 0.00001*spacestep*40
    eps3 = 0.1*spacestep*40
    """This function return the optimized density.

    Parameter:
        cf solvehelmholtz's remarks
        Alpha: complex, it corresponds to the absorbtion coefficient;
        mu: float, it is the initial step of the gradient's descent;
        V_obj: float, it characterizes the volume constraint on the density chi;
        mu1: float, it characterizes the importance of the volume constraint on
        the domain (not really important for our case, you can set it up to 0);
        V_0: float, volume constraint on the domain (you can it up to 1).
    """

    k = 0
    (M, N) = numpy.shape(domain_omega)
    numb_iter = 10
    energy = numpy.zeros((numb_iter+1, 1), dtype=numpy.float64)

    ##########################################################
    #   Différentes distributions initiales
    ##########################################################

    # indices = list(range(2*(len(x)-1)//10, 4*(len(x)-1)//10))

    # indices.extend(list(range(6*(len(x)-1)//10, 8*(len(x)-1)//10)))
    # # print(indices)

    # # budget : percentage of the border we can cover with liners
    # budget = len(indices)/len(x)

    # x_sub = [x[k] for k in indices]
    # y_sub = [y[k] for k in indices]

    # # -- define material density matrix
    # chi = preprocessing._set_chi(M, N, x_sub, y_sub)
    # chi = preprocessing.set2zero(chi, domain_omega)

    ##########################################################

    while k < numb_iter and mu > eps0:
        print('---- iteration number = ', k)
        # print('1. computing solution of Helmholtz problem, i.e., u')

        u = processing.solve_helmholtz(domain_omega, spacestep, wavenumber, f, f_dir, f_neu,
                                       f_rob, beta_pde, alpha_pde, alpha_dir, beta_neu, beta_rob, alpha_rob)

        # p = processing.solve_helmholtz(domain_omega, spacestep, wavenumber, -2 * numpy.conj(u), numpy.zeros(
        #     (M, N)), f_neu, f_rob, beta_pde, alpha_pde, alpha_dir, beta_neu, beta_rob, alpha_rob)

        # # ubis = u/numpy.max(numpy.abs(u))
        # postprocessing._plot_perso_solution(p, chi)

        # print('2. computing solution of adjoint problem, i.e., p')
        p = processing.solve_helmholtz(domain_omega, spacestep, wavenumber, -2 * numpy.conj(u), numpy.zeros(
            (M, N)), f_neu, f_rob, beta_pde, alpha_pde, alpha_dir, beta_neu, beta_rob, alpha_rob)

        # print('3. computing objective function, i.e., energy')

        ene = your_compute_objective_function(
            domain_omega, u, spacestep, mu1, V_0)

        # print('1st energy = ', ene)

        energy[k] = ene

        # print('4. computing parametric gradient')

        Jp = -numpy.real(Alpha*u*p)

        Jp[1:, :] = Jp[:-1, :]

        # postprocessing._plot_perso_solution(Jp, chi*0)
        while ene >= energy[k] and mu > eps0:
            l = 0
            # print('    a. computing gradient descent')
            chi_next = projector(chi-mu*Jp, l, domain_omega)

            # postprocessing._plot_perso_solution(chi_next, chi*0)

            # print('    b. computing projected gradient')

            int0 = integral(chi_next)
            eps2 = eps2_0

            while abs(int0-V_obj) >= eps1:
                # print(int0,V_obj)
                if int0 > V_obj:
                    l -= eps2
                else:
                    l += eps2
                # print(l)
                chi_next = projector(chi-mu*Jp, l, domain_omega)
                int0 = integral(chi_next)
                eps2 /= 2
                # print(V_obj, int, eps2, l)
                # postprocessing._plot_perso_solution(chi_next, chi*0)

            # postprocessing._plot_perso_solution(chi_next, chi*0)
            # print('    c. computing solution of Helmholtz problem, i.e., u')
            alpha_rob = Alpha * chi_next
            u = processing.solve_helmholtz(domain_omega, spacestep, wavenumber, f, f_dir, f_neu,
                                           f_rob, beta_pde, alpha_pde, alpha_dir, beta_neu, beta_rob, alpha_rob)

            # print('    d. computing objective function, i.e., energy (E)')
            ene = your_compute_objective_function(
                domain_omega, u, spacestep, mu1, V_0)

            # postprocessing._plot_perso_solution(u, chi*0)
            # print("Energy = ", ene)

            energy[k+1] = ene
            if ene < energy[k]:
                # The step is increased if the energy decreased
                mu = mu * (1+eps3)
            else:
                # The step is decreased is the energy increased
                mu = mu / 2
        chi = chi_next
        k += 1

    grad = 0
    # print('end. computing solution of Helmholtz problem, i.e., u')

    return chi, energy, u, grad


def your_compute_objective_function(domain_omega, u, spacestep, mu1, V_0):
    """
    This function compute the objective function:
    J(u,domain_omega)= \int_{domain_omega}||u||^2 + mu1*(Vol(domain_omega)-V_0)

    Parameter:
        domain_omega: Matrix (NxP), it defines the domain and the shape of the
        Robin frontier;
        u: Matrix (NxP), it is the solution of the Helmholtz problem, we are
        computing its energy;
        spacestep: float, it corresponds to the step used to solve the Helmholtz
        equation;
        mu1: float, it is the constant that defines the importance of the volume
        constraint;
        V_0: float, it is a reference volume.
    """

    J = 0
    for i in range(M-1):
        for j in range(N-1):
            moy = (u[i, j]+u[i, j+1]+u[i+1, j]+u[i+1, j+1])/4
            J += numpy.abs(moy)**2*spacestep**2
    J += mu1*(preprocessing.volume(domain_omega)-V_0)

    return J


def projector(chi, l, domain):
    n, m = numpy.shape(chi)
    new_chi = numpy.zeros((n, m), dtype='float')
    for i in range(n):
        for j in range(m):
            if preprocessing.is_on_boundary(domain[i, j]) == 'BOUNDARY':
                new_chi[i][j] = max(0, min(chi[i][j]+l, 1))
    return new_chi


def integral(chi):
    global spacestep
    res = 0
    n, m = numpy.shape(chi)
    for j in range(m-1):
        for i in range(n-1):
            res += chi[i][j]*spacestep
    return res


### Affichage de l'énergie en fct de omega

def valeur_finale_energie(energy):
    E_final=energy[0]
    k=1
    while energy[k]!=0 and k<len(energy):
        E_final=energy[k]
        k+=1
    return E_final

if __name__ == '__main__':
    # Plage de fréquences audibles f = 20Hz - 20kHz
    # omega = 60 - 60000
    # wavenumber = 0.18 - 176
    
    # HF : 80 - 176
    # MF : 10 - 100
    # BF : 0.18 - 20

    # Plage de fréquences audibles f = 20Hz - 200Hz
    # omega = 60 - 600
    # wavenumber = 0.18 - 18

    L_wavenumber= numpy.linspace(0.18,18,200)
    print(L_wavenumber)
    L_Energy=[]
    for wavenumber in L_wavenumber:
        # ----------------------------------------------------------------------
        # -- Fell free to modify the function call in this cell.
        # ----------------------------------------------------------------------
        # -- set parameters of the partial differential equation
        kx = -1.0
        ky = -1.0
        c = 340
        omega = wavenumber*c

        # -- set parameters of the geometry
        # N = max(int(7.*wavenumber),20) # number of points along x-axis
        N = 100
        M = 2 * N  # number of points along y-axis
        level = 2  # level of the fractal : limited by N
        spacestep = 1.0 / N  # mesh size

        # ----------------------------------------------------------------------
        # -- Do not modify this cell, these are the values that you will be assessed against.
        # ----------------------------------------------------------------------
        # --- set coefficients of the partial differential equation
        beta_pde, alpha_pde, alpha_dir, beta_neu, alpha_rob, beta_rob = preprocessing._set_coefficients_of_pde(
            M, N)

        # -- set right hand sides of the partial differential equation
        f, f_dir, f_neu, f_rob = preprocessing._set_rhs_of_pde(M, N)

        # -- set geometry of domain
        domain_omega, x, y, _, _ = preprocessing._set_geometry_of_domain(
            M, N, level)

        # ----------------------------------------------------------------------
        # -- Fell free to modify the function call in this cell.
        # ----------------------------------------------------------------------
        # -- define boundary conditions
        # planar wave defined on top
        f_dir[:, :] = 0.0
        f_dir[0, 0:N] = 1.0
        # spherical wave defined on top
        # f_dir[:, :] = 0.0
        # f_dir[0, int(N/2)] = 10.0

        # -- initialize
        alpha_rob[:, :] = - wavenumber * 1j

        # indices = []
        # for i in range(int(len(x)*2/5)):
        #     a = random.randint(0, len(x)-1)
        #     if a not in indices:
        #         indices.append(a)
        # indices.sort()

        # -- define subset of border on which we put the liner
        # modify this to change liners distribution

        indices = list(range(0*(len(x)-1)//10, 3*(len(x)-1)//10))

        # indices.extend(list(range(6*(len(x)-1)//10, 8*(len(x)-1)//10)))
        # print(indices)

        # budget : percentage of the border we can cover with liners
        budget = len(indices)/len(x)

        x_sub = [x[k] for k in indices]
        y_sub = [y[k] for k in indices]

        # -- define material density matrix
        chi = preprocessing._set_chi(M, N, x_sub, y_sub)
        chi = preprocessing.set2zero(chi, domain_omega)

        # -- define absorbing material
        al = alpha_compute.compute(wavenumber*340)
        Alpha = al[0] + 1.0j*al[1]

        # -- this is the function you have written during your project
        #import compute_alpha
        #Alpha = compute_alpha.compute_alpha(...)
        alpha_rob = Alpha * chi

        # -- set parameters for optimization
        S = 0  # surface of the fractal
        for i in range(0, M):
            for j in range(0, N):
                if domain_omega[i, j] == _env.NODE_ROBIN:
                    S += 1
        V_0 = 1  # initial volume of the domain
        V_obj = integral(chi) 
        # numpy.sum(numpy.sum(chi)) / S  # constraint on the density
        # V_obj = numpy.sum(numpy.sum(chi))
        mu = 0.5  # initial gradient step
        mu1 = 10**(-5)  # parameter of the volume functional

        # ----------------------------------------------------------------------
        # -- Do not modify this cell, these are the values that you will be assessed against.
        # ----------------------------------------------------------------------
        # -- compute finite difference solution
        u = processing.solve_helmholtz(domain_omega, spacestep, wavenumber, f, f_dir, f_neu, f_rob,
                                    beta_pde, alpha_pde, alpha_dir, beta_neu, beta_rob, alpha_rob)
        chi0 = chi.copy()
        u0 = u.copy()

        # ----------------------------------------------------------------------
        # -- Fell free to modify the function call in this cell.
        # ----------------------------------------------------------------------
        # -- compute optimization
        energy = numpy.zeros((100+1, 1), dtype=numpy.float64)
        chi, energy, u, grad = your_optimization_procedure(domain_omega, spacestep, wavenumber, f, f_dir, f_neu,
                                                        f_rob, beta_pde, alpha_pde, alpha_dir, beta_neu, beta_rob, alpha_rob, Alpha, mu, chi, V_obj, mu1, V_0)
        # chi, energy, u, grad = solutions.optimization_procedure(domain_omega, spacestep, wavenumber, f, f_dir, f_neu, f_rob,
        #                    beta_pde, alpha_pde, alpha_dir, beta_neu, beta_rob, alpha_rob,
        #                    Alpha, mu, chi, V_obj, mu1, V_0)
        # --- en of optimization
        
        L_Energy.append(valeur_finale_energie(energy))

        chin = chi.copy()
        un = u.copy()

        # -- plot chi, u, and energy
        # postprocessing._plot_uncontroled_solution(u0, chi0)
        # postprocessing._plot_controled_solution(un, chin)
        # err = un - u0
        # postprocessing._plot_error(err)
        # postprocessing._plot_energy_history(energy)

        # print('End.')
        print("WAVENUMBER ",wavenumber," - ENERGY ",valeur_finale_energie(energy))
    
    L_omega=[w*340 for w in L_wavenumber]
    matplotlib.pyplot.plot(L_omega,L_Energy)
    matplotlib.pyplot.xlabel('Omega')
    matplotlib.pyplot.ylabel('Energy')
    matplotlib.pyplot.xscale("log")
    matplotlib.pyplot.show()





        