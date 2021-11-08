
from scipy.optimize import minimize
import math as m
import numpy as np


################
# Caractéristiques de l'onde
################
A, B = 1, 1
k = 0
L = 10
omega = 300


################
# Caractéristiques du matériau
################

phi1 = 0.99  # porosity
alpha_h1 = 1.02  # tortuosity
sigma1 = 14000.0  # resitivity
gamma_p1 = 7.0 / 5.0
rho_01 = 1.2
c_01 = 340.0

###################
# Caractéristiques de l'air
###################

eta_0 = 1
ksi_0 = 1
a_0 = 0


ksi_1 = phi1*gamma_p1/(c_01**2)
a = sigma1*(phi1**2)*gamma_p1/((c_01)**2 * rho_01*alpha_h1)
eta_1 = phi1/alpha_h1
###################


def g(k):
    if k == 0:
        return 1
    else:
        return 0


def f(x, lambda_0):
    return ((lambda_0*eta_0-x)*np.exp(-lambda_0*L) + (lambda_0*eta_0+x)*np.exp(lambda_0*L))


def ki(k, alpha, lambda_0, lambda_1):
    return (g(k)*(((lambda_0*eta_0-lambda_1*eta_1)/f(lambda_1*eta_1, lambda_0)) - (lambda_0*eta_0-alpha)/f(alpha, lambda_0)))


def gamma(k, alpha, lambda_0, lambda_1):
    return (g(k)*(((lambda_0*eta_0+lambda_1*eta_1)/f(lambda_1*eta_1, lambda_0)) - (lambda_0*eta_0+alpha)/f(alpha, lambda_0)))


def e(alpha, k, lambda_0, lambda_1):
    fact = (A+B*(np.abs(k)**2))
    ki0 = ki(k, alpha, lambda_0, lambda_1)
    gamma0 = gamma(k, alpha, lambda_0, lambda_1)
    if k**2 >= ksi_0*(omega**2)/eta_0:
        term3 = 0
        term1 = (1/(2*lambda_0))*((np.abs(ki0)**2)*(1-np.exp(-2*lambda_0*L)) +
                                  np.abs(gamma0**2)*(np.exp(2*lambda_0*L)-1))+2*L*np.real(ki0*np.conj(gamma0))
        term2 = (B*(lambda_0**2)/2)*((np.abs(ki0)**2)*(1-np.exp(-2*lambda_0*L)) +
                                     np.abs(gamma0**2)*(np.exp(2*lambda_0*L)-1)) - 2*B*(lambda_0**2)*L*np.real(ki0*np.conj(gamma0))
    else:
        term1 = L*(np.abs(ki0)**2 + np.abs(gamma0**2)) + complex(0, 1) * \
            (1/lambda_0) * np.imag(ki0*np.conj(gamma0)*(1-np.exp(-2*lambda_0*L)))
        term2 = B*L*(np.abs(lambda_0)**2) * \
            (np.abs(ki0)**2 + np.abs(gamma0**2))
        term3 = complex(0, 1) * B*lambda_0*np.imag(ki0 *
                                                   np.conj(gamma0)*(1-np.exp(-2*lambda_0*L)))
    return (fact*term1 + term2 + term3)


def somme(alphaT):
    alpha = complex(1, 0)*alphaT[0] + complex(0, 1)*alphaT[1]
    s = 0
    for i in range(-10, 11):
        k = i*np.pi/L

        if k**2 > ksi_0*(omega**2)/eta_0:
            lambda_0 = np.sqrt(k**2-ksi_0*(omega**2)/eta_0)
        else:
            lambda_0 = complex(0, np.sqrt((ksi_0*(omega**2)/eta_0) - k**2))
        lambda_1 = complex(np.sqrt((k**2-ksi_1*(omega**2)/eta_1) + np.sqrt(((k**2-ksi_1*(omega**2)/eta_1)**2)+(a*omega/eta_1)**2)) /
                           np.sqrt(2), np.sqrt((-k**2+ksi_1*(omega**2)/eta_1) + np.sqrt(((k**2-ksi_1*(omega**2)/eta_1)**2)+(a*omega/eta_1)**2))/np.sqrt(2))

        s += e(alpha, k, lambda_0, lambda_1)

    return s


def compute(w):
    global omega
    omega = w
    alpha0 = [0, 0]
    alpha_min = minimize(somme, alpha0).x

    return alpha_min