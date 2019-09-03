import numpy as np
import matplotlib.pyplot as plt

params = {'Ms': 955,
          'gamma': 2.95,
          'A': 10**-6}


def unif_fmr(H, params):
    gamma = params['gamma']
    Ms = params['Ms']
    pi_2 = 2*np.pi
    omega = (gamma) * (H * (H + 4*np.pi*Ms))**0.5
    return omega


def spin_wave_fmr(H, params, d):
    gamma = params['gamma']
    Ms = params['Ms']
    A = params['A']
    piMs_4 = 4*np.pi*Ms
    D = 2 * A / Ms
    k = (H / D) ** 0.5

    N = (1 - np.exp(-k*d)) / (k*d)

    one = H + D * k**2 + piMs_4 * N
    two = H + D * k**2
    omega_k = gamma * (one * two) ** 0.5
    return omega_k



def plot_fmr(params):
    H = np.arange(1, 400, 1)
    omega = unif_fmr(H, params)

    # d = 40*10**-7 # cm
    plt.plot(H, omega * 10 ** -3, '-k', alpha=0.3, linewidth=4, label='Uniform FMR')

    d_vect = np.arange(10, 100, 10) * 10**-7

    for d in d_vect:
        omega_k = spin_wave_fmr(H, params, d)
        plt.plot(H, omega_k * 10**-3, label='d = ' + str(int(np.ceil(d*10**7))) + ' nm')
    plt.xlabel('H (Oe)')
    plt.ylabel('Frequency (GHz)')
    plt.legend(loc='lower right')
    plt.show()


plot_fmr(params)

