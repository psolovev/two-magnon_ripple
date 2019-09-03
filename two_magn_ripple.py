import numpy as np
import matplotlib.pyplot as plt

params = {'Ms': 955,
          'gamma': 2.95,
          'Hk': 171.7,
          'A': 10**-6}


def unif_fmr(H, params):
    gamma = params['gamma']
    Ms = params['Ms']
    pi_2 = 2*np.pi
    omega = (gamma) * (H * (H + 4*np.pi*Ms))**0.5
    return omega


def plot_ufmr(params):
    H = np.arange(0, 100, 0.5)
    omega = unif_fmr(H, params)
    plt.plot(H, omega * 10**-3)
    plt.show()


plot_ufmr(params)

