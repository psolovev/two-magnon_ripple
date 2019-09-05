import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 14
plt.rcParams["font.family"] = "Times New Roman"

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


def find_cross(omega, omega_k, H):
    diff = np.abs(omega - omega_k)
    min_diff_ix = np.argmin(diff)
    toler = abs(omega[min_diff_ix] - omega[min_diff_ix-1])

    if min_diff_ix > 0 and diff[min_diff_ix] < toler:
        cross_y = omega_k[min_diff_ix]
        cross_x = H[min_diff_ix]
        return cross_x, cross_y
    else:
        return False, False


def plot_fmr(params):
    H = np.arange(1, 400, 1)
    omega = unif_fmr(H, params)

    # d = 40*10**-7 # cm
    plt.plot(H, omega * 10 ** -3, '-k', alpha=0.3, linewidth=4, label='Uniform FMR')

    d_vect = np.arange(10, 100, 10) * 10**-7

    for d in d_vect:
        omega_k = spin_wave_fmr(H, params, d)
        cross_x, cross_y = find_cross(omega, omega_k, H)
        plt.plot(H, omega_k * 10**-3, label='d = ' + str(int(np.ceil(d*10**7))) + ' nm')
        if cross_x:
            plt.plot(cross_x, cross_y * 10**-3, '*r', markersize=8)
    plt.xlabel('H (Oe)')
    plt.ylabel('Frequency (GHz)')

    legfont = {'fontsize': 10}
    plt.legend(loc='lower right', **legfont)
    plt.savefig('freq_vs_field.png', bbox_inches='tight', dpi=600)


    plt.show()


plot_fmr(params)

