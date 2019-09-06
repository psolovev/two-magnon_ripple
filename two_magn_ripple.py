import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 14
plt.rcParams["font.family"] = "Times New Roman"
plt.rc('pdf', fonttype=42)

params = {'Ms': 955,
          'gamma': 2.95,
          'A': 10**-6}


def unif_fmr(H, params):
    gamma = params['gamma']
    Ms = params['Ms']
    pi_2 = 2*np.pi
    omega = (gamma) * (H * (H + 4*np.pi*Ms))**0.5
    return omega


def n_factor(k, d):
    # N = (1 - np.exp(-k*d)) / (1.0*(k*d))
    N = 1 - k*d / np.exp(1)
    return N


def n_factor_new(k, d):
    # N = (1 - np.exp(-k*d)) / (k*d)
    N = 1 - k*d / np.exp(1)
    return N


def plot_nfactor_vs_d():
    k = 40000
    d_vect = np.arange(20, 100, 0.1) * 10 ** -7
    N = n_factor(k, d_vect)
    N_new = n_factor_new(k, d_vect)
    plt.plot(d_vect, N, '-r')
    plt.plot(d_vect, N_new, '-b')
    plt.show()


def spin_wave_fmr(H, params, d):
    gamma = params['gamma']
    Ms = params['Ms']
    A = params['A']
    piMs_4 = 4*np.pi*Ms
    D = 2 * A / Ms
    k = (H / D) ** 0.5

    N = n_factor(k, d)

    one = H + D * k**2 + piMs_4 * N
    two = H + D * k**2
    omega_k = gamma * (one * two) ** 0.5
    return omega_k, N


def find_cross(omega, omega_k, H):

    omega_k_t = omega_k.copy()
    if np.any(np.isnan(omega_k_t)):
        notnan = omega_k_t[~np.isnan(omega_k_t)]
        omega_k_t[np.isnan(omega_k_t)] = np.min(notnan)

    diff = np.abs(omega - omega_k_t)
    min_diff_ix = np.argmin(diff)
    toler = abs(omega[min_diff_ix] - omega[min_diff_ix-1])

    if min_diff_ix > 0 and diff[min_diff_ix] < toler:
        cross_y = omega_k_t[min_diff_ix]
        cross_x = H[min_diff_ix]
        return cross_x, cross_y
    else:
        return False, False


def plot_fmr(params):
    H = np.arange(1, 500, 1)
    omega = unif_fmr(H, params)
    plt.plot(H, omega * 10 ** -3, '-k', alpha=0.3, linewidth=4, label='Uniform FMR')

    # d_vect = np.arange(10, 100, 10) * 10**-7
    d_vect = np.array([24, 32, 42, 56, 75, 100]) * 10**-7
    # d_vect = np.array([75]) * 10**-7
    N_vect = []
    for d in d_vect:
        omega_k, N = spin_wave_fmr(H, params, d)
        N_vect.append(N)
        cross_x, cross_y = find_cross(omega, omega_k, H)
        plt.plot(H, omega_k * 10**-3, label='d = ' + str(int(np.ceil(d*10**7))) + ' nm')
        if cross_x:
            plt.plot(cross_x, cross_y * 10**-3, '*r', markersize=8)
            print(cross_x, cross_y * 10**-3)
    plt.xlabel('H (Oe)')
    plt.ylabel('Frequency (GHz)')

    legfont = {'fontsize': 10}
    plt.legend(loc='lower right', **legfont)
    # plt.savefig('freq_vs_field.png', bbox_inches='tight', dpi=700)
    plt.savefig('freq_vs_field.pdf', bbox_inches='tight')
    plt.show()

    # plt.plot(d_vect, N_vect, '-o', label='Nfact vs d')
    # plt.legend()
    # plt.show()

plot_fmr(params)

# plot_nfactor_vs_d()
