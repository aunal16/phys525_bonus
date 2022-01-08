from scipy.optimize import fsolve   # numerical solver
from scipy.integrate import quad    # integration
from constants import *             # includes problem constants
import numpy as np
import matplotlib.pyplot as plt     # plotting
import warnings                     # to supress imaginary sqrt warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

INFTY = 1e2  # pseudo infinity value for the integral limits


def integral_func_n(x, x0):
    return np.sqrt(x)/(1 + np.exp(x - x0))


def integral_func_p(y, y0):
    return np.sqrt(y)/(1 + np.exp(y - y0))


def n_concentration_eq(x0):
    return -n + (2 * m_e_eff * kB * T)**(3/2) / (2 * np.pi**2 * hbar**3) * \
           (quad(integral_func_n, 0, INFTY, args=(x0,))[0])


def p_concentration_eq(x0):
    return -p + (2 * m_h_eff * kB * T)**(3/2) / (2 * np.pi**2 * hbar**3) * \
           (quad(integral_func_p, 0, INFTY, args=(x0,))[0])


def rough_plotter(x_start, x_end, func):
    x_line = np.linspace(x_start, x_end)
    y_line = np.zeros((len(x_line)))
    for i in range(len(x_line)):
        y_line[i] = func(x_line[i])
    plt.plot(x_line, y_line)
    plt.axhline(y=0, color='r', linestyle='--')
    plt.axvline(x=-1.4, color='g', linestyle='--')
    plt.show()


def main():
    """ This program plots and solves the equations specified in Sennaroglu Prob. 7.19"""
    print("PHYS 525 / ELEC 425: BONUS PROBLEM #1")
    print("-----------------------------------------")
    print("Implemented by Ahmet Hamdi Unal, 0060167.\nNo rights reserved, feel free to use it!")
    print("-----------------------------------------")
    print(" ----------- Program started! -----------\n")
    # 1ST PLOT: fc(E2) - fv(E1)
    # Plot for a small range of x0 for initial guess (Uncomment the next two lines to see)
    # rough_plotter(0, 10, e_concentration_eq)  # approx. 3.7
    # rough_plotter(-5, 10, p_concentration_eq) # approx. -1.4
    x0 = fsolve(n_concentration_eq, 3.7)[0]
    y0 = fsolve(p_concentration_eq, -1.4)[0]
    print("x0:", x0, "\ny0:", y0)

    # Efn = Ec + kBTx0
    # Efp = Ev - kBTy0
    # E2 = Ec + m_reduc/m_e_eff * (h\nu_0 - Eg)
    # E1 = Ev - m_reduc/m_h_eff * (h\nu_0 - Eg)
    hnu0 = np.linspace(1.4, 1.6, 2000)
    E2_minus_Efn = + m_reduc/m_e_eff * (hnu0 - E_g) - kBT_ev * x0
    E1_minus_Efp = - m_reduc/m_h_eff * (hnu0 - E_g) + kBT_ev * y0

    fcE2 = (np.exp(E2_minus_Efn / kBT_ev) + 1) ** (-1)
    fvE1 = (np.exp(E1_minus_Efp / kBT_ev) + 1) ** (-1)

    fig = plt.figure()
    fig.set_size_inches(15, 4.5, forward=True)
    plt.subplot(1, 3, 1)
    plt.plot(hnu0, fcE2 - fvE1)
    plt.title('Occupation Probabilities Difference')
    plt.xlabel(r'$E = h\nu_0$ (eV)', labelpad=-0)
    plt.ylabel(r'$f_c(E_2)-f_v(E_1)$', labelpad=0)
    plt.grid()

    # 2ND PLOT: sqrt(h\nu0 - Eg)
    # sqrt(h\nu0 - Eg)
    #plt.figure()
    plt.subplot(1, 3, 2)
    plt.plot(hnu0, np.sqrt(hnu0 - E_g))
    plt.xlim([1.4, 1.6])
    plt.title('√Photon Energy - Bandgap Difference')
    plt.xlabel(r'$E = h\nu_0$ (eV)', labelpad=-0)
    plt.ylabel(r'$\sqrt{E-E_g}\quad(\sqrt{eV})$', labelpad=0)
    plt.grid()

    # 3RD PLOT: gamma(\nu_0)
    # gamma(\nu_0)
    gamma = (fcE2 - fvE1) * (2 * m_reduc) ** (3/2) * c**2 / (2 * tau_r * n_s ** 4) * np.sqrt(hnu0 - E_g) / (hnu0)**2
    #plt.figure()
    plt.subplot(1, 3, 3)
    plt.plot(hnu0, gamma * 4.930993265537583 * 1e17)
    plt.xlim([1.4, 1.6])
    plt.title('Differential Gain')
    plt.xlabel(r'$E = h\nu_0$ (eV)', labelpad=-0)
    plt.ylabel(r'$\gamma(\nu_0)\quad(cm^{-1})$', labelpad=-7)
    plt.grid()

    # DETERMINE E value for fc(E2) - fv(E1)=0
    def zero_finder(E):
        E2_minus_Efn_0 = + m_reduc / m_e_eff * (E - E_g) - kBT_ev * x0
        E1_minus_Efp_0 = - m_reduc / m_h_eff * (E - E_g) + kBT_ev * y0
        fcE2_0 = (np.exp(E2_minus_Efn_0 / kBT_ev) + 1) ** (-1)
        fvE1_0 = (np.exp(E1_minus_Efp_0 / kBT_ev) + 1) ** (-1)

        return (fcE2_0 - fvE1_0)

    E_zero = fsolve(zero_finder, 1.475)[0]
    print("\nEnergy value that makes the occupational probabilities difference"+
          "0 is found to be", E_zero, "eV")

    # Show this value on all plots
    for i in range(3):
        plt.subplot(1, 3, i + 1)
        plt.axvline(x=E_zero, color='r', linestyle='--')
        xmin, xmax, ymin, ymax = plt.axis()
        plt.text(x=E_zero + (xmax - xmin)/20, y=ymin + (ymax - ymin)/10, s='E ≈ ' + str(np.round(E_zero, 4)) + ' eV',
                 color = 'w', rotation=90, bbox=dict(facecolor='red', alpha=0.6))

    plt.show()


if __name__ == "__main__":
    main()
