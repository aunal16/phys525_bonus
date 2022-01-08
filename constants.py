from math import pi
n, p = 2.5 * 1e18, 2.5 * 1e18  # cm-3
E_g_J = 1.42 * 1.60218 * 1e-19  # eV * J/eV = J
E_g = 1.42  # eV
n_s = 3.6
T = 300  # K
tau_r = 2 * 1e-9  # s
m_e = 9.1094 * 1e-28  # g
m_e_eff = 0.067 * m_e
m_h_eff = 0.550 * m_e
m_reduc = (1/m_e_eff + 1/m_h_eff)**(-1)
kB = 1.3807 * 1e-16  # cm2 x g x s-2 x K-1
hbar = 1.0546 * 1e-27  # cm2 x g x s-1
h = hbar * 2 * pi
kBT_ev = (8.6173 * 1e-5) * T
c = 2.99792458 * 1e10  # cm s-1
