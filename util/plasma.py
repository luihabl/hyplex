import pandas as pd
import numpy as np
from scipy.special import iv
from scipy.constants import m_e, e, epsilon_0, k, pi

m_xe = 2.17e-25
m_he = 6.67e-27
at_sccm = 4.477962e17

def voltage_self_bias(t_e_ev, v_rf, m_i, alpha):
    l1 = np.log(iv(0, v_rf / t_e_ev))
    l2 = 0.5 * np.log(m_i / (2 * pi * m_e))
    l3 = - np.log(alpha)
    return t_e_ev * (l1 + l2 + l3) 

def phi_const_rho(rho, x, boundaries=[0,0]):
    L = x.max()
    c0 = - rho / (2 * epsilon_0)
    c1 = rho * L / (2 * epsilon_0) - boundaries[0] / L + boundaries[1] / L
    c2 = boundaries[0]
    return c0 * x**2 + c1 * x + c2 

def debye_length(t_e_ev, n_0):
    return np.sqrt(epsilon_0 * e * t_e_ev / (e**2 * n_0))

def plasma_freq_e(n_0):
    return np.sqrt(n_0 * e**2 / (m_e * epsilon_0))

def v_bohm(t_e_ev, m_i):
    return np.sqrt(e * t_e_ev / m_i)

def v_volt(voltage, m):
    return np.sqrt(2 * e * voltage / m)

def mach_from_energy(energy, t_e_ev):
    return np.sqrt(2 * energy / t_e_ev)

def balanced_currents(n_0, t_e_ev, area, m_i, mach):
    j_i = e * n_0 * mach *v_bohm(t_e_ev, m_i)
    j_e = j_i * np.sqrt(m_i / (2 * pi * m_e)) / mach
    return (j_i* area, j_e * area)
    
def particle_per_cell(current, lx, ly, nx, ny, v, n_factor):
    return (current * np.sqrt(lx**2 + ly**2)) / (e * v * nx * ny * n_factor)
    
def particle_per_cell_2(current, dx, nx, ny, v, n_factor):
    return (current * dx) / (e * v * ny * n_factor)
    
    
def physical_space(i, a, L, N):
    return (a * L / (N - 1)**2) * i**2 + ((1 - a) * L / (N - 1)) * i
    
def logical_space(x, a, L, N):
    return ((N - 1) * 2 * x / L) / (1 - a + np.sqrt( (1 - a)**2 + (4*x*a/L) ))
    
