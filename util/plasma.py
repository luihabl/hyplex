import pandas as pd
import numpy as np
from scipy.special import iv
from scipy.constants import m_e, e, epsilon_0, k, pi

m_xe = 2.17e-25
m_he = 6.67e-27
at_sccm = 4.477962e17

def voltage_self_bias(t_ev, v_rf, m_i, alpha):
    l1 = np.log(iv(0, v_rf / t_ev))
    l2 = 0.5 * np.log(m_i / (2 * pi * m_e))
    l3 = - np.log(alpha)
    return t_ev * (l1 + l2 + l3) 

def debye_length(t_ev, n_0):
    return np.sqrt(epsilon_0 * e * t_ev / (e**2 * n_0))

