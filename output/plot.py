import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.constants import m_e, e, epsilon_0
m_i = 6.7e-27


# p_e = pd.read_csv('p_e.csv', header=None).values
# p_i = pd.read_csv('p_i.csv', header=None).values

# mesh = pd.read_csv('mesh.csv', header=None).values.ravel()
# phi = np.mean(pd.read_csv('phi.csv', header=None).values, axis = 1)
# efield = np.mean(pd.read_csv('efield.csv', header=None).values, axis = 1)
dens_i = np.mean(pd.read_csv('dens_i.csv', header=None).values, axis=1)
dens_e = np.mean(pd.read_csv('dens_e.csv', header=None).values, axis=1)

benchmark_data = pd.read_csv('turner_benchmark_results.dat', sep=' ', header=None, comment='#') 
caseA = benchmark_data[0:129]
caseB = benchmark_data[129:130 + 256]
# caseC = benchmark_data[386:387 + 512]

mesh = (6.7e-2) * np.arange(129) / 129


# rho = - 2.56e14 * e

# efield_expected = (rho * mesh / epsilon_0) - (0.5*rho * (6.7e-2) / epsilon_0)
# phi_expected = - ( rho * mesh**2 / (2 * epsilon_0) ) + (rho *  (6.7e-2) / (2 * epsilon_0)) * mesh

# plt.plot(mesh, efield, label='code', color='r')
# plt.plot(mesh, efield_expected, label='theory', color='b')

plt.plot(mesh, dens_i, label='$n_i$ code', color='b')
plt.plot(mesh, dens_e, label='$n_e$ code', color='r')

plt.plot(mesh, caseA[4], label='$n_i$ benchmark',  linestyle='--')
plt.plot(mesh, caseA[1], label='$n_e$ benchmark',  linestyle='--')


#plt.plot(mesh, efield, label='code', color='r')
#plt.plot(mesh, efield_expected, label='theory', color='b')

#plt.plot(caseC[0], caseC[4], label='$n_i$ benchmark', linestyle='--', color='b')
#plt.plot(caseC[0], caseC[1], label='$n_e$ benchmark', linestyle='--', color='r')

plt.xlabel('x [m]')
# plt.ylabel('$E$ [$V/m$]')
plt.legend()
plt.show()




#plt.plot(mesh, phi, label='code', color='r')
#plt.plot(mesh, phi_expected, label='theory', color='b')


plt.show()








# plt.plot(mesh, phi)
# plt.xlabel('x [m]')
# plt.ylabel('$\phi$ [V]')
# plt.show()
# 
# plt.plot(mesh, efield)
# plt.xlabel('x [m]')
# plt.ylabel('$E$ [V/m]')
# plt.show()

# plt.scatter(p_e[:, 0], p_e[:, 3], s=0.25)
# plt.xlabel('x [m]')
# plt.ylabel('$v_{e,x}$ [m/s]')
# plt.show()
# 
# plt.scatter(p_i[:, 0], p_i[:, 3], s=0.25)
# plt.xlabel('x [m]')
# plt.ylabel('$v_{i,x}$ [m/s]')
# plt.show()
# 
# plt.scatter(p_e[:, 0], 0.5 * m_e * (p_e[:, 3]**2 + p_e[:, 4]**2 + p_e[:, 5]**2) / e, s=0.25)
# plt.xlabel('x [m]')
# plt.ylabel('$E_e$ [eV]')
# plt.show()
# 
# 
# plt.scatter(p_i[:, 0], 0.5 * m_i * (p_i[:, 3]**2 + p_i[:, 4]**2 + p_i[:, 5]**2) / e, s=0.25)
# plt.xlabel('x [m]')
# plt.ylabel('$E_i$ [eV]')
# plt.show()





### Linear rho: (theory)
# def phi_lin(rho0, rho1, V0, V1, x):
#     L = x[-1]
#     a = (rho1 - rho0) / L
#     b = rho0
#     c1 = V0
#     c2 = (V1 - V0) / L + (a * L**2) / (6 * epsilon_0) + (b * L) / (2 * epsilon_0)
#     return - (a * np.power(x, 3)) / (6 * epsilon_0) - (b * np.power(x, 2)) / (2 * epsilon_0) + c2 * x + c1
#     
# def ef_lin(rho0, rho1, V0, V1, x):
#     L = x[-1]
#     a = (rho1 - rho0) / L
#     b = rho0
#     c2 = (V1 - V0) / L + (a * L**2) / (6 * epsilon_0) + (b * L) / (2 * epsilon_0)
#     return (a * np.power(x, 2)) / (2 * epsilon_0) + (b * x) / (epsilon_0) - c2 
#     
#     
# phi_lin_mesh = phi_lin(5e-6, -5e-6, -45, 45, mesh)
# e_lin_mesh = ef_lin(5e-6, -5e-6, -45, 45, mesh)
# 
# 
# 
# plt.plot(mesh, phi, label='C++')
# # plt.plot(mesh, phi_lin_mesh, label='Theory')
# plt.xlabel('x [m]')
# plt.ylabel('$\phi$ [V]')
# plt.legend()
# plt.show()
# 
# plt.plot(mesh, efield, label='C++')
# # plt.plot(mesh, e_lin_mesh, label='Theory')
# plt.xlabel('x [m]')
# plt.ylabel('E [V/m]')
# plt.legend()
# plt.show()







