import matplotlib.pyplot as plt
import numpy as np

T = 4
dt = .5
tt = np.transpose(np.arange(0, T+dt, dt))
true_sol = np.exp(tt)

plt.figure()
for k in range(0, 8):
    data = np.loadtxt(f'Data/linear1d_fw0{k}.out')
    plt.plot(tt, data[:], '--', label=f'{k}-Parareal')
plt.plot(tt, true_sol, 's-', label=f'True Solution')
plt.title('Linear 1d ODE Parareal Convergence')
plt.legend()
plt.show()

plt.figure()
for k in range(0, 8):
    data = np.loadtxt(f'Data/linear1d_fw0{k}.out')
    print(data - true_sol)
    plt.semilogy(tt, np.abs(true_sol-data[:]), '*-', label=f'{k}-Parareal')
plt.title('Linear 1d ODE Parareal Convergence')
plt.legend()
plt.show()
