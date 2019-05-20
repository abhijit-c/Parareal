import matplotlib.pyplot as plt
import numpy as np

T = 4
dt = .5
tt = np.transpose(np.arange(0, T+dt, dt))
true_sol = np.exp(-tt/2)

plt.figure()
for k in range(0, 4):
    data = np.loadtxt(f'Data/rk2_{k}.out')
    plt.semilogy(tt, np.abs(data[:]-true_sol), 's-', label=f'rk2, {k}-Parareal')
for k in range(0, 4):
    data = np.loadtxt(f'Data/rk4_{k}.out')
    plt.semilogy(tt, np.abs(data[:]-true_sol), 'x-', label=f'rk4, {k}-Parareal')
#plt.plot(tt, true_sol, 's-', label=f'True Solution')
plt.title('Runge Kutta Parareal Comparison')
plt.xlabel('$t \in [0,4]$ in increments of $1/2$')
plt.ylabel('$1$-norm error')
plt.legend(loc='lower left', bbox_to_anchor=(0, 0),
          ncol=1, fancybox=True, shadow=True)
plt.tight_layout()
plt.show()

plt.figure()

dt = .5
tt = np.transpose(np.arange(0, T+dt, dt))
true_sol = np.exp(-tt/2)
data = np.loadtxt(f'Data/rk2_{3}.out')
plt.semilogy(tt, np.abs(data[:]-true_sol), 's-', label=f'$\Delta t = 1/2$')

dt = .25
tt = np.transpose(np.arange(0, T+dt, dt))
true_sol = np.exp(-tt/2)
data = np.loadtxt(f'Data/rk2half_3.out')
plt.semilogy(tt, np.abs(data[:]-true_sol), 's-', label=f'$\Delta t = 1/4$')

dt = .25
tt = np.transpose(np.arange(0, T+dt, dt))
true_sol = np.exp(-tt/2)
data = np.loadtxt(f'Data/rk2half_3.out')
plt.semilogy(tt, (2**6)*np.abs(data[:]-true_sol), 's-', label=f'$mk(\Delta t = 1/4)$')

plt.title('RK2 3-Parareal Order Confirmation')
plt.xlabel('$t \in [0,4]$')
plt.ylabel('$1$-norm error')
plt.legend(loc='lower left', bbox_to_anchor=(0, 0),
          ncol=1, fancybox=True, shadow=True)
plt.tight_layout()
plt.show()
