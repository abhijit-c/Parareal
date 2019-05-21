import matplotlib.pyplot as plt
import numpy as np

TFINE = .4555

plt.figure()

P = np.transpose(np.arange(1,29))
data = np.loadtxt(f'Data/scalability_k4.txt')
plt.loglog(P, data, 's-', label='4-Parareal')
speedstr = f'Speedup: {TFINE/data[-1]:.2f}'
data = np.loadtxt(f'Data/scalability_k6.txt')
plt.loglog(P, data, 's-', label='6-Parareal')
plt.loglog(P, 1/P, '--', label='Theoretical Speedup')
plt.loglog(P, TFINE*np.ones(P.shape), '--', label='Serial Computation')

plt.title('Strong Scaling Study for Parareal w/ 28 Processors')
plt.xlabel('Number of Processors')
plt.ylabel('Time taken in seconds')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
plt.text(1.0,0.1, speedstr, fontsize=14, bbox=props)
plt.legend(loc='upper right', bbox_to_anchor=(1, 1),
          ncol=1, fancybox=True, shadow=True)

plt.tight_layout()
plt.show()
