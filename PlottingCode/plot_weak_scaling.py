import matplotlib.pyplot as plt
import numpy as np

plt.figure()

P = np.transpose(np.arange(1,26))
data = np.loadtxt(f'Data/weak_scaling.txt')
plt.plot(P, data, 's-', label='Computed Scaling')
plt.plot(P, data[0]*P, label='Linear Scaling')

plt.title('Weak Scaling Study for Parareal')
plt.xlabel('Number of Processors')
plt.ylabel('Time taken in seconds')
plt.legend(loc='upper left', bbox_to_anchor=(0, 1),
          ncol=1, fancybox=True, shadow=True)

plt.tight_layout()
plt.show()
