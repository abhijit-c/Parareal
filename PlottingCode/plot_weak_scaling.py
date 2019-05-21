import matplotlib.pyplot as plt
import numpy as np

plt.figure()

P = np.transpose(np.arange(2,21))
data = np.loadtxt(f'Data/pipeline_weak.txt')
plt.loglog(P, data[:,0], 's-', label='Computed Scaling')
plt.loglog(P, data[:,1], 's-', label='Pipelined Computed Scaling')
plt.loglog(P, data[0,0]*P, 's-', label='Linear Scaling')

plt.title('Weak Scaling Study for Parareal w/ 28 Processors')
plt.xlabel('Number of Processors')
plt.ylabel('Time taken in seconds')
plt.legend(loc='upper left', bbox_to_anchor=(0, 1),
          ncol=1, fancybox=True, shadow=True)

plt.tight_layout()
plt.show()
