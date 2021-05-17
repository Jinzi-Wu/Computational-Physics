import matplotlib.pyplot as plt
import numpy as np

import damper as sd

omega_n = 10


y0 = np.ones([2, 1])
y0[0, 0] = 0
order = 0

fig, ax = plt.subplots(6, 2)

for sigma in np.array([0.1, 1, 10]):
     A = np.array([0, 1, -omega_n ** 2, -2 * sigma * omega_n]).reshape([2, 2])
     temp = sd.backward_euler(y0, A, 0.001)
     x = temp[0]
     y = temp[1]
     x1 = temp[2]
     y1 = temp[3]
     

     ax[int(order), 0].set_title('v-x BE')
     ax[int(order), 0].plot(x, y, 'g-')
     plt.ylabel('v')
     plt.xlabel('x')


     ax[int(order), 1].set_title('x-t BE')
     ax[int(order), 1].plot(x1, y1, 'y-')
     plt.ylabel('x')
     plt.xlabel('t')
    
    
     temp = sd.forward_euler(y0, A, 0.01)
     x = temp[0]
     y = temp[1]
     x1 = temp[2]
     y1 = temp[3]
     

    
     ax[int(order) + 1, 0].set_title('v-x FE')
     ax[int(order) + 1, 0].plot(x, y, 'g-')
     plt.ylabel('v')
     plt.xlabel('x')


     ax[int(order) + 1, 1].set_title('x-t FE')
     ax[int(order) + 1, 1].plot(x1, y1, 'y-')
     plt.ylabel('x')
     plt.xlabel('t')
     
     order += 2







fig.suptitle('backward euler & forward euler: sigma=0.1, sigma=1, sigma=10 ')
fig.savefig('euler.png')
plt.show()

