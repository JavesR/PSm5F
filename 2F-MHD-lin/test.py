import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("test.txt")

x = np.linspace(0,1,257)

plt.plot(x,data[0])
plt.plot(x,data[1])
plt.plot(x,data[2])


# data = np.genfromtxt("./2f_v/initial_profile.dat")

# ar= data[:,0]
# r = data[:,1]
# q = data[:,2]
# s = data[:,3]
# ss = data[:,4]

# plt.plot(r,q)
# plt.plot(r,s)
# plt.plot(r,ss)

plt.show()