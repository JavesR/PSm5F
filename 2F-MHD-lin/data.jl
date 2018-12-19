# using PyPlot
using PyCall

@pyimport numpy as np
@pyimport matplotlib.pyplot as plt

nr = 256
r = range(0.05, stop = 1.0, length = nr+1)

data = np.genfromtxt("../../2f_v/xy0200.dat")

phi = reshape(data[:,3]+im*data[:,4],(nr+1,:))
psi = reshape(data[:,5]+im*data[:,6],(nr+1,:))
vol = reshape(data[:,7]+im*data[:,8],(nr+1,:))
cur = reshape(data[:,9]+im*data[:,10],(nr+1,:))
dphi = reshape(data[:,17]+im*data[:,18],(nr+1,:))
# psi_rhs = reshape(data[:,3]+im*data[:,4],(nr+1,:))
# vol_rhs = reshape(data[:,5]+im*data[:,6],(nr+1,:))


plt.plot(r,real(cur[:,2]))
# plt.plot(r,real(cur[:,2]))


data = np.genfromtxt("./t200.dat")

phi = reshape(data[:,3]+im*data[:,4],(nr+1,:))
psi = reshape(data[:,5]+im*data[:,6],(nr+1,:))
vol = reshape(data[:,7]+im*data[:,8],(nr+1,:))
cur = reshape(data[:,9]+im*data[:,10],(nr+1,:))


# plt.plot(r,real(psi[:,2]),lw=1.0)
plt.plot(r,real(cur[:,2]),ls="--")

plt.grid()

plt.show()
