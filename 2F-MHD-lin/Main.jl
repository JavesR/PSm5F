using PyCall
using BenchmarkTools

@pyimport numpy as np
@pyimport matplotlib.pyplot as plt

include("./FunEqm.jl")
include("./FunDev.jl")
include("./FunScm.jl")

ϵ = 0.25
β = 0.01
a_minor = 80.

n0 = 0.05
nL = 1.0
nr = 256
dn = (nL-n0)/nr
r = n0:dn:nL

dr = a_minor*dn
ar = a_minor*r


q = 1.05 .+ 2.0 .* r.^ 2
s = ar ./ q .* rdiff1(q,dr,nr+1)
ss = rdiff1(s,dr,nr+1)


qmin = findmin(q)[1]
qmax = findmax(q)[1]
lmx,lkm,lkn = mode_selection(qmax,qmin,2,1,2,1)
fkm,fkn,fkp = MakeFk(q,lkm,lkn,ar)



visvol = 1.e-4
vispsi = 1.e-4

dt = 0.001
tmax = 10

n0bc = zeros(ComplexF64,lmx)
nLbc = zeros(ComplexF64,lmx)


function flow(tmax,dt)

    # initialization
    phi = zeros(ComplexF64,(nr+1,lmx))
    psi = zeros(ComplexF64,(nr+1,lmx))

    phi[:,2:end] .= 1.e-5 * sin.(π*(0:nr)/nr) .+ 0im
    psi[:,2:end] .= 1.e-5 * sin.(π*(0:nr)/nr) .+ 0im

    vol =  laplace(phi,fkm,dr,ar)
    cur = -laplace(psi,fkm,dr,ar)


    nt = Int64(tmax/dt)

    timeloop(vol,psi,nt)

end


function timeloop(vol,psi,nt)

    t = 0.0

    print("t=$t \n")
    ft = [  real(reshape(vol,(:,1))) imag(reshape(vol,(:,1))) real(reshape(psi,(:,1))) imag(reshape(psi,(:,1))) ]
    np.savetxt("t$t.dat",ft,fmt="%+16.7e")

    for i in 1:nt

        # vol,psi = Euler(vol,psi)
        # vol,psi = MEuler(vol,psi)
        vol,psi= RK4(vol,psi)

        t = t+dt

        if i%1000 == 0
            print("t=$t \n")
            ft = [  real(reshape(vol,(:,1))) imag(reshape(vol,(:,1))) real(reshape(psi,(:,1))) imag(reshape(psi,(:,1))) ]
            np.savetxt("t$i.dat",ft,fmt="%+16.7e")
        end

    end


end



flow(tmax,dt)















































































































# mutable struct Emq
  
#     ar :: StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
#     q :: Vector{Float64}
#     s :: Vector{Float64}
#     ss :: Vector{Float64}

# end


# mutable struct Var

#     vol :: Matrix{ComplexF64}
#     phi :: Matrix{ComplexF64}
#     psi :: Matrix{ComplexF64}
#     cur :: Matrix{ComplexF64}

# end


# mutable struct DrVar

#     dr_vol :: Matrix{ComplexF64}
#     dr_phi :: Matrix{ComplexF64}
#     dr_psi :: Matrix{ComplexF64}
#     dr_cur :: Matrix{ComplexF64}

# end

# mutable struct DyVar

#     dy_vol :: Matrix{ComplexF64}
#     dy_phi :: Matrix{ComplexF64}
#     dy_psi :: Matrix{ComplexF64}
#     dy_cur :: Matrix{ComplexF64}

# end


