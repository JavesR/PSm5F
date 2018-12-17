using PyCall
using BenchmarkTools

@pyimport numpy as np
@pyimport matplotlib.pyplot as plt

include("./FunEqm.jl")
include("./FunDerv.jl")

ϵ = 0.25
β = 0.001
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


n0bc = [0.0;0.0]
nLbc = [0.0;0.0]


visvol = 1.e-4
vispsi = 1.e-4

dt = 0.005
tmax = 500



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
    for i in 1:nt

        # vol,psi = Euler(vol,psi)
        vol,psi = RK4(vol,psi)

        t = t+dt

        if i%1000 == 0
            print("t=$t \n")
            ft = [  real(reshape(vol,(:,1))) imag(reshape(vol,(:,1))) real(reshape(psi,(:,1))) imag(reshape(psi,(:,1))) ]
            np.savetxt("psi$i.dat",ft,fmt="%+16.7e")
        end

    end


end


function Euler(vol,psi)

    phi = laplace_r(vol,dr,fkm,n0bc,nLbc,ar)
    cur = -laplace(psi,fkm,dr,ar)
    dy_psi = ydiff1(psi,fkm)

    vol1 = vol .+ dt*(fkp.*cur -ϵ*(ss./q+(2.0 .-s).*s ./(ar.*q)).*dy_psi .+ visvol.*laplace(psi,fkm,dr,ar) )
    psi1 = psi .+ dt*(-fkp.*phi .- vispsi .*cur)

    return vol1,psi1

end


function RK4(vol,psi)

    k1vol,k1psi = RHS(vol,             psi            )
    k2vol,k2psi = RHS(vol.+dt*k1vol/2, psi.+dt*k1psi/2)
    k3vol,k3psi = RHS(vol.+dt*k2vol/2, psi.+dt*k2psi/2)
    k4vol,k4psi = RHS(vol.+dt*k3vol,   psi.+dt*k3psi  )

    vol1 = vol .+ dt*(k1vol .+2k2vol .+2k3vol .+k4vol)/6
    psi1 = psi .+ dt*(k1psi .+2k2psi .+2k3psi .+k4psi)/6

    return vol1,psi1

end


function RHS(vol,psi)

    phi = laplace_r(vol,dr,fkm,n0bc,nLbc,ar)
    cur = -laplace(psi,fkm,dr,ar)
    dy_psi = ydiff1(psi,fkm)

    kvol = fkp.*cur -ϵ*(ss./q+(2.0 .-s).*s ./(ar.*q)).*dy_psi .+ visvol.*laplace(psi,fkm,dr,ar)
    kpsi = -fkp.*phi .- vispsi .*cur

    return kvol,kpsi

end

# timeloop()

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


