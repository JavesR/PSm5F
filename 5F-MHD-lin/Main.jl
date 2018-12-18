using PyCall
using BenchmarkTools

@pyimport numpy as np
@pyimport matplotlib.pyplot as plt

include("./FunEqm.jl")
include("./FunDev.jl")
include("./FunScm.jl")

ϵ = 0.25
β = 0.01
τ = 1.0
Γ = 5/3
a_minor = 80.


n0 = 0.05
nL = 1.0
nr = 256
dn = (nL-n0)/nr
r = n0:dn:nL

dr = a_minor*dn
ar = a_minor*r


q = 1.05 .+ 2.0 .* r.^ 2
den0 = 0.8 .+ 0.2.*exp.(-2 .*r.^2)
tem0 = 0.35 .+ 0.65.*(1 .-r.^2).^2

s = ar ./ q .* rdiff1(q,dr,nr+1)
ss = rdiff1(s,dr,nr+1)

ηi = den0 ./tem0 .* rdiff1(tem0,dr,nr+1) ./ rdiff1(den0,dr,nr+1)

aln = a_minor .*rdiff1(den0,dr,nr+1) ./ den0
tn0 = tem0 ./ den0


q_sn = Int64(round(0.2*(nr+1)))
q_en = Int64(round(0.8*(nr+1)))

qmin = findmin( q[q_sn:q_en] )[1]
qmax = findmax( q[q_sn:q_en] )[1]
# lmx,lkm,lkn = mode_selection(qmax,qmin,2,1,2,1)  # model 1
lmx,lkm,lkn = mode_selection(qmax,qmin,50,17)      # model 2
fkm,fkn,fkp = MakeFk(q,lkm,lkn,ar)

println(lkm)
println(lkn)

visden = 1.e-7 .* reshape(lkm.^4,(1,:)) .* ones(Float64,nr+1,lmx)
visvol = 1.e-7 .* reshape(lkm.^4,(1,:)) .* ones(Float64,nr+1,lmx)
visval = 1.e-7 .* reshape(lkm.^4,(1,:)) .* ones(Float64,nr+1,lmx)
vistem = 1.e-7 .* reshape(lkm.^4,(1,:)) .* ones(Float64,nr+1,lmx)
vispsi = 4.e-5 .* ones(Float64,nr+1,lmx)


dt = 0.001
tmax = 300

n0bc = zeros(ComplexF64,lmx)
nLbc = zeros(ComplexF64,lmx)


function flow(tmax,dt)

    # initialization
    den = zeros(ComplexF64,(nr+1,lmx))
    phi = zeros(ComplexF64,(nr+1,lmx))
    val = zeros(ComplexF64,(nr+1,lmx))
    psi = zeros(ComplexF64,(nr+1,lmx))
    tem = zeros(ComplexF64,(nr+1,lmx))

    den[:,2:end] .= 1.e-5 * sin.(π*(0:nr)/nr) .+ 0im
    phi[:,2:end] .= 1.e-5 * sin.(π*(0:nr)/nr) .+ 0im
    val[:,2:end] .= 1.e-5 * sin.(π*(0:nr)/nr) .+ 0im
    psi[:,2:end] .= 1.e-5 * sin.(π*(0:nr)/nr) .+ 0im
    tem[:,2:end] .= 1.e-5 * sin.(π*(0:nr)/nr) .+ 0im



    vol =  laplace(phi,fkm,dr,ar)

    nt = Int64(tmax/dt)

    timeloop(den,vol,val,psi,tem,nt)

end


function timeloop(den,vol,val,psi,tem,nt)

    t = 0.0

    print("t=$t \n")
    ft = [  real(reshape(vol,(:,1))) imag(reshape(vol,(:,1))) real(reshape(psi,(:,1))) imag(reshape(psi,(:,1))) ]
    np.savetxt("t$t.dat",ft,fmt="%+16.7e")

    for i in 1:nt

        # vol,psi = Euler(vol,psi)
        den,vol,val,psi,tem = MEuler(den,vol,val,psi,tem)
        # den,vol,val,psi,tem= RK4(den,vol,val,psi,tem)

        t = t+dt

        if i%1000 == 0
			
			phi = laplace_r(vol,dr,fkm,n0bc,nLbc,ar)
			cur = -laplace(psi,fkm,dr,ar)

            print("t=$t \n")
            ft = [  real(reshape(cur,(:,1))) imag(reshape(cur,(:,1))) real(reshape(psi,(:,1))) imag(reshape(psi,(:,1))) ]
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


