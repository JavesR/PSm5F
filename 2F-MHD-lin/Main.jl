using PyCall
using BenchmarkTools
using Printf

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


# ------q profile---------
# q = 1.05 .+ 2.0 .* r.^ 2
q0=1.5
qa=5.0
lambda=1.7
q=q0.*(1 .+r.^(2*lambda).*((qa/q0)^lambda-1)).^(1/lambda)
# ------------------------

s = ar ./ q .* rdiff1(q,dr,nr+1)
ss = rdiff1(s,dr,nr+1)

qmin = findmin(q)[1]
qmax = findmax(q)[1]
lmx,lkm,lkn = mode_selection(qmax,qmin,2,1,2,1)
fkm,fkn,fkp = MakeFk(q,lkm,lkn,ar)



visvol = 1.e-4
vispsi = 1.e-4

dt = 0.001
tmax = 200

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

    timeloop(vol,phi,psi,cur,nt)

end


function timeloop(vol,phi,psi,cur,nt)

    t = 0.0

    contour_data("t0.dat",vol,phi,psi,cur)

    for i in 1:nt

        # vol,phi,psi,cur = Euler(vol,phi,psi,cur)
        vol,phi,psi,cur = MEuler(vol,phi,psi,cur)
        # vol,phi,psi,cur= RK4(vol,phi,psi,cur)

        t = t+dt

        if i%1000 == 0
            println("t = $t")
            filename = @sprintf "t%g.dat" t
            contour_data(filename,vol,phi,psi,cur)
        end

    end


end



function contour_data(filename::String,
                        vol::Matrix,
                        phi::Matrix,
                        psi::Matrix,
                        cur::Matrix)

    fp = open(filename,"w")
        for l in 1:lmx
            for i in 1:nr+1
                @printf(fp,"\t %d", lkm[l])
                @printf(fp,"\t %d", lkn[l])
                @printf(fp,"\t %1.6e",real(phi[i,l]))
                @printf(fp,"\t %1.6e",imag(phi[i,l]))
                @printf(fp,"\t %1.6e",real(psi[i,l]))
                @printf(fp,"\t %1.6e",imag(psi[i,l]))
                @printf(fp,"\t %1.6e",real(vol[i,l]))
                @printf(fp,"\t %1.6e",imag(vol[i,l]))
                @printf(fp,"\t %1.6e",real(cur[i,l]))
                @printf(fp,"\t %1.6e",imag(cur[i,l]))
                @printf(fp,"\n")
            end
        end
    close(fp)

end



flow(tmax,dt)



