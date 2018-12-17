using PyCall

@pyimport matplotlib.pyplot as plt

x0 = 0.0 # left
xL = 1.0 # right

nx = 100
dx = (xL - x0) / nx
x = x0:dx:xL

Tmax = 1.0
dt = 0.01
nt = Int64(Tmax/dt)


# convective coefficient
α = 1.0

#CFL number
println("CFL number = $(α*dt/dx)")

# initial condition π
u = sin.(2*π .* x)


function center2(f,dr,n)

    df = zeros(eltype(f),n)

    for i in 2:n-1
        df[i] =  (f[i+1]-f[i-1])/(2dr)
    end

    # 1st bc
    # df[1] = (-3f[1].+4f[2].-f[3])/(2dr)
    # df[end] = (f[end-2].-4f[end-1].+3f[end])/(2dr)

    # 2nd bc
    df[1] = (f[2] - f[end-1])/(2dr)
    df[end] = (f[2] - f[end-1])/(2dr)

    return df

end


function center4(f,dr,n)

    df = zeros(eltype(f),n)

    for i in 3:n-2
        df[i] = ( -f[i+2] + 8f[i+1] - 8f[i-1] + f[i-2] )/(12dr)
    end


    # 1st bc
    # df[1] = ( -3f[1] + 4f[2] - f[3] )/(2dr)
    # df[2] = ( -2f[1] - 3f[2] + 6f[3] - f[4] )/(6dr)

    # df[n] = ( -3f[n] + 4f[n-1] - f[n-2] )/(-2dr)
    # df[n-1] = ( -2f[n] - 3f[n-1] + 6f[n-2] - f[n-3] )/(-6dr) 

    # 2nd bc
    df[1] = ( -f[3] + 8f[2] - 8f[n-1] + f[n-2] )/(12dr)
    df[2] = ( -f[4] + 8f[3] - 8f[n] + f[n-1] )/(12dr)

    df[n] = ( -f[3] + 8f[2] - 8f[n-1] + f[n-2] )/(12dr)
    df[n-1] = ( -f[2] + 8f[n] - 8f[n-2] + f[n-3])/(12dr)

    # plt.plot(df)

    # plt.show()

    return df

end



function Euler(u,n,dx,dt,α)

    du = center2(u,dx,n)

    u1 = u .- α*dt*du
    
    return u1

end


function MEuler(u,n,dx,dt,α)

    du = center2(u,dx,n)
    u1 = u .- α*0.5*dt*du

    du1 = center2(u1,dx,n)
    u2 = u .- α*dt*du1

    return u2

end


function RK4(u,n,dx,dt,α)

    k1 = RHS(u, dx, n,α)
    k2 = RHS(u.+dt*k1/2, dx, n,α)
    k3 = RHS(u.+dt*k2/2, dx, n,α)
    k4 = RHS(u.+dt*k3, dx, n,α)

    u1 = u .+ dt*(k1 .+2k2 .+2k3 .+k4)/6

    return u1

end


function RHS(u,dx, n,α)

    # du = center2(u,dx,n)
    du = center4(u,dx,n)
    k = -α.*du

    return k

end


function timeloop(u,n,dx,dt,α,nt)

    t = 0.0

    plt.plot(u,lw=1.5,c="r")

    for i in 1:nt

        # u = Euler(u,n,dx,dt,α)
        # u = MEuler(u,n,dx,dt,α)
        u = RK4(u,n,dx,dt,α)
 
        t = t+dt

        if i%(Int64(nt/10)) == 0
            plt.plot(u,lw=1.0,"--",label="i=$i")
        end

    end


    return

end



timeloop(u,nx+1,dx,dt,α,nt)

plt.grid()
plt.legend()
plt.show()