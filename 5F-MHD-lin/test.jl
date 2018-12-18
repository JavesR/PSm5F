# using Plots

# Domain
x0 = 0.0 # left
xL = 1.0 # right

# number of points
nx = 20

# grid spacing (spatial)
dx = (xL - x0) / nx

# spatial coordinates 
x = x0:dx:xL

# maximum time desired
Tmax = 2.0
# time step
dt = 0.001
# number of points in time
nt = Int64(Tmax/dt)


# diffusion coefficient
α = 1.0
# u: temperature variable 
u = Vector(undef,nx+1)
# initial condition π
u = sin.(π .* x)

# boundary conditions: holds for all time
u[1] = 0.0
u[end]= 0.0

r = α*dt/(dx*dx)

# source term (nonhomogenous forcing term in heat equation)
f(x,t) =  (π^2-1.0)*exp(-t)*sin(π*x) 


function timeloop(nt,nx,u)

    t = 0.0

    fp = open("test.txt","a")
    for j in 1:nx+1
        write(fp,string(u[j])*" ")
    end    
    write(fp,"\n") 

    for i in 1:nt
        if i % 500 == 0     
            for j in 1:nx+1
                write(fp,string(u[j])*" ")
            end    
            write(fp,"\n")
        end
        
        u = FTCS(nx,u,t)
        t = t+dt

    end

    close(fp)

    return 
end




function FTCS(nx,u,t)

    global r

    v = zeros(Float64,nx+1)

    # update (t=n+1)
    for i in 2:nx
        v[i] = u[i] + r*(u[i+1]-2.0*u[i]+u[i-1]) + dt*f(x[i],t)
    end
    
    return v

end



timeloop(nt,nx,u)