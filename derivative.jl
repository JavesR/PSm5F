function tridag1(a,b,c,f,n)

    α = zeros(Float64,n)
    β = zeros(Float64,n)
    γ = zeros(Float64,n)
    x = zeros(Float64,n)
    y = zeros(Float64,n)
    
    α[1] = b[1]
    β[1] = c[1]/α[1]

    for i in 2:n-1
        α[i] = b[i] - a[i]*β[i-1]
        β[i] = c[i]/α[i]
        γ[i] = a[i]
    end

    α[n] = b[n] - a[n]*β[n-1]
    γ[n] = a[n]


    y[1] = f[1]/α[1]
    for i in 2:n
        y[i] = ( f[i] - γ[i]*y[i-1] ) / α[i]
    end

    x[n] = y[n]
    for i in n-1:-1:1
        x[i] = y[i] - β[i]*x[i+1]
    end

    return x
    
end


function tridag2(a,b,c,d,n,n0bc,nLbc)

    e = zeros(Float64,n)
    f = zeros(Float64,n)
    x = zeros(Float64,n)

    nm = n-1
    e[nm] = 0.0
    f[nm] = nLbc

    for i in nm-1:-1:1
        e[i] = -a[i+1]/(b[i+1]+c[i+1]*e[i+1])
        f[i] = (d[i+1]-c[i+1]*f[i+1])/(b[i+1]+c[i+1]*e[i+1])
    end

    for i in 2:n
        x[i] = e[i-1]*x[i-1] + f[i-1]
    end

    x[1] = n0bc

    return x

end

function pade_1(f,dr,n)

    a = dr/3*ones(Float64,n)
    a[1] = 0.
    a[n] = 0.
    
    b = 4dr/3*ones(Float64,n)
    b[1] = 1.
    b[n] = 1.
    
    c = dr/3*ones(Float64,n)
    c[1] = 0.
    c[n] = 0.
    
    u = zeros(eltype(f),n)
    for i in 2:n-1
        u[i] = f[i+1]-f[i-1]
    end

    u[1] = ( -25f[1] + 48f[2] - 36f[3] + 16f[4] -3f[5] )/(12dr)
    u[n] = ( -25f[n] + 48f[n-1] - 36f[n-2] + 16f[n-3] -3f[n-4] )/(-12dr)

    df = tridag1(a,b,c,u,n)

    return df

end


function pade_2(df,dr,n)

    a = dr/3*ones(Float64,n)
    a[1] = 0.
    a[n] = 0.
    
    b = 4dr/3*ones(Float64,n)
    b[1] = 1.
    b[n] = 1.
    
    c = dr/3*ones(Float64,n)
    c[1] = 0.
    c[n] = 0.
    
    u = zeros(eltype(df),n)
    for i in 2:n-1
        u[i] = df[i+1]-df[i-1]
    end
    u[1] = ( -25/12*df[1]+4df[2]-3df[3]+4/3*df[4]-1/4*df[5] )/dr
    u[n] = ( 25/12*df[n]-4df[n-1]+3df[n-2]-4/3*df[n-3]+1/4*df[n-4] )/dr

    d2f = tridag1(a,b,c,u,n)

    return d2f

end



function center2(f,dr,n)

    df = zeros(eltype(f),n)

    for i in 2:n-1
        df[i] =  (f[i+1] -f[i-1]) /(2dr)
    end

    df[1] = (-3f[1]+4f[2]-f[3])/(2dr)
    df[n] = (f[n-2]-4f[n-1]+3f[n])/(2dr)

    return df

end



function center4(f,dr,n)

    df = zeros(eltype(f),n)

    for i in 3:n-2
        df[i] = ( -f[i+2] + 8f[i+1] - 8f[i-1] + f[i-2] )/(12dr)
    end

    # df[1] = ( -3f[1] + 4f[2] - f[3] )/(2dr)
    # df[2] = ( -2f[1] - 3f[2] + 6f[3] - f[4] )/(6dr)
    # df[n] = ( -3f[n] + 4f[n-1] - f[n-2] )/(-2dr)
    # df[n-1] = ( -2f[n] - 3f[n-1] + 6f[n-2] - f[n-3] )/(-6dr) 

    df[1] = ( -25f[1] + 48f[2] - 36f[3] + 16f[4] -3f[5] )/(12dr)
    df[2] = ( -25f[2] + 48f[3] - 36f[4] + 16f[5] -3f[6] )/(12dr)
    df[n] = ( -25f[n] + 48f[n-1] - 36f[n-2] + 16f[n-3] -3f[n-4] )/(-12dr)
    df[n-1] = ( -25f[n-1] + 48f[n-2] - 36f[n-3] + 16f[n-4] -3f[n-5] )/(-12dr)

    
    return df

end

function center4_2(f,dr,n)

    d2f = zeros(eltype(f),n)

    for i in 3:n-2
        d2f[i] = ( -f[i+2] + 16f[i+1] -30f[i] + 16f[i-1] - f[i-2] )/(12dr^2)
    end

    d2f[1] = ( 2f[1] - 5f[2] + 4f[3] -f[4])/(dr^2)
    d2f[2] = ( f[1] - 2f[2] + f[3] )/(dr^2)

    d2f[n] = ( 2f[n] - 5f[n-1] + 4f[n-2] - f[n-3] )/(dr^2)
    d2f[n-1] = ( f[n] - 2f[n-1] + f[n-2] )/(dr^2) 

    return d2f

end




using PyCall

@pyimport matplotlib.pyplot as plt

nr = 10

dr = 1/nr
r = 0:dr:1
phi = sin.(π*r) 
dphi0 = π*cos.(π*r)
d2phi0 = - π^2*sin.(π*r)

# plt.plot(r,dphi0)

dphi1 = pade_1(phi,dr,nr+1)
dphi2 = center2(phi,dr,nr+1)
dphi3 = center4(phi,dr,nr+1)

plt.plot(r,dphi1-dphi0,"o-",ms=4, lw=1, alpha=0.6,label="pade")
plt.plot(r,dphi2-dphi0,"o-",ms=4, lw=1, alpha=0.6,label="cer2")
plt.plot(r,dphi3-dphi0,"o-",ms=4, lw=1, alpha=0.6,label="cer4")



# d2phi1 = pade_2(dphi1,dr,nr+1)
# d2phi2 = center4_2(phi,dr,nr+1)

# plt.plot(r,d2phi1-d2phi0)
# plt.plot(r,d2phi2-d2phi0)

plt.legend()
plt.grid()

plt.show()