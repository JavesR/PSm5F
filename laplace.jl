
function center2(f,dr,n)

    df = zeros(eltype(f),n)

    for i in 2:n-1
        df[i] =  (f[i+1] -f[i-1]) /(2dr)
    end

    df[1] = (-3f[1]+4f[2]-f[3])/(2dr)
    df[n] = (f[n-2]-4f[n-1]+3f[n])/(2dr)

    return df

end


function poisson(f,dr,n0bc,nLbc)

    df = zeros(eltype(f),size(f))
    a = zeros(eltype(f),size(f))
    b = zeros(eltype(f),size(f))
    c = zeros(eltype(f),size(f))

    a .=  1/dr^2 
    b .=  -2/dr^2
    c .=  1/dr^2 

    df[1] = n0bc
    df[end] = nLbc

    tridag(a,b,c,f,df,n0bc,nLbc)

    return df
    
end


function poisson2(f,dr,n0bc,nLbc)

    df = zeros(eltype(f),size(f))
    a = zeros(eltype(f),size(f))
    b = zeros(eltype(f),size(f))
    c = zeros(eltype(f),size(f))

    f_new = zeros(eltype(f),size(f))

    a .=  1/dr^2 
    b .=  -2/dr^2
    c .=  1/dr^2 

    df[1] = n0bc
    df[end] = nLbc

    for i in 2:size(f,1)-1
        f_new[i] = f[i] + dr^2/12*(f[i+1] - 2f[i] + f[i-1])/(dr^2)
    end

    tridag(a,b,c,f_new,df,n0bc,nLbc)

    return df
    
end



function tridag(a,b,c,d,df
        ,n0bc,nLbc)

    e = zeros(ComplexF64,size(d))
    f = zeros(ComplexF64,size(d))
    
    nr = size(d,1)-1

    e[nr] = 0.0
    f[nr] = nLbc

    for i in nr-1:-1:1
        e[i] = -a[i+1] /(b[i+1] +c[i+1]*e[i+1])
        f[i] = (d[i+1] -c[i+1]*f[i+1])/(b[i+1] +c[i+1]*e[i+1])
    end

    for i in 2:nr+1
        df[i] = e[i-1]*df[i-1] + f[i-1]
    end
    
    return nothing

end


using PyCall
@pyimport matplotlib.pyplot as plt

dr = 0.01
r = 0.0:dr:1.0

vol = sin.(π*r)
phi0 = -sin.(π*r)/π^2
phi1 = poisson(vol,dr,0.0,0.0)
phi2 = poisson2(vol,dr,0.0,0.0)


# plt.plot(r,vol)
plt.plot(r,abs.(phi1-phi0),"r")
plt.plot(r,abs.(phi2-phi0),"b--")
# plt.plot(r,phi0)

plt.yscale("log")
plt.grid()
plt.show()