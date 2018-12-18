function rdiff1(f::Vector,dr::Float64,n::Int64)

    # return center2(f,dr,n)
    return center4(f,dr,n)

end


function rdiff1(f::Matrix,dr::Float64)

    return center4(f,dr)

end


function rdiff2(f::Matrix,dr::Float64)

    return center4_2(f,dr)

end


function center2(f::Vector,dr::Float64,n::Int64)

    df = zeros(eltype(f),n)

    for i in 2:n-1
        df[i] =  (f[i+1]-f[i-1])/(2dr)
    end

    # 1st bc
    df[1] = (-3f[1].+4f[2].-f[3])/(2dr)
    df[end] = (f[end-2].-4f[end-1].+3f[end])/(2dr)

    # 2nd bc
    # df[1] = (f[2] - f[end-1])/(2dr)
    # df[end] = (f[2] - f[end-1])/(2dr)

    return df

end


function center4(f::Vector,dr::Float64,n::Int64)

    df = zeros(eltype(f),n)

    for i in 3:n-2
        df[i] = ( -f[i+2] + 8f[i+1] - 8f[i-1] + f[i-2] )/(12dr)
    end

    # 1st bc
    df[1] = ( -25f[1] + 48f[2] - 36f[3] + 16f[4] -3f[5] )/(12dr)
    df[2] = ( -25f[2] + 48f[3] - 36f[4] + 16f[5] -3f[6] )/(12dr)
    df[n] = ( -25f[n] + 48f[n-1] - 36f[n-2] + 16f[n-3] -3f[n-4] )/(-12dr)
    df[n-1] = ( -25f[n-1] + 48f[n-2] - 36f[n-3] + 16f[n-4] -3f[n-5] )/(-12dr)


    # 2nd bc
    # df[1] = ( -f[3] + 8f[2] - 8f[n-1] + f[n-2] )/(12dr)
    # df[2] = ( -f[4] + 8f[3] - 8f[n] + f[n-1] )/(12dr)

    # df[n] = ( -f[3] + 8f[2] - 8f[n-1] + f[n-2] )/(12dr)
    # df[n-1] = ( -f[2] + 8f[n] - 8f[n-2] + f[n-3])/(12dr)

    return df

end


function center4(f::Matrix,dr::Float64)

    df = zeros(eltype(f),size(f))

    for l in 1:size(f,2)
        for i in 3:size(f,1)-2

            df[i,l] = ( -f[i+2,l] + 8f[i+1,l] - 8f[i-1,l] + f[i-2,l] )/(12dr)
        
        end

        df[1,l] = ( -25f[1,l] + 48f[2,l] - 36f[3,l] + 16f[4,l] -3f[5,l] )/(12dr)
        df[2,l] = ( -25f[2,l] + 48f[3,l] - 36f[4,l] + 16f[5,l] -3f[6,l] )/(12dr)
        df[end,l] = ( -25f[end,l] + 48f[end-1,l] - 36f[end-2,l] + 16f[end-3,l] -3f[end-4,l] )/(-12dr)
        df[end-1,l] = ( -25f[end-1,l] + 48f[end-2,l] - 36f[end-3,l] + 16f[end-4,l] -3f[end-5,l] )/(-12dr)
    
    end

    return df

end


function center4_2(f::Matrix,dr::Float64)

    d2f = zeros(eltype(f),size(f))

    for l in 1:size(f,2)
        for i in 3:size(f,1)-2
            d2f[i,l] = ( -f[i+2,l] + 16f[i+1,l] -30f[i,l] + 16f[i-1,l] - f[i-2,l] )/(12dr^2)
        end

        d2f[1,l] = ( 2f[1,l] - 5f[2,l] + 4f[3,l] -f[4,l])/(dr^2)
        d2f[2,l] = ( f[1,l] - 2f[2,l] + f[3,l] )/(dr^2)

        d2f[end,l] = ( 2f[end,l] - 5f[end-1,l] + 4f[end-2,l] - f[end-3,l] )/(dr^2)
        d2f[end-1,l] = ( f[end,l] - 2f[end-1,l] + f[end-2,l] )/(dr^2) 
    end

    return d2f

end


function ydiff1(f::Matrix,fkm::Matrix)

    return f.*fkm
    
end


function laplace(f::Matrix,fkm::Matrix,dr::Float64,
    ar::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}})

    d2f =  rdiff2(f,dr) .+  rdiff1(f,dr) ./ ar 

    return d2f .+ f.*fkm.^2

end


function laplace_r(f::Matrix,dr::Float64,fkm::Matrix,n0bc::Vector,nLbc::Vector,
    ar::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}})

    v = zeros(eltype(f),size(f))

    a =  1/dr^2 .-1/(2dr)./ar
    b =  -2/dr^2 .+fkm.^2
    c =  1/dr^2 .+1/(2dr)./ar


    v[1,:] = n0bc
    v[end,:] = nLbc

    v = tridag(a,b,c,f,nr+1,lmx,n0bc,nLbc)

    return v

end


function tridag(a::Vector,b::Matrix,c::Vector,d::Matrix,n::Int64,lmx::Int64,
        n0bc,nLbc)

    u = zeros(ComplexF64,n,lmx)
    e = zeros(ComplexF64,n,lmx)
    f = zeros(ComplexF64,n,lmx)
    
    nr = n-1

    e[nr,:] .= 0.0
    f[nr,:] = nLbc

    for i in nr-1:-1:1
        e[i,:] = -a[i+1] ./(b[i+1,:] .+c[i+1]*e[i+1,:])
        f[i,:] = (d[i+1,:] .-c[i+1].*f[i+1,:])./(b[i+1,:] .+c[i+1].*e[i+1,:])
    end

    for i in 2:n
        u[i,:] = e[i-1,:].*u[i-1,:] .+ f[i-1,:]
    end
    
    return u

end






function pus(u::Matrix,rhs::Matrix,dt,vis::Matrix,dr,ar,
    fkm::Matrix,n0bc::Vector,nLbc::Vector,beta::Float64)

    v = zeros(eltype(u),size(u))

    a =  -dt*vis.*(1/dr^2 .-1/(2dr)./ar)
    b =beta .-dt*vis.*(-2/dr^2 .+fkm.^2)
    c =  -dt*vis.*(1/dr^2 .+1/(2dr)./ar)
    d = dt*rhs .+ beta*u


    v[1,:] = n0bc
    v[end,:] = nLbc

    v = tridag(a,b,c,d,nr+1,lmx,n0bc,nLbc)

    return v

end

function tridag(a::Matrix,b::Matrix,c::Matrix,d::Matrix,n::Int64,lmx::Int64,
    n0bc,nLbc)

    u = zeros(ComplexF64,n,lmx)
    e = zeros(ComplexF64,n,lmx)
    f = zeros(ComplexF64,n,lmx)

    nr = n-1

    e[nr,:] .= 0.0
    f[nr,:] = nLbc

    for i in nr-1:-1:1
        e[i,:] = -a[i+1,:] ./(b[i+1,:] .+c[i+1,:].*e[i+1,:])
        f[i,:] = (d[i+1,:] .-c[i+1,:].*f[i+1,:])./(b[i+1,:] .+c[i+1,:].*e[i+1,:])
    end

    for i in 2:n
        u[i,:] = e[i-1,:].*u[i-1,:] .+ f[i-1,:]
    end

    return u

end