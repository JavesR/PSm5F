function MakeFk(q::Vector,lkm::Vector,lkn::Vector ,ar)
    
    
    fkm = zeros(ComplexF64,size(ar)[1],size(lkm)[1])
    fkn = zeros(ComplexF64,size(ar)[1],size(lkn)[1])
    fkp = zeros(ComplexF64,size(ar)[1],size(lkn)[1])


    for l in 1:size(fkm)[2]
        for i in 1:size(fkm)[1]
            fkm[i,l] = im*lkm[l]/ar[i]
            fkn[i,l] = im*lkn[l]/ar[i]
            fkp[i,l] = ϵ*(im*lkm[l]/q[i] - im*lkn[l])
        end
    end    

    return fkm,fkn,fkp

end


function mode_selection(qmax::Float64,qmin::Float64,
                        m_sh::Int64,n_sh::Int64,
                        m_max::Int64,n_max::Int64)
    
    flagGAM = 0
    lkm = Vector{Int64}(undef,m_max*n_max)
    lkn = Vector{Int64}(undef,m_max*n_max)

    lkm[1] = 0;lkn[1] = 0
    for l in 1:flagGAM
        lkm[l+1] = l
        lkn[l+1] = 0
    end
    
    l = flagGAM+1

    q_sh = m_sh/n_sh
    
    if q_sh > qmin && q_sh < qmax
        for m in 1:m_max, n in 1:n_max
            q = m/n
            if abs(q-q_sh) < 1.e-5 
                l = l+1
                lkm[l] = m
                lkn[l] = n
            end
        end
    end

    return l,lkm[1:l],lkn[1:l]

end

function mode_selection(qmax::Float64,qmin::Float64,
                        m_max::Int64,n_single::Int64)

    flagGAM = 0
    lkm = Vector{Int64}(undef,m_max)
    lkn = Vector{Int64}(undef,m_max)

    lkm[1] = 0;lkn[1] = 0

    for l in 1:flagGAM
        lkm[l+1] = l
        lkn[l+1] = 0
    end

    l = flagGAM+1

    for m in 1:m_max
        q = m / n_single

        if q > qmin && q < qmax
            l = l+1
            lkm[l] = m
            lkn[l] = n_single
        end
    
    end

    return l,lkm[1:l],lkn[1:l]

end