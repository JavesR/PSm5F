
function Euler( vol::Matrix,
                phi::Matrix,
                psi::Matrix,
                cur::Matrix)

    dy_psi = ydiff1(psi,fkm)

    vol = vol .+ dt*(fkp.*cur -ϵ*(ss./q+(2.0 .-s).*s ./(ar.*q)).*dy_psi .+ visvol.*laplace(psi,fkm,dr,ar) )
    psi = psi .+ dt*(-fkp.*phi .- vispsi .*cur) ./β
    phi = laplace_r(vol,dr,fkm,n0bc,nLbc,ar)
    cur = -laplace(psi,fkm,dr,ar)

    return vol,phi,psi,cur

end


function MEuler(vol::Matrix,
                phi::Matrix,
                psi::Matrix,
                cur::Matrix)

    vol_old = vol[:,:]
    psi_old = psi[:,:]


    dy_psi = ydiff1(psi,fkm)

    volrhs = fkp.*cur  -ϵ*(ss./q+(2.0 .-s).*s ./(ar.*q)).*dy_psi
    psirhs = -fkp.*phi

    vol = pus(vol_old,volrhs,0.5*dt,visvol,dr,ar,fkm,n0bc,nLbc,1.0)
    psi = pus(psi_old,psirhs,0.5*dt,vispsi,dr,ar,fkm,n0bc,nLbc,β)
    phi = laplace_r(vol,dr,fkm,n0bc,nLbc,ar)
    cur = -laplace(psi,fkm,dr,ar)


    dy_psi = ydiff1(psi,fkm)

    volrhs = fkp.*cur  -ϵ*(ss./q+(2.0 .-s).*s ./(ar.*q)).*dy_psi
    psirhs = -fkp.*phi

    vol = pus(vol_old,volrhs,dt,visvol,dr,ar,fkm,n0bc,nLbc,1.0)
    psi = pus(psi_old,psirhs,dt,vispsi,dr,ar,fkm,n0bc,nLbc,β)
    phi = laplace_r(vol,dr,fkm,n0bc,nLbc,ar)
    cur = -laplace(psi,fkm,dr,ar)


    return vol,phi,psi,cur

end



function RK4(vol::Matrix,
            phi::Matrix,
            psi::Matrix,
            cur::Matrix)

    k1vol,k1psi,phi,cur = RHS(vol,             psi            ,phi,cur)
    k2vol,k2psi,phi,cur = RHS(vol.+dt*k1vol/2, psi.+dt*k1psi/2,phi,cur)
    k3vol,k3psi,phi,cur = RHS(vol.+dt*k2vol/2, psi.+dt*k2psi/2,phi,cur)
    k4vol,k4psi,phi,cur = RHS(vol.+dt*k3vol,   psi.+dt*k3psi  ,phi,cur)

    vol = vol .+ dt*(k1vol .+2k2vol .+2k3vol .+k4vol)/6
    psi = psi .+ dt*(k1psi .+2k2psi .+2k3psi .+k4psi)/6
    phi = laplace_r(vol,dr,fkm,n0bc,nLbc,ar)
    cur = -laplace(psi,fkm,dr,ar)

    return vol,phi,psi,cur

end


function RHS(vol::Matrix,
            psi::Matrix,
            phi::Matrix,
            cur::Matrix)

    dy_psi = ydiff1(psi,fkm)

    kvol = fkp.*cur -ϵ*(ss./q+(2.0 .-s).*s ./(ar.*q)).*dy_psi .+ visvol.*laplace(psi,fkm,dr,ar)
    kpsi = (-fkp.*phi .- vispsi .*cur) ./β
    phi = laplace_r(vol,dr,fkm,n0bc,nLbc,ar)
    cur = -laplace(psi,fkm,dr,ar)

    return kvol,kpsi,phi,cur

end