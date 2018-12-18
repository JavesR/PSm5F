
function Euler(vol::Matrix,psi::Matrix)

    phi = laplace_r(vol,dr,fkm,n0bc,nLbc,ar)
    cur = -laplace(psi,fkm,dr,ar)
    dy_psi = ydiff1(psi,fkm)

    vol1 = vol .+ dt*(fkp.*cur -ϵ*(ss./q+(2.0 .-s).*s ./(ar.*q)).*dy_psi .+ visvol.*laplace(psi,fkm,dr,ar) )
    psi1 = psi .+ dt*(-fkp.*phi .- vispsi .*cur) ./β

    return vol1,psi1

end


function MEuler(vol::Matrix,psi::Matrix)

    vol_old = vol[:,:]
    psi_old = psi[:,:]

    phi = laplace_r(vol,dr,fkm,n0bc,nLbc,ar)
    cur = -laplace(psi,fkm,dr,ar)
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

    return vol,psi

end



function RK4(vol::Matrix,psi::Matrix)

    k1vol,k1psi = RHS(vol,             psi            )
    k2vol,k2psi = RHS(vol.+dt*k1vol/2, psi.+dt*k1psi/2)
    k3vol,k3psi = RHS(vol.+dt*k2vol/2, psi.+dt*k2psi/2)
    k4vol,k4psi = RHS(vol.+dt*k3vol,   psi.+dt*k3psi  )

    vol = vol .+ dt*(k1vol .+2k2vol .+2k3vol .+k4vol)/6
    psi = psi .+ dt*(k1psi .+2k2psi .+2k3psi .+k4psi)/6

    return vol,psi

end


function RHS(vol::Matrix,psi::Matrix)

    phi = laplace_r(vol,dr,fkm,n0bc,nLbc,ar)
    cur = -laplace(psi,fkm,dr,ar)
    dy_psi = ydiff1(psi,fkm)

    kvol = fkp.*cur -ϵ*(ss./q+(2.0 .-s).*s ./(ar.*q)).*dy_psi .+ visvol.*laplace(psi,fkm,dr,ar)
    kpsi = (-fkp.*phi .- vispsi .*cur) ./β

    return kvol,kpsi

end