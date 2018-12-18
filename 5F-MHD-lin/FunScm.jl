
function Euler(vol::Matrix,psi::Matrix)

    phi = laplace_r(vol,dr,fkm,n0bc,nLbc,ar)
    cur = -laplace(psi,fkm,dr,ar)
    dy_psi = ydiff1(psi,fkm)

    vol1 = vol .+ dt*(fkp.*cur -ϵ*(ss./q+(2.0 .-s).*s ./(ar.*q)).*dy_psi .+ visvol.*laplace(psi,fkm,dr,ar) )
    psi1 = psi .+ dt*(-fkp.*phi .- vispsi .*cur) ./β

    return vol1,psi1

end


function MEuler(den::Matrix,
                vol::Matrix,
                val::Matrix,
                psi::Matrix,
                tem::Matrix)

    den_old = den[:,:]
    vol_old = vol[:,:]
    val_old = val[:,:]
    psi_old = psi[:,:]
    tem_old = tem[:,:]


    phi = laplace_r(vol,dr,fkm,n0bc,nLbc,ar)
    cur = -laplace(psi,fkm,dr,ar)
    
    denrhs =(den0.*aln.*fkm.*phi
            -den0.*fkp.*val
            +fkp.*cur
            -ϵ*(ss./q+(2.0 .-s).*s ./(ar.*q)).*fkm.*psi )

    volrhs =(-tem0.*aln.*(1.0.+ηi).*fkm.*vol
            +(1 ./den0).*(fkp.*cur-ϵ*(ss./q+(2.0 .-s).*s ./(ar.*q)).*fkm.*psi) )
    
    valrhs =(-fkp.*(tem.+(1+τ).*tn0.*den)
            -β.*tem0.*aln.*(1 .+τ .+ηi).*fkm.*psi ) 

    psirhs =(-fkp.*(phi.-τ.*tn0.*den)
            +β.*tem0.*τ.*aln.*fkm.*psi
            +0.029.*τ^0.5.*abs.(fkp).*(val-cur./den0) )

    temrhs =(tem0.*aln.*ηi.*fkm.*phi
            -(Γ-1).*tem0.*fkp.*val
            -(Γ-1).*(8 .*tem0/π).^0.5.*abs.(fkp).*tem)

    den = pus(den_old,denrhs,0.5*dt,visden,dr,ar,fkm,n0bc,nLbc,1.0)
    vol = pus(vol_old,volrhs,0.5*dt,visvol,dr,ar,fkm,n0bc,nLbc,1.0)
    val = pus(val_old,valrhs,0.5*dt,visval,dr,ar,fkm,n0bc,nLbc,1.0)
    psi = pus(psi_old,psirhs,0.5*dt,vispsi,dr,ar,fkm,n0bc,nLbc,β)
    tem = pus(tem_old,temrhs,0.5*dt,vistem,dr,ar,fkm,n0bc,nLbc,1.0)



    phi = laplace_r(vol,dr,fkm,n0bc,nLbc,ar)
    cur = -laplace(psi,fkm,dr,ar)

    denrhs =(den0.*aln.*fkm.*phi
            -den0.*fkp.*val
            +fkp.*cur
            -ϵ*(ss./q+(2.0 .-s).*s ./(ar.*q)).*fkm.*psi )

    volrhs =(-tem0.*aln.*(1.0.+ηi).*fkm.*vol
            +(1 ./den0).*(fkp.*cur-ϵ*(ss./q+(2.0 .-s).*s ./(ar.*q)).*fkm.*psi) )
    
    valrhs =(-fkp.*(tem.+(1+τ).*tn0.*den)
            -β.*tem0.*aln.*(1 .+τ .+ηi).*fkm.*psi ) 

    psirhs =(-fkp.*(phi.-τ.*tn0.*den)
            +β.*tem0.*τ.*aln.*fkm.*psi
            +0.029.*τ^0.5.*abs.(fkp).*(val-cur./den0) )

    temrhs =(tem0.*aln.*ηi.*fkm.*phi
            -(Γ-1).*tem0.*fkp.*val
            -(Γ-1).*(8 .*tem0/π).^0.5.*abs.(fkp).*tem)

    den = pus(den_old,denrhs,dt,visden,dr,ar,fkm,n0bc,nLbc,1.0)
    vol = pus(vol_old,volrhs,dt,visvol,dr,ar,fkm,n0bc,nLbc,1.0)
    val = pus(val_old,valrhs,dt,visval,dr,ar,fkm,n0bc,nLbc,1.0)
    psi = pus(psi_old,psirhs,dt,vispsi,dr,ar,fkm,n0bc,nLbc,β)
    tem = pus(tem_old,temrhs,dt,vistem,dr,ar,fkm,n0bc,nLbc,1.0)

    return den,vol,val,psi,tem

end



# function RK4(den::Matrix,
#             vol::Matrix,
#             val::Matrix,
#             psi::Matrix,
#             tem::Matrix)

#     k1den,k1vol,k1val,k1psi,k1tem = RHS(den,
#                                         vol,
#                                         val,
#                                         psi,
#                                         tem )

#     k2den,k2vol,k2val,k2psi,k2tem = RHS(den.+dt*k1den/2, 
#                                         vol.+dt*k1vol/2,
#                                         val.+dt*k1val/2,
#                                         psi.+dt*k1psi/2,
#                                         tem.+dt*k1tem/2)

#     k3den,k3vol,k3val,k3psi,k3tem = RHS(den.+dt*k2den/2,
#                                         vol.+dt*k2vol/2,
#                                         val.+dt*k2val/2,
#                                         psi.+dt*k2psi/2,
#                                         tem.+dt*k2tem/2)
                                        
#     k4den,k4vol,k4val,k4psi,k4tem = RHS(den.+dt*k3den,
#                                         vol.+dt*k3vol,
#                                         val.+dt*k3val,
#                                         psi.+dt*k3psi,
#                                         tem.+dt*k3tem)

#     den = den .+ dt*(k1den .+2k2den .+2k3den .+k4den)/6
#     vol = vol .+ dt*(k1vol .+2k2vol .+2k3vol .+k4vol)/6
#     val = val .+ dt*(k1val .+2k2val .+2k3val .+k4val)/6
#     psi = psi .+ dt*(k1psi .+2k2psi .+2k3psi .+k4psi)/6
#     tem = tem .+ dt*(k1tem .+2k2tem .+2k3tem .+k4tem)/6


#     return den,vol,val,psi,tem

# end


# function RHS(den::Matrix,
#             vol::Matrix,
#             val::Matrix,
#             psi::Matrix,
#             tem::Matrix)

#     phi =  laplace_r(vol,dr,fkm,n0bc,nLbc,ar)
#     cur = -laplace(psi,fkm,dr,ar)

#     kden =( den0.*aln.*fkm.*phi
#             -den0.*fkp.*val
#             +fkp.*cur
#             -ϵ*(ss./q+(2.0 .-s).*s ./(ar.*q)).*fkm.*psi 
#             +visden.*laplace(den,fkm,dr,ar) )

#     kvol =( -tem0.*aln.*(1.0.+ηi).*fkm.*vol
#             +(1 ./den0).*(fkp.*cur-ϵ*(ss./q+(2.0 .-s).*s ./(ar.*q)).*fkm.*psi)
#             +visvol.*laplace(vol,fkm,dr,ar) )
    
#     kval =( -fkp.*(tem.+(1+τ).*tn0.*den)
#             -β.*tem0.*aln.*(1 .+τ .+ηi).*fkm.*psi
#             +visval.*laplace(val,fkm,dr,ar) ) 

#     kpsi =( -fkp.*(phi.-τ.*tn0.*den)
#             +β.*tem0.*τ.*aln.*fkm.*psi
#             +0.029.*τ^0.5.*abs.(fkp).*(val-cur./den0)
#             +vispsi.*laplace(psi,fkm,dr,ar) )./β

#     ktem =( tem0.*aln.*ηi.*fkm.*phi
#             -(Γ-1).*tem0.*fkp.*val
#             -(Γ-1).*(8 .*tem0/π).^0.5.*abs.(fkp).*tem
#             +vistem.*laplace(tem,fkm,dr,ar) )

#     return kden,kvol,kval,kpsi,ktem

# end