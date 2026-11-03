module solverBounded_mod
    use rigidBody_mod
    implicit none
    
    real(kind=rp), dimension(3, 3) :: ItensorTotal, ItensorTotalInv ! geometric moment of inertia tensor and location of point forces
    real(kind=rp), dimension(3) :: xCTotal

contains
    
subroutine setXCTotal(bodies, container)
    class(rigidBody_class), dimension(:) :: bodies
    class(container_class) :: container
    
    real(kind = rp) :: wgtdArea = 0.0_rp
    integer :: nBodies, iB

    nBodies = size(bodies, 1)

    xCTotal = (/0.0_rp, 0.0_rp, 0.0_rp/)
    do iB = 1, nBodies
        wgtdArea = wgtdArea + 4.0_rp*bodies(iB)%area
        xCTotal = xCTotal + 4.0_rp*bodies(iB)%area*bodies(iB)%xC
    end do
    wgtdArea = wgtdArea + container%area

    xCTotal = xCTotal + container%area*container%xC

    xCTotal = xCTotal/wgtdArea

end subroutine setXCTotal
    
subroutine setItensorTotal(bodies, container)
!! sets inertia tensor about xCTotal    
    class(rigidBody_class), dimension(:) :: bodies
    class(container_class) :: container
    real(kind=rp), dimension(3) :: delX
    integer :: iB, nBodies

    nBodies = size(bodies, 1)
    ItensorTotal = container%Itensor

    ! shift it to xCTotal:
    delX = container%xC - xCTotal
    ItensorTotal(1, :) = ItensorTotal(1, :) + (dot_product(delX, delX)*(/1.0, 0.0, 0.0/) &
                            - delX*delX(1))*container%area
    
    ItensorTotal(2, 2:3) = ItensorTotal(2, 2:3) + (dot_product(delX, delX)*(/1.0, 0.0/) &
                            - delX(2:3)*delX(2))*container%area
                            
    ItensorTotal(3, 3) = ItensorTotal(3, 3) + (dot_product(delX, delX) &
                            - delX(3)*delX(3))*container%area
    
    ! Now add weighted Itensors of the rigid bodies inside the box:
    do iB = 1, nBodies
        ItensorTotal = ItensorTotal + 4.0_rp*bodies(iB)%Itensor
        delX = bodies(iB)%xC - xCTotal
        ! shift ItensorBody to xCTotal:
        ItensorTotal(1, :) = ItensorTotal(1, :) + 4.0_rp*(dot_product(delX, delX)*(/1.0, 0.0, 0.0/) &
                            - delX*delX(1))*bodies(iB)%area
    
        ItensorTotal(2, 2:3) = ItensorTotal(2, 2:3) + 4.0_rp*(dot_product(delX, delX)*(/1.0, 0.0/) &
                                - delX(2:3)*delX(2))*bodies(iB)%area
                                
        ItensorTotal(3, 3) = ItensorTotal(3, 3) + 4.0_rp*(dot_product(delX, delX) &
                                - delX(3)*delX(3))*bodies(iB)%area

    end do

    ! use symmetry to assign other elements:
    ItensorTotal(2, 1) = ItensorTotal(1, 2)
    ItensorTotal(3, 1) = ItensorTotal(1, 3)
    ItensorTotal(3, 2) = ItensorTotal(2, 3)

    call inv3dMat(ItensorTotal, ItensorTotalInv)

end subroutine setItensorTotal
    
function VOmegaTotal(bodies, container)
!! required by the PR operator
    class(rigidBody_class), dimension(:) :: bodies
    class(container_class) :: container
    real(kind=rp), dimension(3) :: delX
    real(kind=rp), dimension(3, 2) :: VOmegaTotal, VIOmegaBody
    integer :: iB, nBodies
    real(kind=rp) :: wgtdArea

    nBodies = size(bodies, 1)
    wgtdArea = container%area

    VOmegaTotal = container%VIOmega(container%qDL)
    delX = container%xC - xCTotal
    VOmegaTotal(:, 2) = VOmegaTotal(:, 2) + cross(delX, VOmegaTotal(:, 1))

    ! add contributions from the bodies:
    do iB = 1, nBodies
        wgtdArea = wgtdArea + 4.0_rp*bodies(iB)%area
        ! for V:
        VIOmegaBody = bodies(iB)%VIOmega(bodies(iB)%qDL)
        VOmegaTotal(:, 1) = VOmegaTotal(:, 1) + 2.0_rp*VIOmegaBody(:, 1)    
        ! for omega:
        delX = bodies(iB)%xC - xCTotal
        VOmegaTotal(:, 2) = VOmegaTotal(:, 2) + 2.0_rp*VIOmegaBody(:, 2)
        VOmegaTotal(:, 2) = VOmegaTotal(:, 2) + 2.0_rp*cross(delX, VIOmegaBody(:, 1))
    end do
    VOmegaTotal(:, 1) = VOmegaTotal(:, 1)/wgtdArea

    VOmegaTotal(:, 2) = mat3dVec(ItensorTotalInv, VOmegaTotal(:, 2))
    
end function VOmegaTotal
    
function PR(bodies, container)
    class(rigidBody_class), dimension(:) :: bodies
    class(container_class) :: container
    integer :: iB, nBodies, cp, dimPR, ct 
    real(kind=rp), dimension(3) :: delX
    real(kind=rp), dimension(3, 2) :: vOmegaTmp
    ! real(kind=rp), dimension(this%nCollocationPts, 3) :: PR
    real(kind=rp), dimension(:, :), allocatable :: PR
    
    nBodies = size(bodies, 1)
    vOmegaTmp = VOmegaTotal(bodies, container)

    dimPR = container%nCollocationPts
    do iB = 1, nBodies
        dimPR = dimPR + bodies(iB)%nCollocationPts
    end do

    allocate(PR(dimPR, 3))

    ct = 1
    do iB = 1, nBodies
        do cp = 1, bodies(iB)%nCollocationPts
            delX = bodies(iB)%collocationPts(cp, :) - xCTotal
            PR(ct, :) = 2.0_rp*vOmegaTmp(:, 1) + 2.0_rp*cross(vOmegaTmp(:, 2), delX)
            ct = ct + 1
        end do
    end do

    do cp = 1, container%nCollocationPts
        delX = container%collocationPts(cp, :) - xCTotal
        PR(ct, :) = vOmegaTmp(:, 1) + cross(vOmegaTmp(:, 2), delX)
        ct = ct + 1
    end do

end function PR

function PLminusPnTotal(bodies, container)
    class(rigidBody_class), dimension(:) :: bodies
    class(container_class) :: container
    integer :: iB, nBodies, cp, dimPL, ct 
    real(kind=rp) :: totalFlux, totalArea, modY2, c1, c2
    ! real(kind=rp), dimension(this%nCollocationPts, 3) :: PR
    real(kind=rp), dimension(:, :), allocatable :: PLminusPnTotal
    real(kind=rp), dimension(:), allocatable :: qdotNs
    
    nBodies = size(bodies, 1)
    allocate(qdotNs(nBodies + 1))

    dimPL = container%nCollocationPts
    totalFlux = 0.0_rp
    totalArea = container%area
    do iB = 1, nBodies
        dimPL = dimPL + bodies(iB)%nCollocationPts
        qdotNs(iB) = bodies(iB)%qDotNhat(bodies(iB)%qDL)
        totalFlux = totalFlux + qdotNs(iB)
        totalArea = totalArea + bodies(iB)%area
    end do
    qdotNs(nBodies+1) = container%qDotNhat(container%qDL)
    totalFlux = totalFlux + qdotNs(nBodies+1)
    allocate(PLminusPnTotal(dimPL, 3))

    modY2 = 1.0_rp - container%area/totalArea
    c1 = (totalFlux/totalArea)*(2.0_rp + modY2) - qdotNs(nBodies+1)/totalArea 
    c2 = (1.0_rp/(totalArea*modY2))*( 1.0_rp/(1.0_rp + 3.0_rp*modY2) - 1.0_rp )
    c2 = c2*( totalFlux*(1.0_rp + modY2) - qdotNs(nBodies+1) )*(1.0_rp + modY2)
    
    ct = 1
    do iB = 1, nBodies
        do cp = 1, bodies(iB)%nCollocationPts
            
            PLminusPnTotal(ct, :) = bodies(iB)%normalVectors(cp, :)*(c1 + qdotNs(iB)/bodies(iB)%area + c2)     
            PLminusPnTotal(ct, :) = PLminusPnTotal(ct, :) - bodies(iB)%normalVectors(cp, :)*totalFlux/totalArea
            ct = ct + 1
        end do
    end do

    ! for the container:
    c1 = (totalFlux/totalArea)*(1.0_rp + modY2) - qdotNs(nBodies+1)/totalArea 
    c2 = c2*(modY2)/(1.0_rp + modY2)
    do cp = 1, container%nCollocationPts
        PLminusPnTotal(ct, :) = container%normalVectors(cp, :)*(c1 + c2)
        PLminusPnTotal(ct, :) = PLminusPnTotal(ct, :) - container%normalVectors(cp, :)*totalFlux/totalArea
        ct = ct + 1
    end do

end function PLminusPnTotal

function picardIterate(bodies, container) result(picardError)
    class(rigidBody_class), dimension(:) :: bodies
    class(container_class) :: container
    real(kind=rp) :: picardError, oneByNFT
    integer :: iB, nBodies, cp0, dimTotal, jB, ct, iFT
    real(kind=rp), dimension(:,:), allocatable :: PrbQ, cRHS, LdQ, PRQ, PLmPNTQ
    real(kind=rp), dimension(3) :: x0, delX
    logical :: chckSlf

    picardError = 0.0_rp
    nBodies = size(bodies, 1)

    dimTotal = container%nCollocationPts
    do iB = 1, nBodies
        dimTotal = dimTotal + bodies(iB)%nCollocationPts
    end do 

    allocate(PRQ(dimTotal, 3))
    allocate(PLmPNTQ(dimTotal, 3))
        
    ct = 1
    do iB = 1, nBodies
        ! reallocate opeartor arrays: 
        if(allocated(cRHS))  deallocate(cRHS)
        allocate(cRHS(bodies(iB)%nCollocationPts, 3))
        if(allocated(PrbQ))  deallocate(PrbQ)
        allocate(PrbQ(bodies(iB)%nCollocationPts, 3))
        if(allocated(LdQ))  deallocate(LdQ)
        allocate(LdQ(bodies(iB)%nCollocationPts, 3))

        PrbQ = bodies(iB)%Prb(bodies(iB)%qDL)
        PRQ = PR(bodies, container) ! recompute for every body because q(previousBodies) have been updated.
        PLmPNTQ = PLminusPnTotal(bodies, container)
        LdQ = 0.0_rp
        cRHS = 0.0_rp

        do cp0 = 1, bodies(iB)%nCollocationPts
        
            x0 = bodies(iB)%collocationPts(cp0, :) 
            ! compute Ld[qDL] and cRHS at each x0:
            do jB = 1, nBodies
                oneByNFT = (1.0_rp/bodies(jB)%nforcePt)
                do iFT = 1, bodies(jB)%nforcePt
                    delX = x0 - bodies(jB)%xFT(iFT, :)
                    cRHS(cp0, :) = cRHS(cp0, :) -  oneByNFT*( bodies(jB)%Fh/norm2(delX) &
                                    + delX*dot_product(bodies(jB)%Fh, delX)/(norm2(delX)**3) ) &
                                    + oneByNFT*cross(delX, bodies(jB)%Th)*(1.0_rp/(norm2(delX)**3))                    
                end do
    
                if ( jB == iB ) then
                    chckSlf = .true.
                else
                    chckSlf = .false.    
                end if
             
                LdQ(cp0, :) = LdQ(cp0, :) + bodies(jB)%LdX0(bodies(jB)%qDL, x0, cp0, chckSlf)
                
            end do
            ! contribution from container:
            LdQ(cp0, :) = LdQ(cp0, :) - container%LdX0(container%qDL, x0, cp0, .false.)
                
            picardError = picardError + norm2(cRHS(cp0, :) - LdQ(cp0, :) - PrbQ(cp0, :) &
                        + PLmPNTQ(ct-1 + cp0, :) + PRQ(ct-1 + cp0, :) )**2

        end do

        bodies(iB)%qDL = cRHS - LdQ - PrbQ + PLmPNTQ(ct:ct+bodies(iB)%nCollocationPts-1, :) &
                            + PRQ(ct:ct+bodies(iB)%nCollocationPts-1, :) + bodies(iB)%qDL   ! singularity subtraction at the end

        ct = ct + bodies(iB)%nCollocationPts
        
    end do

    ! for container:
    PRQ = PR(bodies, container) ! recompute for every body because q(previousBodies) have been updated.
    PLmPNTQ = PLminusPnTotal(bodies, container)
    if(allocated(cRHS))  deallocate(cRHS)
    allocate(cRHS(container%nCollocationPts, 3))
    if(allocated(LdQ))  deallocate(LdQ)
    allocate(LdQ(container%nCollocationPts, 3))

    cRHS = 0.0_rp
    LdQ = 0.0_rp

    do cp0 = 1, container%nCollocationPts
        
        x0 = container%collocationPts(cp0, :) 
        ! compute Ld[qDL] and cRHS at each x0:
        do jB = 1, nBodies
            oneByNFT = (1.0_rp/bodies(jB)%nforcePt)
            do iFT = 1, bodies(jB)%nforcePt
                delX = x0 - bodies(jB)%xFT(iFT, :)
                cRHS(cp0, :) = cRHS(cp0, :) -  oneByNFT*( bodies(jB)%Fh/norm2(delX) &
                                + delX*dot_product(bodies(jB)%Fh, delX)/(norm2(delX)**3) ) &
                                + oneByNFT*cross(delX, bodies(jB)%Th)*(1.0_rp/(norm2(delX)**3))                    
            end do
            ! add container velocity to cRHS if needed
            
            LdQ(cp0, :) = LdQ(cp0, :) - bodies(jB)%LdX0(bodies(jB)%qDL, x0, cp0, .false.)
            
        end do
        ! contribution from container:
        LdQ(cp0, :) = LdQ(cp0, :) + container%LdX0(container%qDL, x0, cp0, .true.)
            
        picardError = picardError + norm2( -cRHS(cp0, :) - LdQ(cp0, :) - 2.0_rp*container%qDL(cp0, :) &
                        + PLmPNTQ(ct-1 + cp0, :) + PRQ(ct-1 + cp0, :) )**2
        
    end do
    container%qDL = -cRHS - LdQ + PLmPNTQ(ct:ct+container%nCollocationPts-1, :) & ! supplementary flow is negative for container
                     + PRQ(ct:ct+container%nCollocationPts-1, :) - container%qDL   ! singularity subtraction at the end

    picardError =  sqrt(picardError)/dimTotal

    deallocate(PrbQ)
    deallocate(cRHS)
    deallocate(LdQ)
    deallocate(PRQ)
    deallocate(PLmPNTQ)

end function picardIterate

subroutine setqDLFinal(bodies, container)
!! call just after the picard iterates have converged
    class(rigidBody_class), dimension(:) :: bodies
    class(container_class) :: container
    integer :: nBodies, iB, ct, dimTotal
    real(kind=rp), dimension(:,:), allocatable :: PRQ

    nBodies = size(bodies, 1)
    
    dimTotal = container%nCollocationPts
    do iB = 1, nBodies
        dimTotal = dimTotal + bodies(iB)%nCollocationPts
    end do
    allocate(PRQ(dimTotal, 3))
    
    PRQ = PR(bodies, container)

    ct = 1
    do iB = 1, nBodies
        bodies(iB)%qDLFinal = bodies(iB)%qDL - 0.5_rp*PRQ(ct:ct-1+bodies(iB)%nCollocationPts, :) 
        ct = ct + bodies(iB)%nCollocationPts
    end do
    ! for container:
    container%qDLFinal = container%qDL - 0.5_rp*PRQ(ct:ct-1+container%nCollocationPts, :) 
    

end subroutine setqDLFinal


function getuField(bodies, container, x0) result(uX0)
    class(rigidBody_class), dimension(:) :: bodies
    class(container_class) :: container
    integer :: iter, nBodies, iB
    real(kind = rp), dimension(3) :: x0, uX0, uSupp, delX
    real(kind=rp) :: oneByNFT

    nBodies = size(bodies, 1)
    uX0 = 0.0_rp
    
    do iB = 1, nBodies
        oneByNFT = 1.0_rp/bodies(iB)%nforcePt 
        uSupp = 0.0_rp       
        do iter = 1, bodies(iB)%nforcePt
            delX = x0 - bodies(iB)%xFT(iter, :)
            uSupp = uSupp - oneByNFT*( bodies(iB)%Fh/norm2(delX) &
                            + delX*dot_product(bodies(iB)%Fh, delX)/(norm2(delX)**3) ) &
                            + oneByNFT*cross(delX, bodies(iB)%Th)*(1.0_rp/(norm2(delX)**3))                        
        end do

        uX0 = uX0 + uSupp - bodies(iB)%LdX0(bodies(iB)%qDL, x0, -1, .false.) ! cp0 doesn't matter for self=false case 
    end do
    ! add container contribution:
    
    uX0 = uX0 + container%LdX0(container%qDL, x0, -1, .false.) ! cp0 doesn't matter for self=false case 
    
end function getuField


end module solverBounded_mod