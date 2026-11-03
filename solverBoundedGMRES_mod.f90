module solverBoundedGMRES_mod
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

    ! shift it to xCTotal: CAUTION: CHECK IF THIS SHIFT IS CORRECT:
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
    
function VOmegaTotal(bodies, container, q0)
!! required by the PR operator
    class(rigidBody_class), dimension(:) :: bodies
    class(container_class) :: container
    real(kind=rp), dimension(3) :: delX
    real(kind=rp), dimension(:,:) :: q0
    real(kind=rp), dimension(3, 2) :: VOmegaTotal, VIOmegaBody
    integer :: iB, nBodies, ct
    real(kind=rp) :: wgtdArea

    nBodies = size(bodies, 1)
    wgtdArea = container%area

    VOmegaTotal = 0.0_rp
    
    ct = 0
    ! add contributions from the bodies:
    do iB = 1, nBodies
        wgtdArea = wgtdArea + 4.0_rp*bodies(iB)%area
        ! for V:
        VIOmegaBody = bodies(iB)%VIOmega( q0(ct+1 : ct+bodies(iB)%nCollocationPts, :) )
        VOmegaTotal(:, 1) = VOmegaTotal(:, 1) + 2.0_rp*VIOmegaBody(:, 1)    
        ! for omega:
        delX = bodies(iB)%xC - xCTotal
        VOmegaTotal(:, 2) = VOmegaTotal(:, 2) + 2.0_rp*VIOmegaBody(:, 2)
        VOmegaTotal(:, 2) = VOmegaTotal(:, 2) + 2.0_rp*cross(delX, VIOmegaBody(:, 1))
        ct = ct + bodies(iB)%nCollocationPts
    end do

    VOmegaTotal = VOmegaTotal + container%VIOmega( q0(ct+1 : ct+container%nCollocationPts, :) )
    delX = container%xC - xCTotal
    VOmegaTotal(:, 2) = VOmegaTotal(:, 2) + cross(delX, VOmegaTotal(:, 1))

    VOmegaTotal(:, 1) = VOmegaTotal(:, 1)/wgtdArea

    VOmegaTotal(:, 2) = mat3dVec(ItensorTotalInv, VOmegaTotal(:, 2))
    
end function VOmegaTotal
    
function PR(bodies, container, q0)
    class(rigidBody_class), dimension(:) :: bodies
    class(container_class) :: container
    integer :: iB, nBodies, cp, dimPR, ct 
    real(kind=rp), dimension(3) :: delX
    real(kind=rp), dimension(3, 2) :: vOmegaTmp
    real(kind=rp), dimension(:, :) :: q0
    real(kind=rp), dimension(:, :), allocatable :: PR
    
    nBodies = size(bodies, 1)
    vOmegaTmp = VOmegaTotal(bodies, container, q0)

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

function PLminusPnTotal(bodies, container, q0)
    class(rigidBody_class), dimension(:) :: bodies
    class(container_class) :: container
    integer :: iB, nBodies, cp, dimPL, ct, ct0 
    real(kind = rp), dimension(:, :) :: q0
    real(kind=rp) :: totalFlux, totalArea, modY2, c1, c2
    ! real(kind=rp), dimension(this%nCollocationPts, 3) :: PR
    real(kind=rp), dimension(:, :), allocatable :: PLminusPnTotal
    real(kind=rp), dimension(:), allocatable :: qdotNs
    
    nBodies = size(bodies, 1)
    allocate(qdotNs(nBodies + 1))

    dimPL = container%nCollocationPts
    totalFlux = 0.0_rp
    totalArea = container%area

    ct0 = 0
    do iB = 1, nBodies
        dimPL = dimPL + bodies(iB)%nCollocationPts
        qdotNs(iB) = bodies(iB)%qDotNhat( q0(ct0+1 : ct0+bodies(iB)%nCollocationPts, : ) )
        totalFlux = totalFlux + qdotNs(iB)
        totalArea = totalArea + bodies(iB)%area
        ct0 = ct0 + bodies(iB)%nCollocationPts
    end do
    qdotNs(nBodies+1) = container%qDotNhat( q0(ct0+1 : ct0+container%nCollocationPts, :) )
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

subroutine setqDLFinal(bodies, container)
!! call just after the iterates have converged
    class(rigidBody_class), dimension(:) :: bodies
    class(container_class) :: container
    integer :: nBodies, iB, ct, dimTotal
    real(kind=rp), dimension(:,:), allocatable :: PRQ, q0

    nBodies = size(bodies, 1)
    
    dimTotal = container%nCollocationPts
    do iB = 1, nBodies
        dimTotal = dimTotal + bodies(iB)%nCollocationPts
    end do
    allocate(PRQ(dimTotal, 3))
    allocate(q0(dimTotal, 3))
    
    ! ct = 0
    ! do iB = 1, nBodies
    !     q0(ct+1 : ct+bodies(iB)%nCollocationPts, :) = bodies(iB)%qDL
    !     ct = ct + bodies(iB)%nCollocationPts
    ! end do
    ! q0(ct+1 : ct+container%nCollocationPts, :) = container%qDL
    

    ! PRQ = PR(bodies, container, q0)

    ct = 1
    do iB = 1, nBodies
        bodies(iB)%qDLFinal = bodies(iB)%qDL! - 0.5_rp*PRQ(ct:ct-1+bodies(iB)%nCollocationPts, :) 
        ct = ct + bodies(iB)%nCollocationPts
    end do
    ! for container:
    container%qDLFinal = container%qDL !- 0.5_rp*PRQ(ct:ct-1+container%nCollocationPts, :) 
    
end subroutine setqDLFinal

function getuField(bodies, container, x0) result(uX0)
    class(rigidBody_class), dimension(:) :: bodies
    class(container_class) :: container
    integer :: iter, nBodies, iB, cp0
    real(kind = rp), dimension(3) :: x0, uX0, uSupp, delX
    real(kind=rp) :: oneByNFT, dMin, dTmp

    nBodies = size(bodies, 1)
    uX0 = 0.0_rp
    
    do iB = 1, nBodies
        oneByNFT = 1.0_rp/bodies(iB)%nforcePt 
        uSupp = 0.0_rp       
        do iter = 1, bodies(iB)%nforcePt
            delX = x0 - bodies(iB)%xFT(iter, :)
            if(norm2(delX) == 0.0_rp) cycle
            uSupp = uSupp - oneByNFT*( bodies(iB)%Fh/norm2(delX) &
                            + delX*dot_product(bodies(iB)%Fh, delX)/(norm2(delX)**3) ) &
                            + oneByNFT*cross(delX, bodies(iB)%Th)*(1.0_rp/(norm2(delX)**3))                        
        end do
        ! find cp0:
        dMin = norm2(x0 - bodies(iB)%collocationPts(1,:))
        cp0 = 1
        do iter = 2, bodies(iB)%nCollocationPts
            dTmp = norm2(x0 - bodies(iB)%collocationPts(iter, :))
            if(dMin >= dTmp)  then
                dMin = dTmp 
                cp0 = iter
            end if
        end do
        uX0 = uX0 + uSupp - bodies(iB)%LdX0(bodies(iB)%qDLFinal, x0, cp0, .true.) 
    end do
    
    ! add container contribution:
    ! find cp0:
    dMin = norm2(x0 - container%collocationPts(1,:))
    cp0 = 1
    do iter = 2, container%nCollocationPts
        dTmp = norm2(x0 - container%collocationPts(iter, :))
        if(dMin >= dTmp)  then
            dMin = dTmp 
            cp0 = iter
        end if
    end do
    uX0 = uX0 + container%LdX0(container%qDLFinal, x0, cp0, .true.)  &
                + 2.0_rp*container%qDL(cp0, :)   ! term from singularity subtraction for interior flow.
    
end function getuField

! For GMRES iterations:

subroutine setbRHS(bodies, container, bRHS)
    class(rigidBody_class), dimension(:) :: bodies
    class(container_class) :: container
    real(kind=rp) :: oneByNFT
    integer :: iB, nBodies, cp0, jB, ct, iFT
    real(kind=rp), dimension(:,:) :: bRHS
    real(kind=rp), dimension(3) :: x0, delX
    
    nBodies = size(bodies)
    bRHS = 0.0_rp

    ct = 0
    ! **** for bodies ****
    do iB = 1, nBodies
        do cp0 = 1, bodies(iB)%nCollocationPts
            x0 = bodies(iB)%collocationPts(cp0, :) 
            ! compute bRHS at each x0:
            do jB = 1, nBodies
                oneByNFT = (1.0_rp/bodies(jB)%nforcePt)   

                do iFT = 1, bodies(jB)%nforcePt
                
                    delX = x0 - bodies(jB)%xFT(iFT, :)
                    bRHS(cp0 + ct, :) = bRHS(cp0 + ct, :) -  oneByNFT*( bodies(jB)%Fh/norm2(delX) &
                                    + delX*dot_product(bodies(jB)%Fh, delX)/(norm2(delX)**3) ) &
                                    + oneByNFT*cross(delX, bodies(jB)%Th)*(1.0_rp/(norm2(delX)**3))    
                end do
            end do
        end do
        ct = ct + bodies(iB)%nCollocationPts
    end do

    ! **** for container ****
    do cp0 = 1, container%nCollocationPts
        x0 = container%collocationPts(cp0, :) 
        ! compute bRHS at each x0:
        do jB = 1, nBodies
            oneByNFT = (1.0_rp/bodies(jB)%nforcePt)   

            do iFT = 1, bodies(jB)%nforcePt
            
                delX = x0 - bodies(jB)%xFT(iFT, :)
                ! this is minus times the supplementary flow:
                bRHS(cp0 + ct, :) = bRHS(cp0 + ct, :) + oneByNFT*( bodies(jB)%Fh/norm2(delX) &
                                + delX*dot_product(bodies(jB)%Fh, delX)/(norm2(delX)**3) ) &
                                - oneByNFT*cross(delX, bodies(jB)%Th)*(1.0_rp/(norm2(delX)**3))    
            end do
        end do
        ! assign container velocity:
        ! bRHS(cp0 + ct, :) += container%uContainer(cp0, :)
    end do
    
end subroutine setbRHS

function AdotQ(bodies, container, q0)
    class(rigidBody_class), dimension(:) :: bodies
    class(container_class) :: container
    integer :: iB, nBodies, cp0, jB, ct, iFT, ctIn, cp0Tmp, iter
    real(kind=rp), dimension(:,:), allocatable :: PrbQ, LdQ, PnQ, PRQ, PLmPNTQ 
    real(kind=rp), dimension(3) :: x0, delX
    logical :: chckSlf
    real(kind=rp), dimension(:,:), allocatable :: AdotQ
    real(kind=rp), dimension(:,:) :: q0
    real(kind=rp) :: dMin, dTmp

    nBodies = size(bodies, 1)

    allocate(PRQ(size(q0, 1), 3))
    allocate(PLmPNTQ(size(q0, 1), 3))

    PRQ = PR(bodies, container, q0) 
    PLmPNTQ = PLminusPnTotal(bodies, container, q0)
            
    allocate(AdotQ(size(q0, 1), 3))

    ! ***** for bodies *****
    ct = 0
    do iB = 1, nBodies
        ! reallocate opeartor arrays: 
        if(allocated(PrbQ))  deallocate(PrbQ)
        allocate(PrbQ(bodies(iB)%nCollocationPts, 3))
        if(allocated(PnQ))  deallocate(PnQ)
        allocate(PnQ(bodies(iB)%nCollocationPts, 3))
        if(allocated(LdQ))  deallocate(LdQ)
        allocate(LdQ(bodies(iB)%nCollocationPts, 3))

        PrbQ = bodies(iB)%Prb(q0(ct+1:ct+bodies(iB)%nCollocationPts, :))
        ! PnQ = bodies(iB)%Pn(q0(ct+1:ct+bodies(iB)%nCollocationPts, :))
        LdQ = 0.0_rp
        
        do cp0 = 1, bodies(iB)%nCollocationPts
        
            x0 = bodies(iB)%collocationPts(cp0, :) 
            ! compute Ld[q0] at each x0:
            ctIn = 0
            do jB = 1, nBodies
                 
                if ( jB == iB ) then
                    chckSlf = .true.
                    cp0Tmp = cp0
                else
                    chckSlf = .false.    
                    
                    ! find cp0Tmp:
                    dMin = norm2(x0 - bodies(jB)%collocationPts(1,:))
                    cp0Tmp = 1
                    do iter = 2, bodies(jB)%nCollocationPts
                        dTmp = norm2(x0 - bodies(jB)%collocationPts(iter, :))
                        if(dMin >= dTmp)  then
                            dMin = dTmp 
                            cp0Tmp = iter
                        end if
                    end do
                
                end if
                
                LdQ(cp0, :) = LdQ(cp0, :) + bodies(jB)%LdX0(q0(ctIn+1 : ctIn+bodies(jB)%nCollocationPts,:), x0, cp0Tmp, .true.)

                ! LdQ(cp0, :) = LdQ(cp0, :) + bodies(jB)%LdX0(q0(ctIn+1 : ctIn+bodies(jB)%nCollocationPts,:), x0, cp0Tmp, chckSlf)

                ctIn = ctIn + bodies(jB)%nCollocationPts

            end do
            ! contribution from container:
            
            ! find cp0Tmp:
            dMin = norm2(x0 - container%collocationPts(1,:))
            cp0Tmp = 1
            do iter = 2, container%nCollocationPts
                dTmp = norm2(x0 - container%collocationPts(iter, :))
                if(dMin >= dTmp)  then
                    dMin = dTmp 
                    cp0Tmp = iter
                end if
            end do

            ! LdQ(cp0, :) = LdQ(cp0, :) - container%LdX0(q0(ctIn+1 : ctIn+container%nCollocationPts,:), x0, cp0Tmp, .true.) &
            !                 - 2.0_rp*q0(ctIn + cp0Tmp, :)   ! term from singularity subtraction for interior flow.

            LdQ(cp0, :) = LdQ(cp0, :) - container%LdX0(q0(ctIn+1 : ctIn+container%nCollocationPts,:), x0, cp0Tmp, .true.) &
                - 2.0_rp*q0(ctIn + cp0Tmp, :)   ! term from singularity subtraction for interior flow.
            
        end do

        ! AdotQ(ct+1 : ct+bodies(iB)%nCollocationPts, :) = LdQ + PrbQ !&
        !                !- PRQ(ct+1 : ct+bodies(iB)%nCollocationPts, :) - PnQ !- PLmPNTQ(ct+1 : ct+bodies(iB)%nCollocationPts, :) ! singularity subtraction cancels Identity times q term

                       
        AdotQ(ct+1 : ct+bodies(iB)%nCollocationPts, :) = LdQ + PrbQ !&
                    !    - PRQ(ct+1 : ct+bodies(iB)%nCollocationPts, :) - PLmPNTQ(ct+1 : ct+bodies(iB)%nCollocationPts, :) ! singularity subtraction cancels Identity times q term

        ct = ct + bodies(iB)%nCollocationPts
        
    end do   

    ! ***** for container *****
    ! reallocate opeartor arrays: 
    if(allocated(LdQ))  deallocate(LdQ)
    allocate(LdQ(container%nCollocationPts, 3))
    LdQ = 0.0_rp
        
    do cp0 = 1, container%nCollocationPts
    
        x0 = container%collocationPts(cp0, :) 
        ! compute Ld[q0] at each x0:
        ctIn = 0
        do jB = 1, nBodies
            ! find cp0Tmp:
            dMin = norm2(x0 - bodies(jB)%collocationPts(1,:))
            cp0Tmp = 1
            do iter = 2, bodies(jB)%nCollocationPts
                dTmp = norm2(x0 - bodies(jB)%collocationPts(iter, :))
                if(dMin >= dTmp)  then
                    dMin = dTmp 
                    cp0Tmp = iter
                end if
            end do

            LdQ(cp0, :) = LdQ(cp0, :) - bodies(jB)%LdX0(q0(ctIn+1 : ctIn+bodies(jB)%nCollocationPts,:), x0, cp0Tmp, .true.)

            ctIn = ctIn + bodies(jB)%nCollocationPts

        end do
        ! contribution from container:
        LdQ(cp0, :) = LdQ(cp0, :) + container%LdX0(q0(ctIn+1 : ctIn+container%nCollocationPts,:), x0, cp0, .true.)
        
    end do

    ! ct is currently the sum of collocation points of the bodies
    AdotQ(ct+1 : ct+container%nCollocationPts, :) = LdQ &
                    ! - PRQ(ct+1 : ct+container%nCollocationPts, :)  - PLmPNTQ(ct+1 : ct+container%nCollocationPts, :) &
                    + 2.0_rp*q0(ct+1 : ct+container%nCollocationPts, :) ! singularity subtraction cancels Identity times q term

    
end function AdotQ

subroutine arnoldIteration(bodies, container, Qmat, k, qCol, hCol)
    class(rigidBody_class), dimension(:) :: bodies
    class(container_class) :: container
    real(kind=rp), dimension(:,:,:) :: Qmat ! thrid dim is for the 3d vector
    integer :: k, iter ! k is the last assigned column number of Qmat
    real(kind=rp), dimension(:, :) :: qCol ! second dim is for the 3d vector.
    real(kind=rp), dimension(:) :: hCol 

    qCol = AdotQ(bodies, container, Qmat(:, k, :)) ! gives Aq_k

    do iter = 1, k
        hCol(iter) =  dot_product(qCol(:, 1), Qmat(:, iter, 1)) &
                    + dot_product(qCol(:, 2), Qmat(:, iter, 2)) + dot_product(qCol(:, 3), Qmat(:, iter, 3))
        
        qCol = qCol - hCol(iter)*Qmat(:, iter, :)
    end do

    hCol(k+1) = norm2(qCol) ! gives the correct norm even with 2d array
    qCol = qCol/hCol(k+1)

end subroutine arnoldIteration

function givensRotation(x1, x2, cn, sn) result(modX)
    real(kind=rp) :: x1, x2, cn, sn, modX

    modX = sqrt(x1**2 + x2**2)

    cn = x1/modX
    sn = x2/modX

end function givensRotation

subroutine gmres(bodies, container, b, maxIterations, threshold)
    !! b should not be a zero vector, xInit is taken as qDL of bodies
    class(rigidBody_class), dimension(:) :: bodies
    class(container_class) :: container
    real(kind=rp), dimension(:, :) :: b
    real(kind=rp), dimension(:, :), allocatable :: r ! residue
    integer :: maxIterations, iter, iB, nBodies, ct, nBest, kCol
    real(kind=rp) :: threshold, error, bNorm
    real(kind=rp), dimension(:,:,:), allocatable :: Qmat
    real(kind=rp), dimension(:,:), allocatable :: Hmat, xInit, omegaMat, RmatInv
    real(kind=rp), dimension(:), allocatable :: yHat ! last row and first column of omega matrix.
    real(kind=rp) :: cn, sn, rho, sigma, beta

    nBodies = size(bodies)
    allocate(Qmat(size(b, 1), maxIterations + 1, 3))
    allocate(Hmat(maxIterations+1, maxIterations))
    allocate(r(size(b, 1), 3))
    allocate(yHat(maxIterations))
    allocate(xInit(size(b, 1), 3))
    allocate(omegaMat(maxIterations+1, maxIterations+1))
    allocate(RmatInv(maxIterations, maxIterations))
    
    Qmat = 0.0_rp
    Hmat = 0.0_rp
    r = 0.0_rp
    yHat = 0.0_rp
    omegaMat = 0.0_rp
    RmatInv = 0.0_rp

    ! assign xInit:
    ct = 0
    do iB = 1, nBodies
        xInit(ct + 1 : ct + bodies(iB)%nCollocationPts, :) = bodies(iB)%qDL
        ct = ct + bodies(iB)%nCollocationPts
    end do
    xInit(ct + 1 : ct + container%nCollocationPts, :) = container%qDL
    
    r = b - AdotQ(bodies, container, xInit)
    bNorm = norm2(b)
    beta = norm2(r)
    error = beta/bNorm
    
    Qmat(:, 1, :) = r/norm2(r)
    
    omegaMat(1, 1) = 1.0_rp ! omega_0 = 1

    do iter = 1, maxIterations
        
        call arnoldIteration(bodies, container, Qmat, iter, Qmat(:, iter+1, :), Hmat(1:iter+1, iter))

        rho = dot_product(omegaMat(iter, 1:iter), Hmat(1:iter, iter)) ! using old omegaMat here...omega_(iter-1)
        sigma = Hmat(iter+1, iter)
        
        RmatInv(iter, iter) = givensRotation(rho, sigma, cn, sn)**(-1.0)

        if(iter>1) then
            RmatInv(1:iter-1, iter) = matmul(omegaMat(1:iter-1, 1:iter), Hmat(1:iter, iter))
            RmatInv(1:iter-1, iter) = -matmul( RmatInv(1:iter-1, 1:iter-1) , RmatInv(1:iter-1, iter) )*RmatInv(iter, iter)
        end if

        ! update omegaMatrix to omega_iter:
        do kCol = 1, iter
            omegaMat(iter : iter+1, kCol) = omegaMat(iter, kCol)*(/cn, -sn/)         
        end do
        omegaMat(iter : iter+1, iter+1) = (/sn, cn/) 

        error = abs(beta*omegaMat(iter+1, 1))/bNorm
        
        if(myRank==0)    print*, 'error:', error, iter

        nBest = iter

        if(error<threshold) exit

    end do

    if(nBest==maxIterations .and. myRank == 0) print*, 'maxIterations reached with error:', error, nBest

    yHat(1:nBest) = beta*omegaMat(1:nBest, 1) ! this is gamma...multiply by RmatInv now:
    yHat(1:nBest) = matmul(RmatInv(1:nBest, 1:nBest), yHat(1:nBest))

    xInit(:, 1) = xInit(:, 1) + matmul(Qmat(:, 1:nBest, 1), yHat(1:nBest))
    xInit(:, 2) = xInit(:, 2) + matmul(Qmat(:, 1:nBest, 2), yHat(1:nBest))
    xInit(:, 3) = xInit(:, 3) + matmul(Qmat(:, 1:nBest, 3), yHat(1:nBest))
    
    ! assign new qDL using xInit:
    ct = 0
    do iB = 1, nBodies
        bodies(iB)%qDL = xInit(ct + 1 : ct + bodies(iB)%nCollocationPts, :)
        ct = ct + bodies(iB)%nCollocationPts
    end do
    container%qDL = xInit(ct + 1 : ct + container%nCollocationPts, :)
    

    deallocate(Qmat)
    deallocate(Hmat)
    deallocate(r)
    deallocate(yHat)
    deallocate(xInit)
    deallocate(omegaMat)
    deallocate(RmatInv)
    
end subroutine gmres

end module solverBoundedGMRES_mod