program tst
    use timeEvolveNew_mod
    ! use solverBounded_mod
    implicit none
    integer, parameter :: nBodies = 1 ! number of bodies
    integer :: iBody, iter, dimTotal, ix, iy, iPos
    type(rigidBody_class), dimension(nBodies) :: body
    type(container_class) :: container
    real(kind=rp) :: startTym, endTym, error
    real(kind=rp), dimension(:,:), allocatable :: qTst, bRHS
    ! point force for supplementary flow:
    real(kind=rp), dimension(8, 3) :: xFpt ! at center of each sphere
    real(kind=rp), dimension(3) :: x0Tst, u0Tst
    real(kind=rp) :: dx, dy, cTocDist = 3.3333_rp, wgtRel
    character(2) :: fileIndx, nBWrt

    ! setup mpi variables:
    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, worldSize, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierror) ! rank starts from zero, as in C++!
    ! allocate workloads:
    allocate(workloads(worldSize))

    open(unit = 1, file = './meshData/sccStructRBy6/FvaL3p33.txt', status = 'old', action = 'read')
    read(1,*) wgtRel
    close(1)

    ! set force points on centre of each sphere in scc structure:
    ! cTocDist = 3.0_rp

    if(myRank == 0) print*, 'Fval:', wgtRel

    xFpt(1, :) = cTocDist*(/0.0_rp, 0.0_rp, 0.0_rp/)
    xFpt(2, :) = cTocDist*(/1.0_rp, 0.0_rp, 0.0_rp/)
    xFpt(3, :) = cTocDist*(/0.0_rp, 1.0_rp, 0.0_rp/)
    xFpt(4, :) = cTocDist*(/0.0_rp, 0.0_rp, 1.0_rp/)
    xFpt(5, :) = cTocDist*(/1.0_rp, 1.0_rp, 0.0_rp/)
    xFpt(6, :) = cTocDist*(/0.0_rp, 1.0_rp, 1.0_rp/)
    xFpt(7, :) = cTocDist*(/1.0_rp, 0.0_rp, 1.0_rp/)
    xFpt(8, :) = cTocDist*(/1.0_rp, 1.0_rp, 1.0_rp/)
    
    do iBody = 1, nBodies
        call body(iBody)%init('./meshData/sccStructRBy6/sccCPL3p33.txt', './meshData/sccStructRBy6/sccNVecL3p33.txt'&
                                , './meshData/sccStructRBy6/sccWgtL3p33.txt', 'redundant.txt', .false.)
        
        ! xFpt(1,:) = body(iBody)%xC
        ! xFpt(1,:) = 0.0_rp*body(iBody)%xC
        call body(iBody)%setSupplementaryFlow(xFpt) ! set supplementary flow before rotating!
        call body(iBody)%setOrientation((/1.0_rp, 0.0_rp, 0.0_rp/), (/0.0_rp, 1.0_rp, 0.0_rp/), (/0.0_rp, 0.0_rp, 1.0_rp/))
        

        ! call body(iBody)%rotate((/0.0_rp, 1.0_rp, 0.0_rp/), pi/2.0_rp)

        body(iBody)%Fh = wgtRel*(3.0_rp/4.0_rp)*(/0.0_rp, 1.0_rp, 0.0_rp/) ! this is defined to be divided by 8piMu 
        body(iBody)%Th = 1.0_rp*(3.0_rp/4.0_rp)*(/0.0_rp, 0.0_rp, 0.0_rp/) ! this is defined to be divided by 8piMu 

    end do

    call container%init('./meshData/boxExp10/boxCollocationPtsLessTst.txt', './meshData/boxExp10/boxNormalVectorsLessTst.txt'&
                        , './meshData/boxExp10/boxWgtsLessTst.txt', './meshData/boxExp10/boxVelLessTst.txt', .true.)
    
    ! calculate dimTotal:
    dimTotal = container%nCollocationPts
    do iBody = 1, nBodies
        dimTotal = dimTotal + body(iBody)%nCollocationPts
    end do

    ! allocate bRHS:
    allocate(bRHS(dimTotal, 3))

    if(myRank == 0) print*, body(1)%xC

    ! call body(1)%translateTo( (/ 0.0_rp, 18.0_rp, 0.0_rp /) )

    ! call body(2)%translateTo( (/ 1.0783_rp, -2.81559_rp, -0.175302_rp /) )
    ! do iBody = 1, nBodies
    !     call body(iBody)%translateTo( (/ (iBody-nBodies/2 - 0.5)*2.0_rp, (iBody-nBodies/2 - 0.5)**2 * 2.0_rp, 0.0_rp /) )
    ! end do
    
    call body(1)%translateTo( (/ 0.0_rp, 0.0_rp, 0.0_rp /) )
    ! call body(2)%translateTo( (/ 2.0_rp, 1.0_rp, 0.0_rp /) )
                    
    if(myRank == 0) then
        do iter = 1, nBodies        
            write(fileIndx, '(I2)') iter
            call body(iter)%writeMesh('./outputData/sccCrystal/mesh/meshOut'//trim(adjustl(fileIndx))//'.txt')            
        end do       
    endif

    if(myRank == 0) print*, container%xC

    if(myRank == 0) print*, body(1)%area
    ! if(myRank == 0) print*, container%area

    startTym = MPI_Wtime()

    call setXCTotal(body, container)
    call setItensorTotal(body, container)
        
    if(myRank == 0) print*, xCTotal
    
    if(myRank == 0)    print*, 'begin iterartes'
    
    call setbRHS(body, container, bRHS) 
    call gmres(body, container, bRHS, 350, 1.0e-4_rp) ! did 10^-3 before
    
    call setqDLFinal(body, container)
    
    do iter = 1, nBodies
        call body(iter)%updateVel()
        if ( myRank == 0 ) then
            print*, body(iter)%vel, body(iter)%omega
        end if
    end do
    ! call tEvolve(body, container, 0.0_rp, 75.0_rp, './outputData/sphereWall.txt', './outputData/sphereWallVel.txt')
    
    ! ***************************************************
    ! save in a file:
    if(myRank==0) then 
        open(unit = 1, file = "./outputData/sccStructRBy6/velL3p33.txt")
        write(1, *) body(1)%vel
        write(1, *) body(1)%omega
        close(1)
    end if

    ! write velocity field in a file:
    if(myRank==0)   open(unit = 1, file = "./outputData/sccStructRBy6/uFieldL3p33.txt")
    ! for xy plane, x in [-24,24] and y in [-2*aLattice,2*aLattice]:
    dx = 0.2_rp
    dy = 0.2_rp
    do iy = 0, int(4.0_rp*cTocDist/dy)
       do ix = 0, int(49.6_rp/dx)
            x0Tst = (/-24.8_rp + dx*ix, -2.0_rp*cTocDist + dy*iy, 0.0_rp/)
            u0Tst = getuField(body, container, x0Tst) 
            if ( myRank == 0 ) then
                write(1, *) x0Tst, u0Tst/body(1)%vel(2)
            end if
    
       end do 
    end do
    if(myRank==0)   close(1)

    endTym = MPI_Wtime()

    if(myRank == 0)    print*, 'time taken', endTym-startTym

    deallocate(workloads)
    call MPI_FINALIZE(ierror)

end program tst
