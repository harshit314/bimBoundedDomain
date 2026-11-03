program tst
    use timeEvolveNew_mod
    ! use solverBounded_mod
    implicit none
    integer, parameter :: nBodies = 14 ! number of bodies 8 + 6
    integer :: iBody, iter, dimTotal, ix, iy
    type(rigidBody_class), dimension(nBodies) :: body
    type(container_class) :: container
    real(kind=rp) :: startTym, endTym, error
    real(kind=rp), dimension(:,:), allocatable :: qTst, bRHS
    ! point force for supplementary flow:
    real(kind=rp), dimension(1, 3) :: xFpt
    real(kind=rp), dimension(3) :: x0Tst, u0Tst
    real(kind=rp) :: dx, dy, aLattice, y0
    character(2) :: fileIndx, nBWrt
    real(kind=rp), dimension(:, :), allocatable :: mobMat

    ! setup mpi variables:
    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, worldSize, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierror) ! rank starts from zero, as in C++!
    ! allocate workloads:
    allocate(workloads(worldSize))

    do iBody = 1, nBodies
        call body(iBody)%init('./meshData/sphereCPLess.txt', './meshData/sphereNVecLess.txt'&
                                , './meshData/sphereWgtsLess.txt', 'redundantLess.txt', .false.)
        
        xFpt(1,:) = body(iBody)%xC
        call body(iBody)%setSupplementaryFlow(xFpt) ! set supplementary flow before rotating!
                                                
        call body(iBody)%setOrientation((/1.0_rp, 0.0_rp, 0.0_rp/), (/0.0_rp, 1.0_rp, 0.0_rp/), (/0.0_rp, 0.0_rp, 1.0_rp/))
        
        ! call body(iBody)%rotate((/0.0_rp, 1.0_rp, 0.0_rp/), pi/2.0_rp)

        body(iBody)%Fh = 0.0_rp*(3.0_rp/4.0_rp)*(/0.0_rp, 1.0_rp, 0.0_rp/) ! this is defined to be divided by 8piMu 
        body(iBody)%Th = 0.0_rp*(3.0_rp/4.0_rp)*(/0.0_rp, 0.0_rp, 0.0_rp/) ! this is defined to be divided by 8piMu 

    end do

    call container%init('./meshData/boxExp10/boxCollocationPtsLess.txt', './meshData/boxExp10/boxNormalVectorsLess.txt'&
                        , './meshData/boxExp10/boxWgtsLess.txt', './meshData/boxExp10/boxVelLess.txt', .true.)
    
    ! calculate dimTotal:
    dimTotal = container%nCollocationPts
    do iBody = 1, nBodies
        dimTotal = dimTotal + body(iBody)%nCollocationPts
    end do

    ! allocate bRHS:
    allocate(bRHS(dimTotal, 3))

    allocate(mobMat(nBodies*6, nBodies*6))
    mobMat = 0.0_rp

    startTym = MPI_Wtime()

    ! set up bodies:
    aLattice = 6_rp
    y0 = -aLattice/2.0_rp
    call body(1)%translateTo( (/-aLattice/2.0_rp, y0, -aLattice/2.0_rp/) )
    call body(2)%translateTo( (/aLattice/2.0_rp, y0, -aLattice/2.0_rp/) )
    call body(3)%translateTo( (/aLattice/2.0_rp, y0, aLattice/2.0_rp/) )
    call body(4)%translateTo( (/-aLattice/2.0_rp, y0, aLattice/2.0_rp/) )

    call body(5)%translateTo( (/-aLattice/2.0_rp, y0+aLattice, -aLattice/2.0_rp/) )
    call body(6)%translateTo( (/aLattice/2.0_rp, y0+aLattice, -aLattice/2.0_rp/) )
    call body(7)%translateTo( (/aLattice/2.0_rp, y0+aLattice, aLattice/2.0_rp/) )
    call body(8)%translateTo( (/-aLattice/2.0_rp, y0+aLattice, aLattice/2.0_rp/) )
    ! put on face centres:
    call body(9)%translateTo( (/-aLattice/2.0_rp, 0.0_rp, 0.0_rp/) )
    call body(10)%translateTo( (/aLattice/2.0_rp, 0.0_rp, 0.0_rp/) )
 
    call body(11)%translateTo( (/0.0_rp, y0, 0.0_rp/) )
    call body(12)%translateTo( (/0.0_rp, y0+aLattice, 0.0_rp/) )

    call body(13)%translateTo( (/0.0_rp, 0.0_rp, y0/) )
    call body(14)%translateTo( (/0.0_rp, 0.0_rp, y0+aLattice/) )


    call setXCTotal(body, container)
    call setItensorTotal(body, container)
    
    if(myRank == 0) then
        do iter = 1, nBodies        
            write(fileIndx, '(I2)') iter
            call body(iter)%writeMesh('./outputData/fccCrystal/mesh/meshOut'//trim(adjustl(fileIndx))//'.txt')            
        end do       
    endif
    

    if(myRank == 0) print*, xCTotal
    
    if(myRank == 0)    print*, 'begin iterartes'
    
    ! open file to write:
    if(myRank==0) open(unit=1, file="./outputData/mobMatL5p33.txt")
    
    do iBody = 1, nBodies
        do iter = 1, 3
            body(iBody)%Fh = 0.0_rp
            body(iBody)%Fh(iter) = 3.0_rp/4.0_rp 
            
            call setbRHS(body, container, bRHS)
            call gmres(body, container, bRHS, 100, 1.0e-5_rp)
            
            call setqDLFinal(body, container)
            
            do ix = 1, nBodies
                call body(ix)%updateVel()
                mobMat(3*(ix-1) + 1 : 3*(ix-1) + 3, 3*(iBody-1) + iter ) = body(ix)%vel
                mobMat(3*(ix-1) + 1 + 3*nBodies : 3*(ix-1) + 3 + 3*nBodies, 3*(iBody-1) + iter ) = body(ix)%omega
            end do
        end do
        if(myRank==0)   print*, 'iBdy:', iBody
        body(iBody)%Fh = 0.0_rp
    end do

    do iBody = 1, nBodies
        do iter = 1, 3
            body(iBody)%Th = 0.0_rp
            body(iBody)%Th(iter) = 3.0_rp/4.0_rp 
            
            call setbRHS(body, container, bRHS)
            call gmres(body, container, bRHS, 100, 1.0e-6_rp)
            
            call setqDLFinal(body, container)
            
            do ix = 1, nBodies
                call body(ix)%updateVel()
                mobMat(3*(ix-1) + 1 : 3*(ix-1) + 3, 3*(iBody-1) + iter + 3*nBodies ) = body(ix)%vel
                mobMat(3*(ix-1) + 1 + 3*nBodies : 3*(ix-1) + 3 + 3*nBodies, 3*(iBody-1) + iter + 3*nBodies ) = body(ix)%omega
            end do
            
        end do        
        if(myRank==0)   print*, 'iBdy:', iBody
        body(iBody)%Th = 0.0_rp
    end do

    if ( myRank==0 ) then
        write(1,*) mobMat
        close(1)
    end if

    ! ! write velocity field in a file:
    ! if(myRank==0)   open(unit = 1, file = "./outputData/uField.txt")
    ! ! for xy plane, x in [-2.5,2.5] and y in [-2.5,2.5]:
    ! dx = 0.02_rp
    ! dy = 0.02_rp
    ! do iy = 0, int(4.0_rp/dy)
    !    do ix = 0, int(4.0_rp/dx)
    !         x0Tst = (/-2.0_rp + dx*ix, -2.0_rp + dy*iy, 0.0_rp/)
    !         u0Tst = getuField(body, container, x0Tst) 
    !         if ( myRank == 0 ) then
    !             write(1, *) x0Tst, u0Tst - body(1)%vel
    !         end if
    
    !    end do 
    ! end do
    ! if(myRank==0)   close(1)

    ! u0Tst = getuField(body, container, 1.00_rp*(/cos(pi/4.0_rp), sin(pi/4.0_rp), 0.0_rp/)) 
            
    ! if(myRank==0)   print*, 'uF:', u0Tst - body(1)%vel


    endTym = MPI_Wtime()

    if(myRank == 0)    print*, 'time taken', endTym-startTym

    deallocate(workloads)
    call MPI_FINALIZE(ierror)

end program tst
