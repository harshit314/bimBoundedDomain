program tst
    use timeEvolveNew_mod
    ! use solverBounded_mod
    implicit none
    integer, parameter :: nBodies = 8 ! number of bodies
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
    character(len=:), allocatable :: fSv1(:), fSv2(:)
    real(kind=rp), dimension(12) :: dList
    integer :: iOuter
    
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

    call container%init('./meshData/boxExp10/boxCollocationPtsLessTst.txt', './meshData/boxExp10/boxNormalVectorsLessTst.txt'&
                        , './meshData/boxExp10/boxWgtsLessTst.txt', './meshData/boxExp10/boxVelLessTst.txt', .true.)
    
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

    ! create array for filenames:
    allocate(character(len=16) :: fSv1(12))
    allocate(character(len=16) :: fSv2(12))
    dList = [2.5_rp, 3.0_rp, 3.33_rp, 3.5_rp, 4.0_rp, 4.5_rp, 5.0_rp, 5.33_rp, 5.5_rp, 6.0_rp, 6.5_rp, 7.0_rp]

    fSv1 = [character(len=16) :: "mobMatL2p5", "mobMatL3", "mobMatL3p33", "mobMatL3p5", &
                "mobMatL4", "mobMatL4p5", "mobMatL5", "mobMatL5p33", "mobMatL5p5", "mobMatL6", "mobMatL6p5", "mobMatL7"]
    fSv2 = [character(len=16) :: "uFieldL2p5", "uFieldL3", "uFieldL3p33", "uFieldL3p5", &
                "uFieldL4", "uFieldL4p5", "uFieldL5", "uFieldL5p33", "uFieldL5p5", "uFieldL6", "uFieldL6p5", "uFieldL7"]
    
    ! call body(1)%translateTo( 0.0_rp*(/-aLattice/2.0_rp, y0, -aLattice/2.0_rp/) )
    do iOuter = 1, 12
        mobMat = 0.0_rp
        ! set up bodies:
        aLattice = dList(iOuter)
        y0 = -aLattice/2.0_rp
        call body(1)%translateTo( (/-aLattice/2.0_rp, y0, -aLattice/2.0_rp/) )
        call body(2)%translateTo( (/aLattice/2.0_rp, y0, -aLattice/2.0_rp/) )
        call body(3)%translateTo( (/aLattice/2.0_rp, y0, aLattice/2.0_rp/) )
        call body(4)%translateTo( (/-aLattice/2.0_rp, y0, aLattice/2.0_rp/) )

        call body(5)%translateTo( (/-aLattice/2.0_rp, y0+aLattice, -aLattice/2.0_rp/) )
        call body(6)%translateTo( (/aLattice/2.0_rp, y0+aLattice, -aLattice/2.0_rp/) )
        call body(7)%translateTo( (/aLattice/2.0_rp, y0+aLattice, aLattice/2.0_rp/) )
        call body(8)%translateTo( (/-aLattice/2.0_rp, y0+aLattice, aLattice/2.0_rp/) )
        
        call setXCTotal(body, container)
        call setItensorTotal(body, container)
        
        if(myRank == 0) print*, xCTotal, container%area
        
        if(myRank == 0)    print*, 'begin iterartes'

        !**********************************************************
        ! call setbRHS(body, container, bRHS)
        ! call gmres(body, container, bRHS, 100, 1.0e-5_rp)
        
        ! call setqDLFinal(body, container)

        ! do ix = 1, nBodies
        !     call body(ix)%updateVel()
        !     print*, body(ix)%vel, body(ix)%omega
        ! end do

        ! ******************************************************
        ! open file to write:
        if(myRank==0) open(unit=1, file="./outputData/mob8S/" // trim(fSv1(iOuter)) // ".txt")
        
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
            if(myRank==0)   print*, 'iBdy:', iBody, 'velocity part'
            body(iBody)%Fh = 0.0_rp
        end do

        do iBody = 1, nBodies
            do iter = 1, 3
                body(iBody)%Th = 0.0_rp
                body(iBody)%Th(iter) = 3.0_rp/4.0_rp 
                
                call setbRHS(body, container, bRHS)
                call gmres(body, container, bRHS, 100, 1.0e-5_rp)
                
                call setqDLFinal(body, container)
                
                do ix = 1, nBodies
                    call body(ix)%updateVel()
                    mobMat(3*(ix-1) + 1 : 3*(ix-1) + 3, 3*(iBody-1) + iter + 3*nBodies ) = body(ix)%vel
                    mobMat(3*(ix-1) + 1 + 3*nBodies : 3*(ix-1) + 3 + 3*nBodies, 3*(iBody-1) + iter + 3*nBodies ) = body(ix)%omega
                end do
                
            end do        
            if(myRank==0)   print*, 'iBdy:', iBody, 'omega part'
            body(iBody)%Th = 0.0_rp
        end do

        if ( myRank==0 ) then
            write(1,*) mobMat
            close(1)
        end if

        ! write velocity field in a file:
        if(myRank==0)   open(unit = 1, file = "./outputData/mob8S/" // trim(fSv2(iOuter)) // ".txt")
        ! for xy plane, x in [-24,24] and y in [-2*aLattice,2*aLattice]:
        dx = 0.2_rp
        dy = 0.2_rp
        do iy = 0, int(4.0_rp*aLattice/dy)
        do ix = 0, int(49.6_rp/dx)
                x0Tst = (/-24.8_rp + dx*ix, -2.0_rp*aLattice + dy*iy, 0.0_rp/)
                u0Tst = getuField(body, container, x0Tst) 
                if ( myRank == 0 ) then
                    write(1, *) x0Tst, u0Tst/body(1)%vel(2)
                end if
        
        end do 
        end do
        if(myRank==0)   close(1)
    
        if(myRank == 0)    print*, 'iOuter: ', iOuter
    end do
    endTym = MPI_Wtime()

    if(myRank == 0)    print*, 'time taken', endTym-startTym

    deallocate(workloads)
    call MPI_FINALIZE(ierror)

end program tst
