module timeEvolveNew_mod
    use solverBoundedGMRES_mod
    implicit none
    real(kind=rp) :: dt = 0.0001_rp, dtMin = 0.0001_rp, dtMax = 0.2_rp, dtTolerance = 1.0e-4_rp
    real(kind=rp) :: picardTolerance = 1.0e-8_rp
    integer, parameter :: iterMax = 120 
    private :: setqDL

contains
        
    subroutine setqDL(bodies, container)
    ! set the double layer densities after Picard iterates
        class(container_class) :: container
        class(rigidBody_class), dimension(:) :: bodies
        integer :: nBodies, iter, dim
        real(kind=rp), dimension(:,:), allocatable :: bRHS
        real(kind=rp) :: picardError
        ! real(kind=rp), external :: solver

        nBodies = size(bodies, 1)
        dim = container%nCollocationPts
        do iter = 1, nBodies
            dim = dim + bodies(iter)%nCollocationPts
        end do
        allocate(bRHS(dim, 3))

        call setbRHS(bodies, container, bRHS)
        call gmres(bodies, container, bRHS, iterMax, picardTolerance)
        
        ! do iter = 1, iterMax
        !     picardError = 0.0_rp
        !     ! print*, 'iteration no. ', iter
            
        !     picardError = picardIterate(bodies)

        !     if(picardError <= picardTolerance ) then
        !         print*, 'error after iteration', iter, 'is', picardError
        !         exit
        !     endif
        !     if(iter == iterMax .and. myRank == 0) print*, 'picard iterates didnt converge, picardError:', picardError

        ! end do

    end subroutine setqDL

    function yDot(y, bodies, container)
        real(kind=rp), dimension(:, :) :: y
        real(kind=rp), dimension(:, :), allocatable :: yDot
        class(rigidBody_class), dimension(:) :: bodies
        class(container_class) :: container
        integer :: nBodies, iB

        allocate(yDot(size(y, 1), size(y, 2)))
        nBodies = size(bodies)

        do iB = 1, nBodies
            call bodies(iB)%translate(y(iB, 1:3))     
            call bodies(iB)%rotate(y(iB, 4:6), 1.0_rp)  
        end do

        call setqDL(bodies, container)

        do iB = 1, nBodies
            call bodies(iB)%updateVel()
            ! bring bodies back to previous config:
            call bodies(iB)%translate(-y(iB, 1:3))     
            call bodies(iB)%rotate(-y(iB, 4:6), 1.0_rp)

            yDot(iB, 1:3) = bodies(iB)%vel
            yDot(iB, 4:6) = bodies(iB)%omega            
        end do
        
    end function

    function RKF4Step(t, h, bodies, container) result(hNew)
        !! advance time by h and update the positions and orientations of bodies
            class(rigidBody_class), dimension(:) :: bodies
            class(container_class) :: container
            real(kind=rp), dimension(:, :), allocatable :: k1, k2, k3, k4, k5, k6, yVel
            integer :: dim, iB, nBodies
            real(kind=rp) :: t, h, truncError, hNew
            
            dim = 6 ! 6 DOF per body
            nBodies = size(bodies)
    
            allocate(k1(nBodies, dim))
            allocate(k2(nBodies, dim))
            allocate(k3(nBodies, dim))
            allocate(k4(nBodies, dim))
            allocate(k5(nBodies, dim))
            allocate(k6(nBodies, dim))
            allocate(yVel(nBodies, dim))
    
            ! set y:
            ! do iB = 1, nBodies
            !     y(iB, 1:3) = bodies(iB)%xC 
            !     y(iB, 4:6) = bodies(iB)%orient
            ! end do
            ! y=0 ensures calculating yDot at previous configuration 
            
            ! assign k1:
            call setqDL(bodies, container)
            do iB = 1, nBodies
                call bodies(iB)%updateVel()
                k1(iB, 1:3) = bodies(iB)%vel
                k1(iB, 4:6) = bodies(iB)%omega            
            end do
            ! assign k2:
            k2 = yDot( (2.0d-0/9.0d-0)*h*k1, bodies, container)
            ! assign k3:
            k3 = yDot( (1.0d-0/12.0d-0)*h*k1 + (1.0d-0/4.0d-0)*h*k2, bodies, container)
            ! assign k4:
            k4 = yDot( (69.0d-0/128.0d-0)*h*k1 + (-243.0d-0/128.0d-0)*h*k2 + (135.0d-0/64.0d-0)*h*k3, bodies, container)
            ! assign k5:
            k5 = yDot( (-17.0d-0/12.0d-0)*h*k1 + (27.0d-0/4.0d-0)*h*k2 + (-27.0d-0/5.0d-0)*h*k3 &
                        + (16.0d-0/15.0d-0)*h*k4, bodies, container)
            ! assign k6:
            k6 = yDot( (65.0d-0/432.0d-0)*h*k1 + (-5.0d-0/16.0d-0)*h*k2 + (13.0d-0/16.0d-0)*h*k3 &
                        + (4.0d-0/27.0d-0)*h*k4 + (5.0d-0/144.0d-0)*h*k5, bodies, container)
            ! get yNext:
    
            yVel = ((47.0d-0/450.0d-0)*k1 + (12.0d-0/25.0d-0)*k3 + (32.0d-0/225.0d-0)*k4 &
                            + (1.0d-0/30.0d-0)*k5 + (6.0d-0/25.0d-0)*k6)
            
            ! truncation error:
            truncError = h*norm2((1.0d-0/150.0d-0)*k1 + (-3.0d-0/100.0d-0)*k3 + (16.0d-0/75.0d-0)*k4 &
                                    + (1.0d-0/20.0d-0)*k5 + (-6.0d-0/25.0d-0)*k6)

            hNew = 0.9d-0*h*(dtTolerance/truncError)**(1.0d-0/5.0d-0)

            if(truncError < dtTolerance) then
                ! update the bodies position and orientation:
                do iB = 1, nBodies
                    call bodies(iB)%translate(h*yVel(iB, 1:3))   
                    call bodies(iB)%rotate(yVel(iB, 4:6), h)
                    bodies(iB)%vel = yVel(iB,1:3)
                    bodies(iB)%omega = yVel(iB,4:6)
                end do
                
                if(myRank==0) then
                    write(100, *) t+h, (bodies(iB)%xC, iB = 1, nBodies), (bodies(iB)%pz, iB = 1, nBodies)
                    write(101, *) t+h, (bodies(iB)%vel, iB = 1, nBodies), (bodies(iB)%omega, iB = 1, nBodies)
                    call flush(100)
                    call flush(101)
                end if

                ! increase t by h:
                if(hNew>=dtMax) hNew = dtMax
                t = t + h
            end if
    end function RKF4Step

    subroutine tEvolveAndWrite(tI, tF, bodies, container, filePos, fileVel)
    !! evolve the particles from tI to tF and write position at each time step.
        class(rigidBody_class), dimension(:) :: bodies
        class(container_class) :: container
        real(kind=rp) :: tI, tF, tCurr, dtCurr
        character(*) :: filePos, fileVel
        integer :: iB, nBodies

        tCurr = tI
        nBodies = size(bodies) 

        ! write data to a file:
        if(myRank==0) then
            open(unit=100, file=filePos)
            open(unit=101, file=fileVel)    
            ! write(100, *) tCurr, (bodies(iB)%xC, iB = 1, nBodies), (bodies(iB)%orient, iB = 1, nBodies)
            ! write(101, *) tCurr, (bodies(iB)%vel, iB = 1, nBodies), (bodies(iB)%omega, iB = 1, nBodies)
        end if 

        dtCurr = dt
        do while(tCurr<=tF)
            
            ! call RK4Step(tCurr, bodies)
            
            dtCurr = RKF4Step(tCurr, dtCurr, bodies, container)
            
            if(myRank==0)    print*, 'tCurr:', tCurr, 'dtNext:', dtCurr
            
            if(tCurr+dtCurr>tF) then
                dtCurr = RKF4Step(tCurr, tF-tCurr, bodies, container)
                exit
            end if

        end do
        
        if(myRank==0) then
            close(100)
            close(101)
        end if
    
    end subroutine tEvolveAndWrite

end module timeEvolveNew_mod