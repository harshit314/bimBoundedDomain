module container_mod
    use matrix_mod
    ! use mpi_f08
    implicit none
    include "mpif.h"
    ! define container class:
    type :: container_class

        real(kind=rp), dimension(:,:), allocatable :: normalVectors, collocationPts, uContainer ! normalVectors(i,1:3) is the normalized normal vector of ith element.
        real(kind=rp), dimension(:), allocatable :: wgts ! wgts associated with the quadrature
        integer :: nCollocationPts
        integer :: myStartIndx, myEndIndx ! for mpi
        
        real(kind=rp), dimension(3) :: xC
        real(kind=rp) :: area
        real(kind=rp), dimension(3, 3) :: Itensor, ItensorInv ! geometric moment of inertia tensor and location of point forces
        real(kind=rp), dimension(:,:), allocatable :: qDL, qDLFinal ! define a double layer density on the surface. qDL(k) is a 3d vector corresponding to collocation point k.
    
    contains

        procedure  :: setMeshData

        procedure :: init
        procedure :: assignWorkload

        procedure :: setArea
        procedure :: setXC
        procedure :: setItensor
        
        procedure :: writeMesh
        ! BIM related functions:
        procedure :: VIOmega
        procedure :: qDotNhat
    
        ! procedure :: Prb
        procedure :: Pn
        procedure :: LdX0
        
    end type container_class

contains

    subroutine init(this, cpFile, normalVecFile, wgtsFile, velFile, isContainer)
    !! initializes all the data types defined in the class
        class(container_class) :: this
        character(*) :: cpFile, normalVecFile, wgtsFile, velFile
        logical :: isContainer

        call this%setMeshData(cpFile, normalVecFile, wgtsFile, velFile, isContainer)
        call this%setArea()
        call this%setXC()
        call this%setItensor()
        ! allocate the double layer density:
        allocate(this%qDL(this%nCollocationPts, 3))
        this%qDL = 0.0_rp
        allocate(this%qDLFinal(this%nCollocationPts, 3))
        this%qDLFinal = 0.0_rp

        ! for mpi:
        call this%assignWorkload()

    end subroutine init

    subroutine assignWorkload(this)
    !! divides the collocation points on each node. This evaluates integrals faster.
        class(container_class) :: this
        integer :: i
        
        do i = 1, worldSize
            workloads(i) = int(this%nCollocationPts/worldSize) ! defined in global module
            if( i-1 < mod(this%nCollocationPts, worldSize) ) workloads(i) = workloads(i) + 1
        end do
        
        this%myStartIndx = 1
        do i = 1, myRank ! myRank starts from 0
            this%myStartIndx = this%myStartIndx + workloads(i)
        end do
        this%myEndIndx = this%myStartIndx + (workloads(myRank+1) - 1)

        ! print*, 'myRank', myRank, 'myStart', this%myStartIndx, 'myEnd', this%myEndIndx

    end subroutine assignWorkload

    subroutine setMeshData(this, cpFile, normalVecFile, wgtsFile, velFile, isContainer)
    !! set the collocation points, normalVectors and quadrature points
        class(container_class) :: this
        integer :: io
        logical :: isContainer
        character(*) :: cpFile, normalVecFile, wgtsFile, velFile
        open(unit = 1, file = cpFile, status = 'old', action = 'read', iostat=io)
        if(io /= 0 .and. myRank == 0) print*, 'Error reading collocationPts file'  
        open(unit = 2, file = normalVecFile, status = 'old', action = 'read', iostat=io) 
        if(io /= 0 .and. myRank == 0) print*, 'Error reading normalVector file'  
        open(unit = 3, file = wgtsFile, status = 'old', action = 'read', iostat=io) 
        if(io /= 0 .and. myRank == 0) print*, 'Error reading weights file'
        if ( isContainer ) then
            open(unit = 4, file = velFile, status = 'old', action = 'read', iostat=io) 
            if(io /= 0 .and. myRank == 0) print*, 'Error reading velocity file'        
        end if
        
        ! allocate collocationPts
        this%nCollocationPts = 0
        do 
            read(1, *, iostat=io) 
            if ( io /= 0 ) exit
            this%nCollocationPts = this%nCollocationPts + 1
        end do
        allocate(this%collocationPts(3, this%nCollocationPts))! fortran reads column wise. Need to take transpose later
        rewind(1)
        
        ! allocate normalVectors and weights
        allocate(this%normalVectors(3, this%nCollocationPts))! fortran reads column wise. Need to take transpose later
        allocate(this%wgts(this%nCollocationPts))
        allocate(this%uContainer(3, this%nCollocationPts))! fortran reads column wise. Need to take transpose later
        this%uContainer = 0.0_rp

        read(1, *) this%collocationPts
        this%collocationPts = transpose(this%collocationPts)
        read(2, *) this%normalVectors
        this%normalVectors = transpose(this%normalVectors)

        ! read quadrature weights:
        read(3, *) this%wgts  

        ! read velocity file:
        if ( isContainer ) then
            read(4, *) this%uContainer
            this%uContainer = transpose(this%uContainer)            
        end if

        if(myRank == 0)    print*, 'Collocation points loaded:', this%nCollocationPts

        close(1)
        close(2)
        close(3)
        if(isContainer)    close(4)

    end subroutine setMeshData

    subroutine setArea(this)
    !! sets the area of the rigid body. 
        class(container_class) :: this
        integer :: cp
        
        this%area = 0.0_rp
        
        do cp = 1, this%nCollocationPts
            this%area = this%area + this%wgts(cp)
        end do
    end subroutine setArea

    subroutine setXC(this)
    !! Make sure to call setArea before setting xC
        class(container_class) :: this
        real(kind=rp), dimension(3) :: x
        integer :: cp
        
        this%xC = (/0.0_rp, 0.0_rp, 0.0_rp/)

        do cp = 1, this%nCollocationPts
            x = this%collocationPts(cp, :)
            this%xC = this%xC + x*this%wgts(cp)
        end do

        this%xC = this%xC/this%area

    end subroutine setXC

    subroutine setItensor(this)
    !! sets inertia tensor about the geometric center    
        class(container_class) :: this
        real(kind=rp), dimension(3) :: delX
        integer :: cp
        
        this%Itensor = 0.0_rp
        
        do cp = 1, this%nCollocationPts
            delX = this%collocationPts(cp, :) - this%xC
            this%Itensor(1, 1:3) = this%Itensor(1, 1:3) + this%wgts(cp)*dot_product(delX, delX)*(/1.0, 0.0, 0.0/) &
                                    - delX*delX(1)*this%wgts(cp)
            this%Itensor(2, 2:3) = this%Itensor(2, 2:3) + this%wgts(cp)*dot_product(delX, delX)*(/1.0, 0.0/) &
                                    - delX(2:3)*delX(2)*this%wgts(cp)
            this%Itensor(3, 3) = this%Itensor(3, 3) + this%wgts(cp)*dot_product(delX, delX) &
                                    - delX(3)*delX(3)*this%wgts(cp)
        end do

        ! use symmetry to assign other elements:
        this%Itensor(2, 1) = this%Itensor(1, 2)
        this%Itensor(3, 1) = this%Itensor(1, 3)
        this%Itensor(3, 2) = this%Itensor(2, 3)

        call inv3dMat(this%Itensor, this%ItensorInv)

    end subroutine setItensor

    subroutine writeMesh(this, fileName)
        class(container_class) :: this
        integer :: iter, io
        character(*) :: fileName
        
        open(unit = 1, file = fileName, iostat = io)
        if(io /= 0) print*, 'Cant open mesh file to write'
        
        do iter = 1, this%nCollocationPts
            write(1, *) this%collocationPts(iter,:)
        end do
        close(1)

    end subroutine writeMesh

    function VIOmega(this, q0)
    !! returns array of V*area and Itensor.Omega obtained by applying rigid body projection to qDL of the body.
        class(container_class) :: this
        real(kind=rp), dimension(3) :: delX
        real(kind=rp), dimension(6) :: VOmegaPartial, VOmegaFull
        real(kind=rp), dimension(3, 2) :: VIOmega
        integer :: cp
        real(kind=rp), dimension(:, :) :: q0 

        ! communication latency may win over faster computation. Test it without mpi
        ! VPartial = 0.0_rp
        VOmegaPartial = 0.0_rp
        do cp = this%myStartIndx, this%myEndIndx
        ! do cp = 1, this%nCollocationPts
            delX = this%collocationPts(cp, :) - this%xC
            VOmegaPartial(1:3) = VOmegaPartial(1:3) + this%wgts(cp)*q0(cp, :)
            VOmegaPartial(4:6) = VOmegaPartial(4:6) + cross(delX, q0(cp, :))*this%wgts(cp)
        end do
        
        call MPI_Allreduce(VOmegaPartial, VOmegaFull, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierror)
        
        VIOmega(:, 1) = VOmegaFull(1:3)
        VIOmega(:, 2) = VOmegaFull(4:6)
        
    end function VIOmega

    function qDotNhat(this, q0)
    !! returns Pn (normal projection) acting on qDL of the body
        class(container_class) :: this
        real(kind=rp) :: qDotNhat
        integer :: cp
        real(kind=rp), dimension(:, :) :: q0 

        ! calculate <n, q> using mpi:
        ! communication latency may win over faster computation. Test it without mpi
        qDotNhat = 0.0_rp
        ! do cp = this%myStartIndx, this%myEndIndx
        do cp = 1, this%nCollocationPts
            qDotNhat = qDotNhat + this%wgts(cp)*dot_product(q0(cp, :), this%normalVectors(cp, :))
        end do
        
        ! call MPI_Allreduce(qDotnPartial, qDotn, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierror)
    
    end function qDotNhat

    function Pn(this, q0)
    !! returns Pn (normal projection) acting on qDL of the body
        class(container_class) :: this
        real(kind=rp), dimension(this%nCollocationPts, 3) :: Pn
        integer :: cp
        real(kind=rp) :: qDotn!, qDotnPartial
        real(kind=rp), dimension(:, :) :: q0 

        qDotn = this%qDotNhat(q0)
        
        ! once <n, q> is calculated, get the projections at collocation points.
        do cp = 1, this%nCollocationPts
            Pn(cp, :) = this%normalVectors(cp,:)*qDotn/this%area
        end do

    end function Pn

    function LdX0(this, q0, x0, cp0, self)
    !! returns the double layer operator Ld acting on qDL of the body at x0
        class(container_class) :: this
        real(kind=rp), dimension(3) :: nHat, dx
        real(kind=rp), dimension(3) :: qdotTdotnHat, tmp, x0, LdX0
        integer :: cp, cp0
        logical :: self
        real(kind=rp), dimension(:, :) :: q0 

        ! this is order(nCollocationPts^2) calculations for each x0 in collocationPts... very expensive. If only the Kernel was separable! (See FMM)
        tmp = 0.0_rp
        do cp = this%myStartIndx, this%myEndIndx
            nHat = this%normalVectors(cp, :)
            dx = this%collocationPts(cp, :) - x0    
            if(cp == cp0 .and. self .eqv. .true.  )   then
                cycle
            endif

            ! if(cp0==-1)   print*, 'yo', q0(cp0, 1)

            if ( self .eqv. .true. ) then
                ! doing singularity subtraction:
                qdotTdotnHat = -6.0_rp*dx*dot_product( q0(cp, :) &
                                                        - q0(cp0, :), dx)*dot_product(nHat, dx)/(norm2(dx)**5)
            else
                qdotTdotnHat = -6.0_rp*dx*dot_product( q0(cp, :) , dx)*dot_product(nHat, dx)/(norm2(dx)**5)
            end if
            
            tmp = tmp + this%wgts(cp)*qdotTdotnHat/(4.0_rp*pi)
        end do

        call MPI_Allreduce(tmp, LdX0, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierror)
        
    end function LdX0

end module container_mod