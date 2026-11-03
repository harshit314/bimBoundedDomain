module rigidBody_mod
    ! use mpi_f08 ! recommended over the mpif.h
    use container_mod
    implicit none
    ! include 'mpif.h'
    !private -> all the variables and procedures in surfaceMesh_mod will become private
    ! define rigid body class:
    type, extends(container_class) :: rigidBody_class

        real(kind=rp), dimension(3) :: pz, px, py, vel = (/0.0_rp, 0.0_rp, 0.0_rp/), omega = (/0.0_rp, 0.0_rp, 0.0_rp/)
        real(kind=rp), dimension(3) :: Fh = (-3.0_rp/4.0_rp)*(/0.0_rp, -1.0_rp, 0.0_rp/), Th = (/0.0_rp, 0.0_rp, 0.0_rp/)
        real(kind=rp), dimension(:, :), allocatable :: xFT
        integer :: nforcePt ! no. of point forces and torques needed for supplementary flow

    contains

        procedure :: setSupplementaryFlow
        procedure :: setOrientation

        procedure :: translateTo
        procedure :: translate
        procedure :: rotate
        
        ! BIM related functions:
        procedure :: VOmega
        procedure :: Prb

        procedure :: updateVel

    end type rigidBody_class

contains

    subroutine setSupplementaryFlow(this, posPts)
    !! places the point forces and torques at posPts, dividing force/torque magnitudes among posPts equally. 
    !! If needed, change forces and torques from default values before calling.
        class(rigidBody_class) :: this
        real(kind=rp), dimension(:,:) :: posPts

        this%nforcePt = size(posPts, 1)
        if(.not. allocated(this%xFT))   allocate(this%xFT(this%nforcePt, 3))
        
        this%xFT = posPts
        
    end subroutine setSupplementaryFlow

   subroutine setOrientation(this, px, py, pz)
    !! set the body fixed coordinate system as per the imported mesh.
        class(rigidBody_class) :: this
        real(kind=rp), dimension(3) :: px, py, pz

        this%px = px
        this%py = py
        this%pz = pz
        
    end subroutine setOrientation

    ! functions that transform the rigid body. 
    ! Must update normalVectors and collocationPts

    subroutine translate(this, dx)
    !! translate the rigid body to xTo position. Normal vectors dont have to be updated
        class(rigidBody_class) :: this
        real(kind=rp), dimension(3) :: dx
        integer :: cp, iFT
        
        ! shift the global coordinates:
        do cp = 1, this%nCollocationPts
            this%collocationPts(cp, :) = this%collocationPts(cp, :) + dx
        end do
        
        this%xC = this%xC + dx
         
        do iFT = 1, this%nforcePt
            this%xFT(iFT, :) = this%xFT(iFT, :) + dx
        end do
        
    end subroutine translate
    
    subroutine translateTo(this, xTo)
    !! translate the rigid body to xTo position. Normal vectors dont have to be updated
        class(rigidBody_class) :: this
        real(kind=rp), dimension(3) :: xTo, dx
        
        dx = xTo - this%xC
        call this%translate(dx)

    end subroutine translateTo

    subroutine rotate(this, w, theta)
    !! rotate the rigid body along w axis by |w|*theta about its center xC. w need not be unit vector.
        class(rigidBody_class) :: this
        real(kind=rp), dimension(3) :: w, xI, xF
        real(kind=rp) :: theta
        integer :: cp, iFT

        ! shift the collocationPts and normalVectors:
        do cp = 1, this%nCollocationPts
            ! bring origin to body's center:
            xI = this%collocationPts(cp, :) - this%xC
            ! now rotate:
            xF = rotate3dVec(xI, w, theta)
            this%collocationPts(cp, :) = this%xC + xF 
            ! rotate normal vectors:
            this%normalVectors(cp, :) = rotate3dVec(this%normalVectors(cp,:), w, theta)
        end do
            
        call this%setItensor()
        ! update orientation vector:
        this%px = rotate3dVec(this%px, w, theta) 
        this%py = rotate3dVec(this%py, w, theta) 
        this%pz = rotate3dVec(this%pz, w, theta) 
        
        do iFT = 1, this%nforcePt
            ! bring origin to body's center:
            xI = this%xFT(iFT, :) - this%xC
            this%xFT(iFT, :) = rotate3dVec(xI, w, theta) + this%xC
        end do

    end subroutine rotate

    function VOmega(this, q0)
    !! returns array of V and Omega obtained by applying rigid body projection to qDL of the body
        class(rigidBody_class) :: this
        real(kind=rp), dimension(3) :: delX
        real(kind=rp), dimension(6) :: VOmegaPartial, VOmegaFull
        real(kind=rp), dimension(3, 2) :: VOmega
        integer :: cp
        real(kind=rp), dimension(:,:) :: q0

        ! communication latency may win over faster computation. Test it without mpi
        ! VPartial = 0.0_rp
        VOmegaPartial = 0.0_rp
        VOmegaFull = 0.0_rp
        do cp = this%myStartIndx, this%myEndIndx
        ! do cp = 1, this%nCollocationPts
            delX = this%collocationPts(cp, :) - this%xC
            VOmegaPartial(1:3) = VOmegaPartial(1:3) + this%wgts(cp)*q0(cp, :)/this%area
            VOmegaPartial(4:6) = VOmegaPartial(4:6) + cross(delX, q0(cp, :))*this%wgts(cp)
        end do
        
        call MPI_Allreduce(VOmegaPartial, VOmegaFull, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierror)
        
        VOmega(:, 1) = VOmegaFull(1:3)
        VOmega(:, 2) = mat3dVec(this%ItensorInv, VOmegaFull(4:6))
        
    end function VOmega

    function Prb(this, q0)
    !! returns Prb acting on qDL of the body
        class(rigidBody_class) :: this
        real(kind=rp), dimension(3) :: delX
        real(kind=rp), dimension(this%nCollocationPts, 3) :: Prb
        real(kind=rp), dimension(:, :) :: q0 

        real(kind=rp), dimension(3, 2) :: Vgen
        integer :: cp
        
        Vgen = this%VOmega(q0)
        do cp = 1, this%nCollocationPts
            delX = this%collocationPts(cp, :) - this%xC
            Prb(cp, :) = Vgen(:, 1) + cross( Vgen(:, 2), delX )
        end do

    end function Prb

    subroutine updateVel(this)
    !! Set the velocity of the body using the double layer density
        class(rigidBody_class) :: this
        real(kind=rp), dimension(3, 2) :: vTmp
        
        vTmp = this%VOmega(this%qDLFinal)

        this%vel = vTmp(:, 1)
        this%omega = vTmp(:, 2)

    end subroutine updateVel

end module rigidBody_mod