module matrix_mod
    use global_mod
    implicit none
    
contains
    
    function cross(v1, v2)
    !! compute the cross product of 3d vectors v1 and v2
        real(kind=rp), dimension(3) :: v1, v2, cross
        cross = (/v1(2)*v2(3) - v1(3)*v2(2), &
                    v1(3)*v2(1) - v1(1)*v2(3), v1(1)*v2(2) - v1(2)*v2(1)  /)         
    end function cross

    subroutine inv3dMat(matA, matAInv)
    !! compute the inverse of a 3x3 matrix matA. Returns |A|AInv for singular matrix:
        real(kind=rp), dimension(3, 3) :: matA, matAInv
        real(kind=rp), dimension(3) :: vec1, vec2, vec3
        real(kind=rp) :: det
        ! set ItensorInv now:
        vec1 = matA(1, :)
        vec2 = matA(2, :)
        vec3 = matA(3, :)
        
        det = dot_product(vec1, cross(vec2, vec3) )

        if (det == 0.0_rp .and. myRank == 0) print*, 'Singular matrix!'

        matAInv(1, :) = cross(vec2, vec3)
        matAInv(2, :) = cross(vec3, vec1)
        matAInv(3, :) = cross(vec1, vec2)
        
        if (det /= 0.0_rp) matAInv = transpose(matAInv)/det

    end subroutine inv3dMat

    function rotate3dVec(xI, w, theta) result(xF)
    !! rotate the 3d vector xI along w axis by |w|*theta. w need not be unit vector.
        real(kind=rp), dimension(3) :: w, xI, xF
        real(kind=rp) :: theta, w1Sq, w2Sq, w3Sq, wSq, phi, modW
        
        modW = norm2(w)
        wSq = modW**2
        w1Sq = w(1)**2
        w2Sq = w(2)**2
        w3Sq = w(3)**2
        phi = modW*theta
        
        if(phi /= 0.0_rp) then
            ! solve dr/dt = wxr to get these:
            xF(1) = ( xI(1)*( w1Sq + (w2Sq+w3Sq)*cos(phi) ) + xI(3)*( w(1)*w(3)*(1.0_rp-cos(phi)) &
                                + w(2)*modW*sin(phi) ) + xI(2)*( w(1)*w(2)*(1.0_rp-cos(phi)) - w(3)*modW*sin(phi) ) )/wSq;

            xF(2) = ( xI(2)*( w2Sq + (w1Sq+w3Sq)*cos(phi) ) + xI(3)*( w(2)*w(3)*(1.0_rp-cos(phi)) &
                                - w(1)*modW*sin(phi) ) + xI(1)*( w(1)*w(2)*(1.0_rp-cos(phi)) + w(3)*modW*sin(phi) ) )/wSq;

            xF(3) = ( xI(3)*( w3Sq + (w1Sq+w2Sq)*cos(phi) ) + xI(2)*( w(2)*w(3)*(1.0_rp-cos(phi)) &
                                + w(1)*modW*sin(phi) ) + xI(1)*( w(1)*w(3)*(1.0_rp-cos(phi)) - w(2)*modW*sin(phi) ) )/wSq;
        
        else
            xF = xI
        endif

    end function rotate3dVec

    function mat3dVec(Amat, vec)
        real(kind=rp), dimension(3), intent(in) :: vec
        real(kind=rp), dimension(3) :: mat3dVec
        real(kind=rp), dimension(3,3) :: Amat

        mat3dVec(1) = dot_product(Amat(1, :), vec)
        mat3dVec(2) = dot_product(Amat(2, :), vec)
        mat3dVec(3) = dot_product(Amat(3, :), vec)
        
    end function

    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition. Good for large Matrix A
    ! function inv(A) result(Ainv)
    !     real(kind=rp), dimension(:,:), intent(in) :: A
    !     real(kind=rp), dimension(size(A,1),size(A,2)) :: Ainv
    
    !     real(kind=rp), dimension(size(A,1)) :: work  ! work array for LAPACK
    !     integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    !     integer :: n, info
    
    !     ! External procedures defined in LAPACK
    !     external DGETRF
    !     external DGETRI
    
    !     ! Store A in Ainv to prevent it from being overwritten by LAPACK
    !     Ainv = A
    !     n = size(A,1)
    
    !     ! DGETRF computes an LU factorization of a general M-by-N matrix A
    !     ! using partial pivoting with row interchanges.
    !     call DGETRF(n, n, Ainv, n, ipiv, info)
    
    !     if (info /= 0) then
    !     stop 'Matrix is numerically singular!'
    !     end if
    
    !     ! DGETRI computes the inverse of a matrix using the LU factorization
    !     ! computed by DGETRF.
    !     call DGETRI(n, Ainv, n, ipiv, work, n, info)
    
    !     if (info /= 0) then
    !     stop 'Matrix inversion failed!'
    !     end if
    ! end function inv

end module matrix_mod

