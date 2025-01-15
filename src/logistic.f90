!===============================================================================
!===============================================================================
! logistic model with constant mortality
!===============================================================================
!===============================================================================

!===============================================================================
! Calculates the analytical solution of the logistic model with ct mortality
!===============================================================================

       SUBROUTINE logistic(nspec, ntimes, B0,                                   & 
                           K, r, m, times, tendperturb, B)
          
       IMPLICIT NONE
       
       INTEGER, INTENT(IN) :: nspec, ntimes
       
       DOUBLE PRECISION, INTENT(IN) :: K(nspec)
       DOUBLE PRECISION, INTENT(IN) :: r(nspec), m(nspec)
       DOUBLE PRECISION, INTENT(IN) :: times(ntimes), tendperturb

       DOUBLE PRECISION, INTENT(INOUT) :: B0(nspec)

       DOUBLE PRECISION, INTENT(OUT) :: B(nspec, ntimes)

       INTEGER :: N, J
       DOUBLE PRECISION :: t0, D, rn
       LOGICAL pastEnd
    
         t0  = times(1)
         DO N = 1, nspec     ! for all species
           
           rn = r(N) - m(N)    ! net rate of increase
           D  = B0(N)          ! initial value
           pastEnd = .FALSE.
           
           DO J = 1, ntimes
              B(N,J) = rn*K(N)*D/(r(N)*D +                                      &
                       (rn*K(N)-r(N)*D)*Dexp(-rn*(times(J) - t0)))
              IF (times(J) > tendperturb .AND. .NOT. pastEnd ) THEN
                pastEnd = .TRUE.
                D  = rn*K(N)*D/(r(N)*D +                                        &
                     (rn*K(N)-r(N)*D)*Dexp(-rn*(tendperturb - t0)))
                rn = r(N)  ! mortality is 0
                t0 = tendperturb
              ENDIF
            END DO  ! with J    
         END DO     ! with N

       END SUBROUTINE logistic
       
       