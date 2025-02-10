!===============================================================================
!===============================================================================
! logistic model with constant mortality
!===============================================================================
!===============================================================================

!===============================================================================
! Calculates the analytical solution of the logistic model with ct mortality
! for specific times
!===============================================================================

       SUBROUTINE logistic_time(nspec, ntimes, B0, K, r, m,                     & 
                           times, tstartperturb, tendperturb, B)
          
       IMPLICIT NONE
       
       INTEGER, INTENT(IN) :: nspec, ntimes
       
       DOUBLE PRECISION, INTENT(IN) :: K(nspec)
       DOUBLE PRECISION, INTENT(IN) :: r(nspec), m(nspec)
       DOUBLE PRECISION, INTENT(IN) :: times(ntimes)
       DOUBLE PRECISION, INTENT(IN) :: tstartperturb, tendperturb

       DOUBLE PRECISION, INTENT(INOUT) :: B0(nspec)

       DOUBLE PRECISION, INTENT(OUT) :: B(nspec, ntimes)

       INTEGER :: N, J
       DOUBLE PRECISION :: t0, D, rn
       LOGICAL :: pastEnd
       LOGICAL :: pastStart
    
         t0  = times(1)
         DO N = 1, nspec     ! for all species
           
           D  = B0(N)          ! initial value
           
           pastEnd   = .FALSE.
           pastStart = .FALSE.
           
           IF (tstartperturb <= times(1)) pastStart = .TRUE.
           
           DO J = 1, ntimes
           
              ! pass the start of fishing         
              IF (times(J) > tstartperturb .AND. .NOT. pastStart ) THEN
                pastStart = .TRUE.
                
                ! new initial condition
                D  = rn*K(N)*D/(r(N)*D +                                        &
                     (rn*K(N)-r(N)*D)*Dexp(-rn*(tstartperturb - t0)))
                t0 = tstartperturb
              ENDIF
              
              ! pass the end of fishing         
              IF (times(J) > tendperturb .AND. .NOT. pastEnd ) THEN
                pastEnd = .TRUE.
                D  = rn*K(N)*D/(r(N)*D +                                        &
                     (rn*K(N)-r(N)*D)*Dexp(-rn*(tendperturb - t0)))
                t0 = tendperturb
              ENDIF
              
              IF (times(J) <= tstartperturb .OR. times(J) > tendperturb) THEN
                 rn = r(N)
              ELSE 
                 rn = r(N) - m(N)    ! net rate of increase
              END IF

              B(N,J) = rn*K(N)*D/(r(N)*D +                                      &
                       (rn*K(N)-r(N)*D)*Dexp(-rn*(times(J) - t0)))
    
            END DO  ! with J    
         END DO     ! with N

       END SUBROUTINE logistic_time
       
       