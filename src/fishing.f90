!===============================================================================
! Calculates the analytical solution of the logistic model with events
!===============================================================================

       SUBROUTINE logistictrawl(nspec, nevent, ntimes, B0,                     & 
                                K, r, d, times, events, B, dTrawl)
          
       IMPLICIT NONE
       
       INTEGER, INTENT(IN) :: nspec, nevent, ntimes
       DOUBLE PRECISION, INTENT(IN) :: K(nspec)
       DOUBLE PRECISION, INTENT(IN) :: r(nspec), d(nspec)
       DOUBLE PRECISION, INTENT(IN) :: times(ntimes)
       DOUBLE PRECISION, INTENT(IN) :: events(nevent)

       DOUBLE PRECISION, INTENT(INOUT) :: B0(nspec)

       DOUBLE PRECISION, INTENT(OUT) :: B(nspec, ntimes)
       DOUBLE PRECISION, INTENT(OUT) :: dTrawl(nspec, nevent)

       INTEGER :: I, J, N, nT, nE
       DOUBLE PRECISION :: tevent, tnext, tt, t0, denom
    
        nT     = 0             ! index to start of times in loop
        t0     = times(1)      ! time in loop
        tevent = events(1)     ! first event   

        nE     = 1             ! index to first event  
        IF (tevent == t0) THEN  ! start with an event
          nE   = 2
          DO N = 1, nspec
            B0(N)        = B0(N)*(1.d0 -d(N))
            dTrawl(N, 1) = B0(N)*d(N)
          END DO
        END IF 
                
        DO I = nE, nevent
          IF (I < nevent) THEN
            tnext = events(I)
          ELSE
            tnext = times(ntimes)+1.d0
          END IF
          
          DO J = nT+1, ntimes
            IF (times(j) < tnext) THEN
              DO N = 1, nspec     ! for all species
            ! gives NaN for B and K=0
                Denom  = B0(N)+(K(N)-B0(N))*Dexp(-r(N)*(times(j)-t0))
                B(N,J) = B0(N)*K(N)/Denom
              END DO
            ELSE 
              EXIT
            END IF
          END DO   ! with J
          
          nT   = J-1
          t0     = tnext
          tevent = tnext
          ! new initial condition after trawling 
          DO N = 1, nspec
            B0(N)        = B(N, nT)*(1.d0 -d(N))
            dTrawl(N, I) = B(N, nT)*d(N)
          END DO
        END DO  
       END SUBROUTINE logistictrawl
       
!===============================================================================
! Calculates the analytical solution of the logistic model with ct mortality
!===============================================================================

       SUBROUTINE logistic(nspec, ntimes, B0,                     & 
                           K, r, m, times, B, mort)
          
       IMPLICIT NONE
       
       INTEGER, INTENT(IN) :: nspec, ntimes
       DOUBLE PRECISION, INTENT(IN) :: K(nspec)
       DOUBLE PRECISION, INTENT(IN) :: r(nspec), m(nspec)
       DOUBLE PRECISION, INTENT(IN) :: times(ntimes)

       DOUBLE PRECISION, INTENT(INOUT) :: B0(nspec)

       DOUBLE PRECISION, INTENT(OUT) :: B(nspec, ntimes)
       DOUBLE PRECISION, INTENT(OUT) :: mort(nspec, ntimes)

       INTEGER :: N
       DOUBLE PRECISION :: t(ntimes), D, rn
    
         t    = times-times(1)
         DO N = 1, nspec     ! for all species
           rn = r(N)-m(N)    ! net rate of increase
           D  = B0(N)
           B(N,:) = rn*K(N)*D/(r(N)*D + (rn*K(N)-r(N)*D)*Dexp(-rn*(t)))
           mort(N,:) = m(N) * B(N,:)
         END DO   ! with n

       END SUBROUTINE logistic
       