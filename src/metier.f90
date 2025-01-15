!===============================================================================
!===============================================================================
! Runs the perturbation model for several metiers (analytic solution)
!===============================================================================
!===============================================================================

!===============================================================================
! Input: sar        : per species and metier, 
!        depletion d: per species and metier
!===============================================================================

       SUBROUTINE logisticmetier(nspec, nmetier, ntimes, B0,                    & 
                                 K, r, d, sar, times, tendperturb, B)
          
       IMPLICIT NONE
       
       
       ! Number of taxa, number of metiers, length of times
       INTEGER, INTENT(IN) :: nspec, nmetier,  ntimes
       
       ! taxon-specific parameters
       DOUBLE PRECISION, INTENT(IN) :: K(nspec)    ! carrying capacity
       DOUBLE PRECISION, INTENT(IN) :: r(nspec)    ! rate of increase

       ! metier-specific parameters
       DOUBLE PRECISION, INTENT(IN) :: tendperturb(nmetier) ! last trawling
       
       ! taxon x metier specific
       DOUBLE PRECISION, INTENT(IN) :: d(nspec, nmetier)    ! depletion fraction
       DOUBLE PRECISION, INTENT(IN) :: sar(nspec, nmetier)  ! swept area ratio
      
       ! times
       DOUBLE PRECISION, INTENT(IN) :: times(ntimes)  ! output times

       ! output
       DOUBLE PRECISION, INTENT(INOUT) :: B0(nspec)
       DOUBLE PRECISION, INTENT(OUT)   :: B(nspec, ntimes)
       
       ! local variables
       INTEGER :: I, J, N
       DOUBLE PRECISION :: dtevent(nmetier)     ! time between trawling for this metier
       DOUBLE PRECISION :: tnextmetier(nmetier) ! next trawling time for this metier
       DOUBLE PRECISION :: tprev                ! time of last trawl
       DOUBLE PRECISION :: dt                   ! time between trawls
                           
       DOUBLE PRECISION :: denom, tnext
       
        DO N = 1, nspec
        
        ! initialisation
        ! ----------------
        
         tprev = times(1)
         tnext = times(ntimes) +1.d0
         
         DO I = 1, nmetier
           dtevent(I)       = 1d0/sar(N, I)
           tnextmetier(I)   = times(1) + 0.5D0 * dtevent(I) 
           IF (tnextmetier(I) < tnext) tnext = tnextmetier(I)
         END DO
         
         ! save initial condition
         B(N, 1) = B0(N)
         
         ! loop over all times         
         DO J = 2, ntimes
           
           ! take care of the event(s) first
           
           DO WHILE (times(J) > tnext)    
              
              dt = tnext - tprev

             ! new initial condition before the event(s)
              Denom  = B0(N) + (K(N)-B0(N))*Dexp(-r(N)*dt)
              B0(N)  = B0(N) * K(N)/Denom   
             
             ! which metier is trawling (can be more than one)
              DO I = 1, nmetier
                
                IF (tnextmetier(I) == tnext) THEN
                
                  ! new condition after the event
                  B0(N)  = B0(N)*(1.d0 -d(N, I)) 
                  
                  ! update the metier
                  tnextmetier(I)  = tnextmetier(I) + dtevent(I) 
                  
                  ! check if still trawling
                  IF (tnextmetier(I) > tendperturb(I)) THEN
                     tnextmetier(I) = times(ntimes) + 1.d0 
                  END IF
                  
                END IF  ! tnextmetier == tnext
              
              END DO    ! I
              
              tprev = tnext
              
              ! next trawl
              tnext = times(ntimes) +1.d0
         
              DO I = 1, nmetier
                IF (tnextmetier(I) < tnext) THEN
                  tnext = tnextmetier(I)
                END IF  
              END DO
           
           END DO   ! with trawling

           ! Estimate density at times(j)
            Denom  = B0(N) + (K(N)-B0(N))*Dexp(-r(N)*(times(j)-tprev))
            B(N,J) = B0(N)*K(N)/Denom

          END DO   ! with J
        END DO  ! with N
       
       END SUBROUTINE logisticmetier

