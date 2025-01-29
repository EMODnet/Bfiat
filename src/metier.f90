!===============================================================================
!===============================================================================
! Estimates densities at selected time for the perturbation model 
! with several metiers (analytic solution)
!===============================================================================
!===============================================================================

!===============================================================================
! Input: sar        : per species and metier, 
!        depletion d: per species and metier
!===============================================================================

       SUBROUTINE metier_time   (nspec, nmetier, ntimes, B0,                    & 
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
       
       END SUBROUTINE metier_time


!===============================================================================
! analytical solution of the metier model with events
! Same sar for all species -> input of event times
!===============================================================================

       SUBROUTINE metier_times_events(nspec, nmetier, nevent, ntimes, B0,       &
                                K, r, d, times, events, metier,                 &
                                B, dTrawl)
          
       IMPLICIT NONE
       
       
       ! Number of taxa, number of metiers, length of events
       INTEGER, INTENT(IN) :: nspec, nmetier,  nevent, ntimes
       
       ! taxon-specific parameters
       DOUBLE PRECISION, INTENT(IN) :: K(nspec)    ! carrying capacity
       DOUBLE PRECISION, INTENT(IN) :: r(nspec)    ! rate of increase

       ! taxon x metier specific
       DOUBLE PRECISION, INTENT(IN) :: d(nspec, nmetier)    ! depletion fraction

        ! output times and trawling
       DOUBLE PRECISION, INTENT(IN) :: times(ntimes)  ! output times
       DOUBLE PRECISION, INTENT(IN) :: events(nevent) ! timing of events
       INTEGER, INTENT(IN)          :: metier(nevent) ! metier of this event

       ! initial condition
       DOUBLE PRECISION, INTENT(INOUT) :: B0(nspec)

       ! output varialbes
       DOUBLE PRECISION, INTENT(OUT)   :: B(nspec, ntimes)
       
       ! all the trawling-induced mortalities
       DOUBLE PRECISION, INTENT(OUT)   :: dTrawl(nspec, nevent)

       INTEGER :: I, J, N, nT, nE, iE
       DOUBLE PRECISION :: tnext, tt, t0, denom
    
        nT     = 0              ! index to start of times in loop
        t0     = times(1)       ! time in loop

        nE     = 1              ! index to first event  
        
        IF (events(1) == t0) THEN  ! start with an event => update initial cond.
          nE   = 2
          
          DO N = 1, nspec
            iE = metier(1)
            B0(N)        = B0(N)*(1.d0 -d(N, iE))   ! depletion on initial condition
            dTrawl(N, 1) = B0(N)*d(N, iE)
          END DO
        
        END IF 
               
        ! Go from event to event - inbetween the solution can be estimated analytically         
        DO I = nE, nevent
          
          IF (I < nevent) THEN          ! there is a next event
            tnext = events(I)
          ELSE                          ! no next event - set it = last time+1
            tnext = times(ntimes)+1.d0  
          END IF
          
          DO J = nT+1, ntimes           ! loop over all times
            
            DO N = 1, nspec             ! for all species
            ! gives NaN for B and K=0
              denom  = B0(N)+(K(N)-B0(N))*Dexp(-r(N)*(times(j)-t0))
              B(N,J) = B0(N)*K(N)/Denom
            
            END DO
            IF (times(j) > tnext) EXIT  ! first take care of event
          
          END DO   ! with J
          
          nT     = J - 1

          ! new initial condition after trawling 
          DO N = 1, nspec
              iE = metier(I)
              Denom  = B0(N)+(K(N)-B0(N))*Dexp(-r(N)*(tnext-t0))
              B0(N)  = B0(N)*K(N)/Denom
              IF (times(j) > tnext) dTrawl(N, I) = B0(N)*d(N, iE)
          
              B0(N)  = B0(N)*(1.d0 -d(N, iE))
          END DO
          t0     = tnext
        
        END DO  ! with I
       END SUBROUTINE metier_times_events
       
