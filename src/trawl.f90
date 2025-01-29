!===============================================================================
!===============================================================================
! Density of the perturbation model at selected events for several species
! One metier
!===============================================================================
!===============================================================================

       SUBROUTINE perturb_event(nspec, nevent, B0, sar,                         & 
                               K, r, d, eventnr, B, Bend, Bmean)
          
       IMPLICIT NONE
       
       
       ! Number of taxa in parameter data, number of events
       INTEGER, INTENT(IN) :: nspec, nevent
       
       DOUBLE PRECISION, INTENT(IN) :: K(nspec)   ! carrying capacity
       DOUBLE PRECISION, INTENT(IN) :: r(nspec)   ! rate of increase
       DOUBLE PRECISION, INTENT(IN) :: d(nspec)   ! depletion fraction
       DOUBLE PRECISION, INTENT(IN) :: B0(nspec)  ! initial condition
       DOUBLE PRECISION, INTENT(IN) :: sar(nspec) ! swept area ratio
       
       INTEGER, INTENT(IN) :: eventnr(nevent)    ! events to keep (sorted)

       ! output: density for each speces at requested events
       DOUBLE PRECISION, INTENT(OUT) :: B(nspec, nevent) ! 
       DOUBLE PRECISION, INTENT(OUT) :: Bend(nspec, nevent) ! 
       DOUBLE PRECISION, INTENT(OUT) :: Bmean(nspec, nevent)
       
       INTEGER :: I, N, maxevents, nev
       DOUBLE PRECISION :: p, dens, expon, denom
    
        maxevents = eventnr(nevent)  ! last event is max event nr (i.e. sorted)
        
        ! for every species
        
        DO N = 1, nspec
          
          expon  = Dexp(-r(N)/sar(N))
          p      = 1.d0-d(N)          ! proportion remaining
          dens   = B0(N)              ! first fishing
          nev    = 1
         
         ! Go from event to event         
          DO I = 1, maxevents
              
              ! save only if this event is requested
              DO WHILE (I >= eventnr(nev))   
              
               ! the density AFTER fishing
                B(N, nev) = p*dens  

                Denom     = p*dens + (K(N) - p*dens) * expon
                
               ! the mean density in the fishing interval
                Bmean(N, nev) = K(N) + sar(N)*K(N)/r(N)*Dlog(Denom/K(N))
               
               ! the density at the end of the fishing interval
                Bend(N, nev) = p*dens*K(N) / Denom
                
                nev = nev + 1
                IF (nev > nevent) EXIT
              END DO 
              
              IF (nev > nevent) EXIT

              ! New density
              Denom  = p*dens + (K(N) - p*dens) * expon
              dens   = p*dens*K(N) / Denom
          END DO   ! with I
          
         END DO    ! with N

       END SUBROUTINE perturb_event

!===============================================================================
! Density at events, one event per species
!===============================================================================

       SUBROUTINE perturb_event2(nspec, B0, sar,                                & 
                               K, r, d, eventnr, B)
          
       IMPLICIT NONE
       
       
       ! Number of taxa in parameter data - one event per species
       INTEGER, INTENT(IN) :: nspec
       
       DOUBLE PRECISION, INTENT(IN) :: K(nspec)   ! carrying capacity
       DOUBLE PRECISION, INTENT(IN) :: r(nspec)   ! rate of increase
       DOUBLE PRECISION, INTENT(IN) :: d(nspec)   ! depletion fraction
       DOUBLE PRECISION, INTENT(IN) :: B0(nspec)  ! initial condition
       DOUBLE PRECISION, INTENT(IN) :: sar(nspec) ! swept area ratio
       
       INTEGER, INTENT(IN) :: eventnr(nspec)      ! events at which to estimate

       ! output: density for each speces at requested events
       DOUBLE PRECISION, INTENT(OUT) :: B(nspec) ! 
       
       INTEGER :: I, N, nev
       DOUBLE PRECISION :: p, dens, expon, denom
    
        ! for every species
        
        DO N = 1, nspec
          
          expon  = Dexp(-r(N)/sar(N))
          p      = 1.d0-d(N)          ! proportion remaining
          dens   = B0(N)              ! first fishing
          nev    = eventnr(N)         ! nr of events to estimate
         
         ! Go from event to event         
          DO I = 1, nev
              ! New density
              Denom  = p*dens + (K(N) - p*dens) * expon
              dens   = p*dens*K(N) / Denom
          END DO   ! with I
          B(N) = p*dens  ! save the density AFTER fishing

         END DO    ! with N

       END SUBROUTINE perturb_event2


!===============================================================================
!===============================================================================
! Density of the perturbation model at selected events for several species
! Several metiers
!===============================================================================
!===============================================================================

       SUBROUTINE metier_event(nspec, nmetier, nevent, B0, sar,                 & 
                     K, r, d, eventnr, B, Bend, Bmean, times, tend)
          
       IMPLICIT NONE
       
       
       ! Number of taxa in parameter data, number of events
       INTEGER, INTENT(IN) :: nspec, nmetier, nevent
       
       DOUBLE PRECISION, INTENT(IN) :: K(nspec)            ! carrying capacity
       DOUBLE PRECISION, INTENT(IN) :: r(nspec)            ! rate of increase
       DOUBLE PRECISION, INTENT(IN) :: d(nspec, nmetier)   ! depletion fraction
       DOUBLE PRECISION, INTENT(IN) :: B0(nspec)           ! initial condition
       DOUBLE PRECISION, INTENT(IN) :: sar(nspec, nmetier) ! swept area ratio
       
       INTEGER, INTENT(IN) :: eventnr(nevent)    ! events to keep (sorted)

       ! output: density for each speces at requested events
       DOUBLE PRECISION, INTENT(OUT) :: B(nspec, nevent), Bend(nspec, nevent) ! 
       DOUBLE PRECISION, INTENT(OUT) :: Bmean(nspec, nevent)
       DOUBLE PRECISION, INTENT(OUT) :: times(nspec, nevent)
       DOUBLE PRECISION, INTENT(OUT) :: tend(nspec, nevent)   
       
       ! local variables
       INTEGER :: I, J, N, maxevents, nev, inetxt
       DOUBLE PRECISION :: p, dens, expon, denom
           ! local variables
       DOUBLE PRECISION :: dtevent(nmetier)     ! time between trawling for this metier
       DOUBLE PRECISION :: tnextmetier(nmetier) ! next trawling time for this metier
       DOUBLE PRECISION :: tnext, tprev, dt, tt     ! time of next trawl
       LOGICAL :: saveEnd 
       DOUBLE PRECISION, PARAMETER :: tiny = 1.d-10
       maxevents = eventnr(nevent)  ! last event is max event nr (i.e. sorted)
        
        ! for every species
        
        DO N = 1, nspec
          
          tprev = 0.d0
          tnext = 1d0/sar(N, 1)
          
          DO I = 1, nmetier
           dtevent(I)       = 1d0/sar(N, I)
           tnextmetier(I)   = 0.5D0 * dtevent(I) 
           
           IF (tnextmetier(I) < tnext) THEN
             tnext = tnextmetier(I)
           END IF
          
          END DO
         
          nev     = 1
          dens    = B0(N)              ! first fishing
 
         ! Go from event to event         
          DO J = 1, maxevents
          
             ! which metier is trawling (can be more than one)
             DO I = 1, nmetier
                
                IF (tnextmetier(I) <= tnext + tiny) THEN
                
                  ! new condition after the event
                  dens  = dens*(1.d0 -d(N, I)) 
                  
                  ! update the metier
                  tnextmetier(I)  = tnextmetier(I) + dtevent(I) 
                  
                END IF  ! tnextmetier == tnext
              
             END DO    ! I
             
             tnext = tnextmetier(1)
             
             DO I = 2, nmetier
               IF (tnextmetier(I) < tnext) THEN
                tnext = tnextmetier(I)
               END IF
             END DO
              
              ! save only if this event is requested
              dt = tnext - tprev

              saveEnd = .FALSE.
              IF (J == eventnr(nev))  THEN 
                B(N, nev) = dens  ! save the density AFTER fishing
                
                Denom = ((K(N)-dens)*Dexp(-r(N)*dt)+dens)/K(N)
                Denom  = DLOG(Denom)
                Bmean(N, nev) =  K(N) + 1.d0/dt*K(N)/r(N)*Denom
                saveEnd = .TRUE.
              END IF 
              
             ! new initial condition before the event(s)
             Denom  = dens + (K(N)-dens)*Dexp(-r(N)*dt)
             dens   = dens * K(N)/Denom   
              
              IF (saveEnd ) THEN
                                 
                Bend(N, nev) = dens
                times(N, nev) = tprev                                 
                tend (N, nev) = tnext                                
                nev = nev + 1
                IF (nev > nevent) EXIT
              END IF
              
             IF (nev > nevent) EXIT
             
             tprev  = tnext

          END DO   ! with J
          
         END DO    ! with N

       END SUBROUTINE metier_event

