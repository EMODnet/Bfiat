
!===============================================================================
!===============================================================================
! Calculates the solution at each event
!===============================================================================
!===============================================================================

       SUBROUTINE eventdensity(nspec, nevent, B0, sar,                         & 
                               K, r, d, eventnr, B)
          
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
                B(N, nev) = p*dens  ! save the density AFTER fishing
                nev = nev + 1
                IF (nev > nevent) EXIT
              END DO 
              
              IF (nev > nevent) EXIT

              ! New density
              Denom  = p*dens + (K(N) - p*dens) * expon
              dens   = p*dens*K(N) / Denom
          END DO   ! with I
          
         END DO    ! with N

       END SUBROUTINE eventdensity

!===============================================================================
! Calculates the solution, one event per species
!===============================================================================

       SUBROUTINE eventdensity2(nspec, B0, sar,                                & 
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

       END SUBROUTINE eventdensity2

