!===============================================================================
! Finds the steady-state solution of the perturbation model with events
!===============================================================================

       SUBROUTINE perturb_steady(nspec, sar, K, r, d, atol, rtol, densini,     & 
                               steadytrawl, steadybefore, steadyafter,         &
                               steadymean, steadytimes)
          
       IMPLICIT NONE
       
       
       ! Number of taxa in parameter data
       INTEGER, INTENT(IN) :: nspec
      
       ! taxon parameters 
       DOUBLE PRECISION, INTENT(IN) :: K(nspec)   ! carrying capacity
       DOUBLE PRECISION, INTENT(IN) :: r(nspec)   ! rate of increase
       DOUBLE PRECISION, INTENT(IN) :: d(nspec)   ! depletion fraction
       
       ! fishing parameters 
       DOUBLE PRECISION, INTENT(IN) :: sar(nspec) ! swept area ratio 
       DOUBLE PRECISION, INTENT(IN) :: densini(nspec) ! initial condition
       DOUBLE PRECISION, INTENT(IN) :: atol, rtol ! absolute& relative tolerance

        
       DOUBLE PRECISION, INTENT(OUT) :: steadytimes (nspec)
       DOUBLE PRECISION, INTENT(OUT) :: steadymean  (nspec) 
       DOUBLE PRECISION, INTENT(OUT) :: steadybefore(nspec) 
       DOUBLE PRECISION, INTENT(OUT) :: steadyafter (nspec) 
       INTEGER, INTENT(OUT)          :: steadytrawl (nspec) 
       
       INTEGER :: I, N
       INTEGER, PARAMETER :: maxevents = 1e6
       
       DOUBLE PRECISION :: ttrawl, dens, prevd, delt, denom, p, expon, Di
    
        ! initialisation
        steadytrawl(:) = 0
        steadytimes(:) = 0.d0
        
        ! for every species
        
        DO N = 1, nspec
          
          ttrawl = 0.d0               ! time of trawling
          expon  = Dexp(-r(N)/sar(N))
          p      = (1.d0-d(N))        ! proportion remaining
          dens   = densini(N)*p       ! first fishing
          
         ! Go from event to event         
          DO I = 1, maxevents
              
              prevd  = dens     ! keep previous value of dens     

              ! New density
              denom  = p*dens + (K(N) - p*dens) * expon
              dens   = p*dens*K(N) / denom
              
              ! check convergence
              delt   = DABS(dens - prevd)
              
              IF (dens <= atol .OR. delt < dens*rtol + atol) THEN
                
                ! mean of density before and after fishing
                 steadybefore(N) = prevd
                 steadyafter (N) = prevd*(1.d0-d(N))  
                ! average density in interval 
                 Di = steadyafter(N) 
                 denom = ((K(N)-Di) * expon +Di)/K(N)
                 steadymean(N)  = K(N) + sar(N)*K(N)/r(N)* log(denom)
                 steadytimes(N) = ttrawl
                 steadytrawl(N) = I-1                 
                 EXIT
              END IF
              
              ttrawl = ttrawl + 1d0/sar(N)  ! next trawling event
          
          END DO   ! with I
          
         END DO    ! with N
         
       END SUBROUTINE perturb_steady

!===============================================================================
! analytical solution of the perturbation model with events
! Same sar for all species -> input of event times
!===============================================================================

       SUBROUTINE perturb_times(nspec, nevent, ntimes, B0,                     & 
                                K, r, d, times, events, B, dTrawl)
          
       IMPLICIT NONE
       
       
       ! Number of taxa in parameter data, number of events, length of times
       INTEGER, INTENT(IN) :: nspec, nevent, ntimes
       
       ! taxon parameters 
       DOUBLE PRECISION, INTENT(IN) :: K(nspec)  ! carrying capacity
       DOUBLE PRECISION, INTENT(IN) :: r(nspec)  ! rate of increase
       DOUBLE PRECISION, INTENT(IN) :: d(nspec)  ! depletion fraction
       
       ! output times and trawling
       DOUBLE PRECISION, INTENT(IN) :: times(ntimes)  ! output times
       DOUBLE PRECISION, INTENT(IN) :: events(nevent) ! timing of events
       
       ! initial condition
       DOUBLE PRECISION, INTENT(INOUT) :: B0(nspec)

       ! output varialbes
       DOUBLE PRECISION, INTENT(OUT)   :: B(nspec, ntimes)
       
       ! all the trawling-induced mortalities
       DOUBLE PRECISION, INTENT(OUT)   :: dTrawl(nspec, nevent)

       INTEGER :: I, J, N, nT, nE
       DOUBLE PRECISION :: tnext, tt, t0, denom
    
        nT     = 0              ! index to start of times in loop
        t0     = times(1)       ! time in loop

        nE     = 1              ! index to first event  
        
        IF (events(1) == t0) THEN  ! start with an event => update initial cond.
          nE   = 2
          
          DO N = 1, nspec
            B0(N)        = B0(N)*(1.d0 -d(N))   ! depletion on initial condition
            dTrawl(N, 1) = B0(N)*d(N)
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
              Denom  = B0(N)+(K(N)-B0(N))*Dexp(-r(N)*(tnext-t0))
              B0(N)  = B0(N)*K(N)/Denom
              IF (times(j) > tnext) dTrawl(N, I) = B0(N)*d(N)
          
              B0(N)  = B0(N)*(1.d0 -d(N))
          END DO
          t0     = tnext
        
        END DO  ! with I
       END SUBROUTINE perturb_times
       
!===============================================================================
! analytical solution of the perturbation model with events
! Input of one sar per species 
!===============================================================================

       SUBROUTINE perturb_times2(nspec, ntimes, B0, K, r, d, sar,               & 
                                 times, tstartperturb, tendperturb, B)
          
       IMPLICIT NONE
       
       
       ! Number of taxa in parameter data, number of events, length of times
       INTEGER, INTENT(IN) :: nspec,  ntimes
       
       DOUBLE PRECISION, INTENT(IN) :: K(nspec)    ! carrying capacity
       DOUBLE PRECISION, INTENT(IN) :: r(nspec)    ! rate of increase
       DOUBLE PRECISION, INTENT(IN) :: d(nspec)    ! depletion fraction
       DOUBLE PRECISION, INTENT(IN) :: sar(nspec)  ! swept area ratio
       DOUBLE PRECISION, INTENT(IN) :: tstartperturb ! first trawling
       DOUBLE PRECISION, INTENT(IN) :: tendperturb ! last trawling
      
       DOUBLE PRECISION, INTENT(IN) :: times(ntimes)  ! output times

       DOUBLE PRECISION, INTENT(INOUT) :: B0(nspec)

       DOUBLE PRECISION, INTENT(OUT)   :: B(nspec, ntimes)
       
       INTEGER :: J, N
       DOUBLE PRECISION :: dtevent, tnext, tprev, denom
    
         
        DO N = 1, nspec
        
        ! start with an event => update initial cond.
        
         tprev   = times(1)
         
         IF (sar(N) > 0.d0) THEN 
            dtevent = 1.d0 / sar(N)           ! time between events 
            
            IF (tstartperturb == times(1)) THEN
              tnext   = times(1)  + dtevent     ! next event
              B0(N)   = B0(N)*(1.d0-d(N))       ! reduce density
            
            ELSE
              tnext   = tstartperturb
            ENDIF  
            
         ELSE 
         
            dtevent = times(ntimes) + 1.d0    ! fictively after last times
            tnext  = dtevent
         
         END IF
         
         ! loop over all times except the first
         B(N, 1) = B0(N)

         IF (tnext > tendperturb) tnext = times(ntimes) + 1.d0 
         
         DO J = 2, ntimes
           
           DO WHILE (times(J) > tnext)    ! take care of the event(s) first
           
              tprev  = tnext
              
              ! new condition before the event
              
              IF (tnext > tstartperturb) THEN
                Denom  = B0(N)+(K(N)-B0(N))*Dexp(-r(N)*dtevent)
                B0(N)  = B0(N)*K(N)/Denom   
              ELSE
                Denom  = B0(N)+(K(N)-B0(N))*Dexp(-r(N)*(tstartperturb - times(1)))
                B0(N)  = B0(N)*K(N)/Denom   
              END IF              

              ! new condition after the event
              B0(N)  = B0(N)*(1.d0 -d(N)) 
              tnext  = tnext + dtevent 
              
              IF (tnext > tendperturb) tnext = times(ntimes) + 1.d0 
           END DO

           ! Estimate density at times(j)
            Denom  = B0(N)+(K(N)-B0(N))*Dexp(-r(N)*(times(j)-tprev))
            B(N,J) = B0(N)*K(N)/Denom

          END DO   ! with J
        END DO  ! with N
       
       END SUBROUTINE perturb_times2

