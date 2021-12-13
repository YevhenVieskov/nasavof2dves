!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETBCOND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE SETBCOND (U, V, P, TEMP, FLAG, imax, jmax, wW, wE, wN, wS,FS)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Setting the boundary conditions at the boundary strip.          
!  The flags wW,wE,wN, and wS can have the values:                 
!
!  1 = slip               2 = no-slip                             
!  3 = outflow            4 = periodic                            
!
!  Moreover, no-slip conditions are set at internal obstacle cells
!  by default.                                                    
!
!  For temperature, adiabatic boundary conditions are set.        
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: imax, jmax, wW, wE, wN, wS
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, V, P, TEMP
	  REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: FS
!
!  local variables
!
      INCLUDE 'defs.h'

      INTEGER i, j

      DO j = 0,jmax+1 
         !
         ! western and eastern boundary
         !
         IF ( wW == 1 ) THEN        ! free-slip
            U(0,j) = 0.0              ! u = 0
            V(0,j) = V(1,j)           ! dv/dn = 0
			FS(0,J)=FS(1,J)

         ELSE IF (wW == 2 ) THEN    ! no-slip
            U(0,j) = 0.0              ! u = 0  
            V(0,j) = (-1.0)*V(1,j)    ! v=0 at the boundary by averaging
			FS(0,J)=FS(1,J)
         ELSE IF (wW == 3) THEN     ! outflow
            U(0,j) = U(1,j) 
            V(0,j) = V(1,j) 
			FS(0,J)=FS(1,J)
         ELSE IF (wW == 4 ) THEN    ! periodic
            U(0,j) = U(imax-1,j) 
            V(0,j) = V(imax-1,j)      ! left and right cells
            V(1,j) = V(imax,j)        ! are overlapping
            P(1,j) = P(imax,j) 
			FS(1,J)=FS(imax,J) 
         END IF

         TEMP(0,j) = TEMP(1,j)      ! dT/dn = 0 

         IF (wE == 1 ) THEN         ! free-slip
            U(imax,j) = 0.0          
            V(imax+1,j) = V(imax,j) 
			FS(imax+1,J)=FS(imax,J)  
         ELSE IF (wE == 2 ) THEN    ! no-slip
            U(imax,j) = 0.0
            V(imax+1,j) = (-1.0)*V(imax,j)
			FS(imax+1,J)=FS(imax,J) 
         ELSE IF (wE == 3) THEN     ! outflow
            U(imax,j) = U(imax-1,j)
            V(imax+1,j) = V(imax,j)
			FS(imax+1,J)=FS(imax,J) 
         ELSE IF (wE == 4 ) THEN    ! periodic
            U(imax,j) = U(1,j)
            V(imax+1,j) = V(2,j)
			FS(imax+1,J)=FS(1,J) 
         END IF

         TEMP(imax+1,j) = TEMP(imax,j)

      END DO

      DO i=0,imax+1
         !
         ! northern and southern boundary
         !
         IF (wN == 1 ) THEN
            V(i,jmax) = 0.0
            U(i,jmax+1) = U(i,jmax)
			FS(I,jmax+1)=FS(I,jmax)
         ELSE IF (wN == 2 ) THEN 
            V(i,jmax) = 0.0
            U(i,jmax+1) = (-1.0)*U(i,jmax)
			FS(I,jmax+1)=FS(I,jmax)
         ELSE IF (wN == 3) THEN
            V(i,jmax) = V(i,jmax-1)
            U(i,jmax+1) = U(i,jmax)
			FS(I,jmax+1)=FS(I,jmax)
         ELSE IF (wN == 4 ) THEN
            V(i,jmax) = V(i,1)
            U(i,jmax+1) = U(i,2)
			FS(I,jmax)=FS(I,1)
         END IF

         TEMP(i,0) = TEMP(i,1) 

         IF (wS == 1 ) THEN
            V(i,0) = 0.0
            U(i,0) = U(i,1)
			FS(I,0)=FS(I,1)
         ELSE IF (wS == 2 ) THEN
            V(i,0) = 0.0
            U(i,0) = (-1.0)*U(i,1)
			FS(I,0)=FS(I,1)
         ELSE IF (wS == 3) THEN
            V(i,0) = V(i,1)
            U(i,0) = U(i,1)
			FS(I,0)=FS(I,1)
         ELSE IF (wS == 4 ) THEN
            V(i,0) = V(i,jmax-1)
            U(i,0) = U(i,jmax-1)
            U(i,1) = U(i,jmax)
            P(i,1) = P(i,jmax)
			FS(I,1)=FS(I,jmax)
         END IF

         TEMP(i,jmax+1) = TEMP(i,jmax) 

      END DO
!
! set the boundary values at inner obstacle cells 
!                  (only no-slip) 
!
      DO i = 1,imax
         DO j = 1,jmax

            IF (IAND(FLAG(i,j),C_X) /= 0) THEN ! mask C_X=000f filters the
                                               ! obstacle cells adjacent to
                                               ! fluid cells 
             SELECT CASE (FLAG(i,j))

               CASE (B_N)
		       V(i,j)   = 0.0
                       U(i,j)   = -U(i,j+1)
                       U(i-1,j) = -U(i-1,j+1)
                       TEMP(i,j) = TEMP(i,j+1)
               CASE (B_E)
		       U(i,j)   = 0.0
                       V(i,j)   = -V(i+1,j)
                       V(i,j-1) = -V(i+1,j-1)
                       TEMP(i,j) = TEMP(i+1,j)
               CASE (B_S)
		       V(i,j-1) = 0.0
                       U(i,j)   = -U(i,j-1)
                       U(i-1,j) = -U(i-1,j-1)
                       TEMP(i,j) = TEMP(i,j-1)
               CASE (B_W)
		       U(i-1,j) = 0.0
                       V(i,j)   = -V(i-1,j)
                       V(i,j-1) = -V(i-1,j-1)
                       TEMP(i,j) = TEMP(i-1,j)
               CASE (B_NE)
		       V(i,j)   = 0.0
                       U(i,j)   = 0.0
                       V(i,j-1) = -V(i+1,j-1)
                       U(i-1,j) = -U(i-1,j+1)
                       TEMP(i,j) = 0.5*(TEMP(i,j+1)+TEMP(i+1,j))
               CASE (B_SE)
		       V(i,j-1) = 0.0
                       U(i,j)   = 0.0
                       V(i,j)   = -V(i+1,j)
                       U(i-1,j) = -U(i-1,j-1)
                       TEMP(i,j) = 0.5*(TEMP(i,j-1)+TEMP(i+1,j))
               CASE (B_SW)
		       V(i,j-1) = 0.0
                       U(i-1,j) = 0.0
                       V(i,j)   = -V(i-1,j)
                       U(i,j)   = -U(i,j-1)
                       TEMP(i,j) = 0.5*(TEMP(i,j-1)+TEMP(i-1,j))
               CASE (B_NW)
		       V(i,j)   = 0.0
                       U(i-1,j) = 0.0
                       V(i,j-1) = -V(i-1,j-1)
                       U(i,j)   = -U(i,j+1)
                       TEMP(i,j) = 0.5*(TEMP(i,j+1)+TEMP(i-1,j))
	       CASE DEFAULT

             END SELECT

            END IF

         END DO
      END DO

 DO j = 1,jmax 
  DO i = 1,imax
  IF(FLAG(I,J)>C_B.AND.FLAG(I,J)>C_E) THEN
    DIV(I,J)=(U(I,J)-U(I-1,J))/DX+(V(I,J)-V(I,J-1))/DY
  ELSE
    DIV(I,J)=0.0
  END IF
  END DO
 END DO
	  
	  RETURN

   END SUBROUTINE SETBCOND
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%% SETSPECBCOND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE SETSPECBCOND (problem, U, V, TEMP, imax, jmax, UI, VI)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Setting specific boundary conditions, depending on "problem" 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype

      IMPLICIT NONE

      CHARACTER (LEN=30) :: problem
      REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, V, TEMP
      INTEGER, INTENT(IN) :: imax, jmax
      REAL(RP), INTENT(IN) :: UI, VI
!
!  local variables
!
      INTEGER :: i,j

      IF ( (problem=="drop") .OR. (problem=="dam") ) THEN

         RETURN

!-----------------------------------------------------------
! Driven Cavity: U = 1.0 at the upper boundary              
!-----------------------------------------------------------
  
      ELSE IF (problem=="dcavity") THEN

         DO i = 0,imax
            U(i,jmax+1) = 2.0 - U(i,jmax)   
         END DO
         RETURN

!-----------------------------------------------------------------
! Flow past a backward facing step, with or without free boundary 
!                  U = 1.0 at the left boundary                   
!-----------------------------------------------------------------
      
      ELSE IF ( (problem=="backstep") .OR. (problem=="wave") ) THEN
       
         DO j = jmax/2+1,jmax
            U(0,j) = 1.0
         END DO
         RETURN

!--------------------------------------------------------------
! Flow past an obstacle: U = 1.0 at left boundary              
!--------------------------------------------------------------
        
      ELSE IF ( (problem=="plate") .OR. (problem=="circle")) THEN

         V(0,0) = 2*VI-V(1,0)
         DO j = 1,jmax
            U(0,j) = UI
            V(0,j) = 2*VI-V(1,j)
         END DO
         RETURN

!---------------------------------------------------------------------
! Inflow for injection molding: U = 1.0 in the mid of left boundary   
!---------------------------------------------------------------------

      ELSE IF (problem=="molding") THEN

         DO j = INT(0.4*jmax) + 1, INT(0.6*jmax)
            U(0,j) = 1.0       
         END DO
         RETURN

!------------------------------------------------------------------
! natural convection or fluidtrap: left T = 0.5 right T = -0.5     
!                          upper and lower wall adiabatic          
!------------------------------------------------------------------

      ELSE IF ((problem=="convection") .OR. (problem=="fluidtrap")) THEN

         DO j = 0,jmax+1
            TEMP(0,j) = 2*(0.5)-TEMP(1,j)           ! left wall heated  
            TEMP(imax+1,j) = 2*(-0.5)-TEMP(imax,j)  ! right wall cooled 
         END DO
       
         DO i = 0,imax+1
	    TEMP(i,0) = TEMP(i,1)
	    TEMP(i,jmax+1) = TEMP(i,jmax)           ! adiabatic walls 
         END DO
         RETURN

!----------------------------------------------------
! Rayleigh-Benard flow: top T = -0.5 bottom T = 0.5  
!                       left and right adiabatic     
!----------------------------------------------------

     ELSE IF (problem=="rayleigh") THEN

        DO j = 0,jmax+1
           TEMP(0,j) = TEMP(1,j)         
           TEMP(imax+1,j) = TEMP(imax,j)       !adiabatic walls 
        END DO
  
        DO i = 0,imax+1
           TEMP(i,0) = 2*(0.5)-TEMP(i,1)           !lower wall heated 
           TEMP(i,jmax+1) = 2*(-0.5)-TEMP(i,jmax)  !upper wall cooled 
        END DO

        RETURN

     ELSE

        WRITE (6,*) ' problem is not defined'
        stop ' setspecbcond'

     END IF

   END SUBROUTINE SETSPECBCOND
