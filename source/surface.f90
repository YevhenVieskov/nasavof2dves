!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%% INIT_PARTICLES %%%%%%%%%%%%%%%%%%%%%%%%%%% 
!
!   SUBROUTINE INIT_PARTICLES (N, imax, jmax, delx, dely, ppc, problem, &
 !                             U, V, Partlines)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!
! Initialize particles for free boundary problems to define fluid domain 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!
!      USE nrtype
!      IMPLICIT NONE

!      INTEGER, INTENT(OUT) :: N
!      INTEGER, INTENT(IN) :: imax, jmax, ppc
!      REAL(RP), INTENT(IN) :: delx, dely
!      REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: U, V
!      TYPE(PARTICLELINE), DIMENSION(:), INTENT(INOUT) :: Partlines
!      CHARACTER(LEN=30), INTENT(IN) :: problem
!
! local variables
!
!      INCLUDE 'interfaces.h'
!      INTEGER :: i, j, ip, jp
!      REAL(RP) :: x, y
!      REAL(RP) :: height, rad, mpx, mpy, vstart
!
! Initialization of some parameters 
!
!      IF (problem == "dam") THEN
!         N = 1
!      ELSE IF (problem == "drop") THEN
!         N = 2
!         height = 1./2.*jmax*dely      ! Height of the basin   
!         rad    = 0.1*jmax*dely        ! Radius of the drop    
!         mpx    = 0.5*imax*delx        ! Mid point of the drop 
!         mpy    = 2./3.*jmax*dely
!         vstart = -2.0                 ! Initial velocity of the drop 
!      END IF

!      DO i = 1,N
!         Partlines(i)%length = 0
!         Partlines(i)%Particles%x = -1.
!         Partlines(i)%Particles%y = -1.
!      END DO
!
! Set the particles 
!
!      DO i = 1,imax
!         DO j = 1,jmax

!            DO  ip = 1,ppc
!               x = (i-1)*delx+(ip-.5)/REAL(ppc)*delx
!               DO  jp = 1,ppc
!                  y = (j-1)*dely+(jp-.5)/REAL(ppc)*dely
    
!                  IF (problem == "dam") THEN
!                     IF (x < 0.2*imax*delx) THEN
!                        CALL SET_PART (Partlines(1), x, y)
!                     END IF
!                  ELSE IF (problem == "drop") THEN
!                     IF (y < height) THEN
!                        CALL SET_PART (Partlines(1), x, y)
!                     ELSE IF ((x-mpx)**2 + (y-mpy)**2 <= rad**2) THEN
!                        CALL SET_PART (Partlines(2), x, y)
!                        V(i,j) = vstart
!                     END IF
!                  END IF
!
!               END DO
!            END DO
!
!         END DO
!      END DO

!      RETURN

!   END SUBROUTINE INIT_PARTICLES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET_PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   SUBROUTINE SET_PART (Partline, x, y)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!
! Add particle to "Partline" at (x,y)                              
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
!      USE nrtype
!      IMPLICIT NONE
 
!      REAL(RP), INTENT(IN) :: x, y
!      TYPE(PARTICLELINE), INTENT(INOUT), TARGET :: Partline
!
!  local variables
!
!      TYPE(PARTICLE), POINTER :: part

!      ALLOCATE(part)    ! 
!      part%x = x        ! create particle at (x,y)
!      part%y = y        ! 

!      part%next => Partline%Particles%next   ! add it to "Partline"  
!      Partline%Particles%next => part        ! in the first position  
!      Partline%length = Partline%length + 1  ! after the dummy

!   END SUBROUTINE SET_PART
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MARK_CELLS %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!
!   SUBROUTINE MARK_CELLS (FLAG, imax, jmax, delx, dely, ifull, isurf, &
!                          N, Partlines)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!
! Mark the cells of the fluid domain                            
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!
!      USE nrtype
!      IMPLICIT NONE
!      INTEGER, INTENT(IN) :: imax, jmax
!      INTEGER, INTENT(OUT) :: ifull, isurf, N
!      REAL(RP), INTENT(IN) :: delx, dely
!      INTEGER(I2B), DIMENSION(0:,0:), INTENT(OUT) :: FLAG
!      TYPE(PARTICLELINE),DIMENSION(:),INTENT(INOUT),TARGET :: Partlines
!
!  local variables
!
!      INCLUDE 'defs.h'
!      INTEGER :: i, j, k
 !     REAL(RP) :: x, y
 !     TYPE(PARTICLE), POINTER :: part, temp, help
!
!  Set all cells which are not obstacle cells to empty cells 
!
!      DO i = 0,imax+1
!         DO j = 0,jmax+1
!            IF (FLAG(i,j) >= C_F) THEN
!               FLAG(i,j) = IAND(IOR(FLAG(i,j),C_E),NOT(C_NSWO))
!            END IF
!         END DO
!      END DO
!
!  Mark cells containing particles as fluid cells (loop over particles) 
!
!      DO k = 1,N

!         part => Partlines(k)%Particles

!         DO

!            IF (.NOT.(ASSOCIATED(part%next))) EXIT

!            temp => part%next

!            x = temp%x
!            y = temp%y

!            i = INT(x/delx) + 1;
!            j = INT(y/dely) + 1;

!            IF (FLAG(i,j) < C_F) THEN !!!!!!!error
               !
               ! delete particles that have moved into obstacle cells 
               ! NOTE:  predecessor pointer doesn't advance if a
               !        deletion is involved
               !
!               help => temp%next
!               DEALLOCATE(part%next)
!               part%next => help
!               Partlines(k)%length = Partlines(k)%length - 1
!            ELSE
!               FLAG(i,j) = IAND(FLAG(i,j),NOT(C_E))
!            END IF

 !           part => part%next
           
!         END DO
!      END DO
!
! Mark surface cells 
!
!      ifull = 0
!      isurf = 0

!      DO j = 1,jmax
!         DO i = 1,imax

!            IF ( (IAND(FLAG(i,j),C_F) /= 0) .AND. (FLAG(i,j) < C_E) ) THEN
!         
!               IF (IAND(FLAG(i-1,j),C_E) /= 0) THEN
!                  FLAG(i,j) = IOR(FLAG(i,j),C_W)
!               END IF

!               IF (IAND(FLAG(i+1,j),C_E) /= 0) THEN
!                  FLAG(i,j) = IOR(FLAG(i,j),C_O)
!               END IF

!               IF (IAND(FLAG(i,j-1),C_E) /= 0) THEN
!                  FLAG(i,j) = IOR(FLAG(i,j),C_S)
!               END IF

!               IF (IAND(FLAG(i,j+1),C_E) /= 0) THEN
!                  FLAG(i,j) = IOR(FLAG(i,j),C_N)
!               END IF

!               IF (FLAG(i,j) < C_O) THEN  
!                  ifull = ifull + 1
!               ELSE
!                  isurf = isurf + 1
!               END IF

!            END IF
      
!         END DO
!      END DO
!
!  DIAGNOSTIC:  Output geometry of the fluid domain
!
!     WRITE (6,*) ' '
!     WRITE (6,*) ' Geometry of the fluid domain'
!     WRITE (6,*) ' '
!     DO j = jmax+1,0,-1
!        WRITE (6,'(1x,200i5)') (FLAG(i,j), i=0,imax+1)
!     END DO

!      RETURN

!   END SUBROUTINE MARK_CELLS
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET_UVP_SURFACE %%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE SET_UVP_SURFACE (U, V, P, FLAG, GX, GY, imax, jmax, &
                               Re, delx, dely, delt)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Set boundary values at free surface                     
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imax, jmax
      REAL(RP), INTENT(IN) :: GX, GY, Re, delx, dely, delt
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, V, P
!
! local variables
!
      INCLUDE 'defs.h'
      INTEGER :: i, j
!
! Set velocity values in empty cells to zero 
!
      DO j = 1,jmax
         DO i = 1,imax-1
            IF ( (IAND(FLAG(i,j),  C_E) /= 0) .AND. &
                 (IAND(FLAG(i+1,j),C_E) /= 0) ) THEN
               U(i,j) = 0.0
            END IF
         END DO
      END DO

      DO  j = 1,jmax-1 
         DO i = 1,imax
            IF ( (IAND(FLAG(i,j),  C_E) /= 0) .AND. &
                 (IAND(FLAG(i,j+1),C_E) /= 0)) THEN
               V(i,j) = 0.0
            END IF
         END DO
      END DO

      DO j = 1,jmax
         DO i = 1,imax
          !
          ! treat only surface cells 
          !
          IF ( (IAND(FLAG(i,j),C_E) == 0) .OR. (FLAG(i,j) < C_O) ) THEN

           SELECT CASE (IAND(FLAG(i,j),C_NSWO)) ! mask NSWO_E=0x0f00    
	                                        ! filters surface cells 
             CASE (C_N)   
               V(i,j) = V(i,j-1)-dely/delx*(U(i,j)-U(i-1,j))
               IF (IAND(FLAG(i-1,j+1),C_E) /= 0) THEN
                 U(i-1,j+1) = U(i-1,j)-dely/delx*(V(i,j)-V(i-1,j))
               END IF

             CASE (C_S)  
               V(i,j-1) = V(i,j)+dely/delx*(U(i,j)-U(i-1,j))
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN 
                 U(i-1,j-1) = U(i-1,j)+dely/delx*(V(i,j-1)-V(i-1,j-1))
               END IF
 
	     CASE (C_O)
               U(i,j) = U(i-1,j)-delx/dely*(V(i,j)-V(i,j-1))
               IF (IAND(FLAG(i+1,j-1),C_E) /= 0) THEN 
                 V(i+1,j-1) = V(i,j-1)-delx/dely*(U(i,j)-U(i,j-1))
               END IF

	     CASE (C_W) 
               U(i-1,j) = U(i,j)+delx/dely*(V(i,j)-V(i,j-1))
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN 
                 V(i-1,j-1) = V(i,j-1)+delx/dely*(U(i-1,j)-U(i-1,j-1))
               END IF

             CASE (C_NO) 
               U(i,j) = U(i-1,j)
               V(i,j) = V(i,j-1)
               IF (IAND(FLAG(i-1,j+1),C_E) /= 0) THEN
                 U(i-1,j+1) = U(i-1,j)-dely/delx*(V(i,j)-V(i-1,j))
               END IF
               IF (IAND(FLAG(i+1,j+1),C_E) /= 0) THEN
                 U(i,j+1) = U(i,j)
                 V(i+1,j) = V(i,j)
               END IF
               IF (IAND(FLAG(i+1,j-1),C_E) /= 0) THEN
                 V(i+1,j-1) = V(i,j-1)-delx/dely*(U(i,j)-U(i,j-1))
               END IF

             CASE (C_NW)
               U(i-1,j) = U(i,j)
               V(i,j)   = V(i,j-1)
               IF (IAND(FLAG(i-1,j+1),C_E) /= 0) THEN
                 U(i-1,j+1) = U(i-1,j)
                 V(i-1,j)   = V(i,j)
               END IF
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 V(i-1,j-1) = V(i,j-1)+delx/dely*(U(i-1,j)-U(i-1,j-1))
               END IF

             CASE (C_SW)
               U(i-1,j) = U(i,j)
               V(i,j-1) = V(i,j)
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 U(i-1,j-1) = U(i-1,j)
                 V(i-1,j-1) = V(i,j-1)
               END IF

             CASE (C_SO)
               U(i,j)   = U(i-1,j)
               V(i,j-1) = V(i,j)
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 U(i-1,j-1) = U(i-1,j)+dely/delx*(V(i,j-1)-V(i-1,j-1))
               END IF
               IF (IAND(FLAG(i+1,j-1),C_E) /= 0) THEN
                 U(i,j-1)   = U(i,j)
                 V(i+1,j-1) = V(i,j-1)
               END IF

             CASE (C_WO)
               U(i,j)   = U(i,j)   + delt*GX
               U(i-1,j) = U(i-1,j) + delt*GX
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 V(i-1,j-1) = V(i,j-1)+delx/dely*(U(i-1,j)-U(i-1,j-1))
               END IF
               IF (IAND(FLAG(i+1,j-1),C_E) /= 0) THEN
                 V(i+1,j-1) = V(i,j-1)-delx/dely*(U(i,j)-U(i,j-1))
               END IF

             CASE (C_NS)
               V(i,j)   = V(i,j)   + delt*GY
               V(i,j-1) = V(i,j-1) + delt*GY
               IF (IAND(FLAG(i-1,j+1),C_E) /= 0) THEN
                 U(i-1,j+1) = U(i-1,j)-dely/delx*(V(i,j)-V(i-1,j))
               END IF
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 U(i-1,j-1) = U(i-1,j)+dely/delx*(V(i,j-1)-V(i-1,j-1))
               END IF

             CASE (C_NWO)
               V(i,j) = V(i,j-1)-dely/delx*(U(i,j)-U(i-1,j))
               U(i,j) = U(i,j) + delt*GX
               U(i-1,j) = U(i-1,j) + delt*GX
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 V(i-1,j-1) = V(i,j-1)+delx/dely*(U(i-1,j)-U(i-1,j-1))
               END IF
               IF (IAND(FLAG(i+1,j-1),C_E) /= 0) THEN
                 V(i+1,j-1) = V(i,j-1)-delx/dely*(U(i,j)-U(i,j-1))
               END IF
               IF (IAND(FLAG(i-1,j+1),C_E) /= 0) THEN
                 V(i-1,j)   = V(i,j)
                 U(i-1,j+1) = U(i-1,j)
               END IF
               IF (IAND(FLAG(i+1,j+1),C_E) /= 0) THEN
                 V(i+1,j) = V(i,j)
                 U(i,j+1) = U(i,j)    
               END IF

             CASE (C_NSW)
               U(i-1,j) = U(i,j)+delx/dely*(V(i,j)-V(i,j-1))
               V(i,j)   = V(i,j) + delt*GY
               V(i,j-1) = V(i,j-1) + delt*GY
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 V(i-1,j-1)  = V(i,j-1)
                 U(i-1,j-1)  = U(i-1,j)
               END IF
               IF (IAND(FLAG(i-1,j+1),C_E) /= 0) THEN
                 V(i-1,j)   = V(i,j)
                 U(i-1,j+1) = U(i-1,j)
               END IF

             CASE (C_SWO)
               V(i,j-1) = V(i,j)+dely/delx*(U(i,j)-U(i-1,j))
               U(i,j)   = U(i,j) + delt*GX
               U(i-1,j) = U(i-1,j) + delt*GX
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 U(i-1,j-1) = U(i-1,j)
                 V(i-1,j-1) = V(i,j-1)
               END IF
               IF (IAND(FLAG(i+1,j-1),C_E) /= 0) THEN
                 U(i,j-1)    = U(i,j)
                 V(i+1,j-1)  = V(i,j-1)
               END IF

             CASE (C_NSO)
               U(i,j)   = U(i-1,j)-delx/dely*(V(i,j)-V(i,j-1))
               V(i,j)   = V(i,j) + delt*GY
               V(i-1,j) = V(i-1,j) + delt*GY
               IF (IAND(FLAG(i-1,j+1),C_E) /= 0) THEN
                 U(i-1,j+1) = U(i-1,j)-dely/delx*(V(i,j)-V(i-1,j))
               END IF
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 U(i-1,j-1) = U(i-1,j)+dely/delx*(V(i,j-1)-V(i-1,j-1))
               END IF
               IF (IAND(FLAG(i+1,j-1),C_E) /= 0) THEN
                 U(i,j-1)   = U(i,j)
                 V(i+1,j-1) = V(i,j-1)
               END IF
               IF (IAND(FLAG(i+1,j+1),C_E) /= 0) THEN
                 U(i,j+1)    = U(i,j)
                 V(i+1,j)    = V(i,j)
               END IF

             CASE (C_NSWO)
               U(i,j)   = U(i,j)   + delt*GX
               U(i-1,j) = U(i-1,j) + delt*GX
               V(i,j)   = V(i,j)   + delt*GY
               V(i,j-1) = V(i,j-1) + delt*GY
               IF (IAND(FLAG(i-1,j+1),C_E) /= 0) THEN
                 U(i-1,j+1) = U(i-1,j)
                 V(i-1,j)   = V(i,j)
               END IF
               IF (IAND(FLAG(i+1,j+1),C_E) /= 0) THEN
                 U(i,j+1) = U(i,j)
                 V(i+1,j) = V(i,j)
               END IF
               IF (IAND(FLAG(i-1,j-1),C_E) /= 0) THEN
                 U(i-1,j-1)  = U(i-1,j)
                 V(i-1,j-1)  = V(i,j-1)
               END IF
               IF (IAND(FLAG(i+1,j-1),C_E) /= 0) THEN
                 U(i,j-1)    = U(i,j)
                 V(i+1,j-1)  = V(i,j-1)
               END IF

	     CASE DEFAULT

           END SELECT

          END IF

         END DO
      END DO
! 
! Second loop DO  pressure boundary values 
!

      DO j = 1,jmax 
         DO i = 1,imax
       
            IF (.NOT.((IAND(FLAG(i,j),C_E)/=0).OR.(FLAG(i,j)<C_O))) THEN

              SELECT CASE (IAND(FLAG(i,j),C_NSWO))  

	        CASE (C_N)
                  P(i,j) = 2./Re/dely*(V(i,j)-V(i,j-1))
	        CASE (C_S)
                  P(i,j) = 2./Re/dely*(V(i,j)-V(i,j-1))  
	        CASE (C_O)
                  P(i,j) = 2./Re/delx*(U(i,j)-U(i-1,j)) 
	        CASE (C_W)
                  P(i,j) = 2./Re/delx*(U(i,j)-U(i-1,j)) 
                CASE (C_NO)
                  P(i,j) = 1./Re/2.* &
                     ( (U(i,j)+U(i-1,j)-U(i,j-1)-U(i-1,j-1))/dely &
                   +   (V(i,j)+V(i,j-1)-V(i-1,j)-V(i-1,j-1))/delx )
                CASE (C_NW)
                  P(i,j) = -1./Re/2.* &
                     ( (U(i,j)+U(i-1,j)-U(i,j-1)-U(i-1,j-1))/dely &
                   +   (V(i+1,j)+V(i+1,j-1)-V(i,j)-V(i,j-1))/delx)
                CASE (C_SW)
                  P(i,j) = 1./Re/2.* &
                     ( (U(i,j+1)+U(i-1,j+1)-U(i,j)-U(i-1,j))/dely &
                   +   (V(i+1,j)+V(i+1,j-1)-V(i,j)-V(i,j-1))/delx)
                CASE (C_SO)
                  P(i,j) = -1./Re/2.* &
                     ( (U(i,j+1)+U(i-1,j+1)-U(i,j)-U(i-1,j))/dely &
                   +   (V(i,j)+V(i,j-1)-V(i-1,j)-V(i-1,j-1))/delx)
                CASE (C_WO)
                  P(i,j) = 0.
                CASE (C_NS)
                  P(i,j) = 0.
                CASE (C_NWO)
                  P(i,j) = 0. 
                CASE (C_NSW)
                  P(i,j) = 0. 
                CASE (C_SWO)
                  P(i,j) = 0. 
                CASE (C_NSO)
                  P(i,j) = 0. 
                CASE (C_NSWO)
                  P(i,j) = 0. 
                CASE DEFAULT

              END SELECT
                P(i,j)=0.0
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

   END SUBROUTINE SET_UVP_SURFACE
